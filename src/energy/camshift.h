// Copyright (C) 2011 by Anders Steen Christensen, Wouter Boomsma, Simon Olsson, Lars Bratholm
//
// This file is part of PHAISTOS
//
// PHAISTOS is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// PHAISTOS is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with PHAISTOS.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef CAMSHIFT
#define CAMSHIFT

#include "energy/energy_term.h"
#include <boost/type_traits/is_base_of.hpp>
#include <iostream>
#include <vector>
#include <math.h>
#include <cmath>
#include <string.h>

#include <boost/random.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/uniform_smallint.hpp>

#include <boost/random/lognormal_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/math/distributions/bernoulli.hpp>
#include "protein/iterators/pair_iterator_chaintree.h" 
#include "protein/chain_fb.h"
#include "protein/definitions.h"
#include "energy/energy.h"
#include "energy/camshift.h"
#include "camshift_data.h"
#include "camshift_predictor.h"
#include <time.h>

namespace phaistos {

//! Class from which Camshift energy terms inherit their predictions
template <typename DERIVED_CLASS>
class TermCamshiftBase: public EnergyTermCommon<DERIVED_CLASS, ChainFB>, 
                        public CamshiftBackendBase {

     //! Class from which Camshift energy terms inherit their predictions
     typedef phaistos::EnergyTermCommon<DERIVED_CLASS, ChainFB> EnergyTermCommon;     

public:

     //! Class variable that holds a pointer to the global random number engine
     RandomNumberEngine *random_number_engine;

     //! Keeps track of whether nuisance parameters have been updated
     bool none_move;

     //! The values of sigma for each atom type in the current MC step
     std::vector<double> weights;

     //! The values of sigma for each atom type in the previous MC step
     std::vector<double> weights_previous;

     //! Counter that prints out 
     unsigned int n_print;

     //! A table containing the chemical shift from the user-defined input file
     std::vector< std::vector<double> > chemshifts_from_file;

     //! A table that contains all chemical shifts predicted from the current structure 
     std::vector< std::vector<double> > protein_predicted_cs;

     //! A table that contains all chemical shifts predicted from the previous structure 
     std::vector< std::vector<double> > protein_predicted_cs_previous;

     //! Local copy of thread id.
     unsigned int thread_id;

     //! Returns samples according to the Normal distribution.·
     //! \param m Mean value.
     //! \param std standard deviation
     //! \param rne A random number engine object
     //! \return A double drawn from the defined normal distribution
     double rnom(const double m, const double std, RandomNumberEngine *rne) {
          boost::normal_distribution<> nom(m, std);
          boost::variate_generator<RandomNumberEngine&, boost::normal_distribution<> > norml( *rne, nom );
          return norml();
     }


     //! Returns a random integer from {min, min+1, ... max-1, max}.
     //! \param min Lowest possible value in the range 
     //! \param max highest possible value in the range
     //! \param rne A random number engine object
     //! \return An integer drawn uniformly from the define range.
     int rand_int(const int min, const int max, RandomNumberEngine *rne) {
          boost::uniform_smallint<> distribution(min, max);
          boost::variate_generator<RandomNumberEngine&, boost::uniform_smallint<> > generator(*rne, distribution);
          return generator();
     }


     //! Update weights (from an unbiased distribution)
     //! \param weights_this The weights vector to be updated.
     //! \param rne A random number engine object
     //! \return a new weight vector.
     std::vector<double> update_weights(const std::vector<double> &weights_this, 
                                        RandomNumberEngine *rne) {

          std::vector<double> weights_new = weights_this;

          const double sigma = 0.05;
          const double mu = 0.0;
          const double increment = rnom(mu, sigma, rne); 

          const int i = rand_int(0, 5, rne);

          weights_new[i] += increment;

          if (weights_new[i] < 0.0) weights_new[i] += std::fabs(2.0 * increment);

          return weights_new;
     }


     void print_weights() {

          std::vector<double> chi_sq = vector_utils::make_vector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
          std::vector<int> n_bins = vector_utils::make_vector(0, 0, 0, 0, 0, 0);

          for (unsigned int i = 0; i <  std::min(this->protein_predicted_cs.size(), this->chemshifts_from_file.size()); i++) {
               for (unsigned int j = 0; j < 6; j++) {

                    if ((std::fabs(this->protein_predicted_cs[i][j]) > 0.0001)
                         && (std::fabs(this->chemshifts_from_file[i][j]) > 0.0001)
                         && !(std::isnan(this->chemshifts_from_file[i][j]))
                         && !(std::isnan(this->protein_predicted_cs[i][j]))) {
                        
                         const double diff  = this->protein_predicted_cs[i][j] -  this->chemshifts_from_file[i][j];
                         chi_sq[j] += diff * diff;
                         n_bins[j] += 1;
                    }
               }
          }

          std::vector<double> rmsd = vector_utils::make_vector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

          for (unsigned int j = 0; j < 6; j++) {
               rmsd[j] = chi_sq[j] / n_bins[j];
          }

          this->n_print += 1;

          if (this->n_print > 500) {
              std::cout << "# CHEMSHIFT -- THREAD# " << this->thread_id
                        << "  -- WEIGHTS: " << this->weights
                        << "  -- RMSD: " << rmsd
                        << std::endl;

              n_print = 0;
         }
     }


     double log_gauss(const std::vector< std::vector<double > > &cs1_this, 
                      const std::vector< std::vector<double > > &cs2_this, 
                      const std::vector<double> &sigma_this) {

          std::vector<double> chi_sq = vector_utils::make_vector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
          std::vector<int> n_bins = vector_utils::make_vector(0, 0, 0, 0, 0, 0);

          for (unsigned int i = 0; i <  std::min(cs1_this.size(), cs2_this.size()); i++) {
               for (unsigned int j = 0; j < 6; j++) {

                   if ((std::fabs(cs1_this[i][j]) > 0.0001)
                    && (std::fabs(cs2_this[i][j]) > 0.0001)
                    && !(std::isnan(cs1_this[i][j]))
                    && !(std::isnan(cs2_this[i][j]))) {

                        const double diff = cs1_this[i][j] - cs2_this[i][j];
                        chi_sq[j] += diff * diff;
                        n_bins[j] += 1;

                    }
               }
          }

          double energy = 0.0;

          for (unsigned int i = 0; i < 6; i++) {
                energy += 0.5 * chi_sq[i] / (sigma_this[i] * sigma_this[i])
                        + (double)(n_bins[i] + 1) * std::log(Math<double>::sqrt_2_pi* sigma_this[i]);
          }

          return energy;

    }


     double log_cauchy(const std::vector< std::vector<double > > &cs1_this, 
                       const std::vector< std::vector<double > > &cs2_this, 
                       const std::vector<double> &gamma_this) {

          double energy = 0.0;

          for (unsigned int j = 0; j < 6; j++) {
               for (unsigned int i = 0; i <  std::min(cs1_this.size(), cs2_this.size()); i++) {

                   if ((std::fabs(cs1_this[i][j]) > 0.0001)
                    && (std::fabs(cs2_this[i][j]) > 0.0001)
                    && !(std::isnan(cs1_this[i][j]))
                    && !(std::isnan(cs2_this[i][j]))) {

                         const double diff = cs1_this[i][j] - cs2_this[i][j];

                         energy += std::log(gamma_this[j]) 
                                 + std::log(1.0 + (diff * diff) / (gamma_this[j] * gamma_this[j]));

                    }
               }

               // Extra contribution from Jeffreys prior
               energy += std::log(gamma_this[j]);
          }

          return energy;

     }


     //! Calculate the energy using the MD-Energy function from Robustelli et al., Structure 18, 923–933.
     //! The flat_bottom_potential is set via a commandline option.
     //! \param predicted_cs A list of predicted chemical shifts
     //! \param experimental_cs A list of calculated chemical shifts
     //! \return energy The associated energy based on the Robustelli et al. MD energy
     double get_energy_robustelli(const std::vector< std::vector<double > > &predicted_cs, 
                                  const std::vector< std::vector<double > > &experimental_cs) {

          const double flat_bottom_tolerance = 0.4;
          double cs_energy = 0.0;

          for (unsigned int i = 0; i < predicted_cs.size(); i++) {

               unsigned int residue_type_index = 0;

               if ((*this->chain)[i].residue_type == PRO) {
                    residue_type_index = 2;
               } else if ((*this->chain)[i].residue_type == GLY) {
                    residue_type_index = 1;
               }

               for (unsigned int j = 0; j < 6; j++) {
                    // Chemical shifts with no predicted or experimental data are either 0.0 or NaN, and we need to skip these!
                    if ((std::fabs(predicted_cs[i][j]) > 0.0001)
                        && (std::fabs(experimental_cs[i][j]) > 0.0001)
                        && !(std::isnan(experimental_cs[i][j]))
                        && !(std::isnan(predicted_cs[i][j]))) {

                         const double epsilon = camshift_constants::flat_bottom_limit[residue_type_index][j];
                         const double x_zero = camshift_constants::end_harmonic_potential[j];
                         const double beta = camshift_constants::scale_harmonic_potential[j];
                         const double delta_cs = std::fabs(predicted_cs[i][j] - experimental_cs[i][j]);

                         if (delta_cs < epsilon * flat_bottom_tolerance) {
                              // Add nothing to energy

                         } else if (delta_cs < x_zero) {
                              cs_energy += std::pow((delta_cs - flat_bottom_tolerance * epsilon)/beta, 2.0);

                         } else {

                              cs_energy += std::pow((x_zero - flat_bottom_tolerance * epsilon)/beta, 2.0) 
                                   + std::tanh(2.0 * (x_zero - flat_bottom_tolerance * epsilon) * (delta_cs - x_zero) / (20.0 * beta * beta));
                         }
                    } 
               }
          }

          return cs_energy;
     }



     //! Local settings class
     const class Settings: public EnergyTerm<ChainFB>::Settings {

     public:

          //! The NMR-STAR filename
          std::string star_filename;

          //! Type of energy function
          int energy_type;

          //! Whether to sample weights
          bool sample_weights;

          //! Constructor
          Settings(std::string star_filename="",
                   int energy_type=1,
                   bool sample_weights=false)
              : star_filename(star_filename),
                energy_type(energy_type),
                sample_weights(sample_weights) {}

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "star-filename:" << settings.star_filename << "\n";
               o << "energy-type:" << settings.energy_type << "\n";
               o << "sample-weights:" << settings.sample_weights << "\n";
               o << static_cast<const typename EnergyTerm<ChainFB>::Settings>(settings);
               return o;
          }                    
     } settings;


     //! Runs the CamShift predictor on a protein structure. This is a virtual
     //! function, and the implementation is found in the derived classes, TermCamshift and TermCamshiftCached.
     //! \param chain The chain object
     //! \return A Matrix containing six chemical shifts for each residue
     virtual std::vector< std::vector<double> > predict(phaistos::ChainFB& chain)=0;


     //! Constructor
     TermCamshiftBase(ChainFB *chain, std::string name, const Settings &settings=Settings(),
                      RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, name, settings, random_number_engine),
            CamshiftBackendBase(chain),
            random_number_engine(random_number_engine),
            settings(settings) {

          // Reset none-move check
          this->none_move = false;

          // Set some initial weights
          this->weights = vector_utils::make_vector(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
          this->weights_previous = this->weights;

          // Reset print-counter
          this->n_print = 0;


          // Check that filename exists
          if (!file_exists(settings.star_filename)) {
               std::cerr << "ERROR (camshift): File \"" << settings.star_filename << "\" not found. Exiting\n";
               assert(false);
          }

          // Parse chemical shifts file into class variable 
          this->chemshifts_from_file = get_chemshifts_from_file(settings.star_filename);

          // Make initial chemical shift prediction
          this->protein_predicted_cs = predict_base(*(this->chain));
          this->protein_predicted_cs_previous = this->protein_predicted_cs;

     }

     // Copy constructor
     TermCamshiftBase(const TermCamshiftBase &other, 
                      RandomNumberEngine *random_number_engine,
                      int thread_index, ChainFB *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain),
            CamshiftBackendBase(chain),
          random_number_engine(other.random_number_engine),
          none_move(other.none_move),
          weights(other.weights),
          weights_previous(other.weights_previous),
          n_print(other.n_print),
          chemshifts_from_file(other.chemshifts_from_file),
          protein_predicted_cs(other.protein_predicted_cs),
          protein_predicted_cs_previous(other.protein_predicted_cs_previous),
          settings(other.settings) {

          this->thread_id = thread_index;
     } 


     //! Calculate the energy and updates the chemical shift table and weight parameters.
     //! \param move_info A MoveInfo object
     //! \return The energy associated with the chemical shift prediction.
     double evaluate(MoveInfo *move_info=NULL) {

          // Set none_move to false as default
          this->none_move = false;

          // Check if current move was a none-move
          if ((move_info) && (settings.sample_weights)) {
               if (move_info->modified_angles.empty() == true) {
                    this->none_move = true;
               }
          }

          // If this move is a none-move, update weights
          if (this->none_move) {
              this->weights = update_weights(this->weights, this->random_number_engine);

          // Otherwise update chemical shifts
          } else {
               this->protein_predicted_cs = predict(*(this->chain));
    
          } 


         // Return the energy according to the energy model
         return log_gauss(this->protein_predicted_cs, 
                          this->chemshifts_from_file,
                          this->weights);


     }


     //! Reject move and backup new parameters
     //! IMPORTANT: The CamShift cached term does NOT inherit this function
     void accept() {

          // If the move was a none-move, backup weights
          if (this->none_move) {
               this->weights_previous = this->weights;

          // Else it was a physical move, and backup chemical shifts
          } else {
               this->protein_predicted_cs_previous = this->protein_predicted_cs;
          }

          // Print weights
          print_weights();
     }


     //! Reject move and roll back new parameters
     //! IMPORTANT: The CamShift cached term does NOT inherit this function
     void reject() {

          // If the move was a none-move, backup weights
          if (this->none_move) {
               this->weights = this->weights_previous;

          // Else it was a physical move, and roll back chemical shifts
          } else {
               this->protein_predicted_cs = this->protein_predicted_cs_previous;
          }

          // Print weights
          print_weights();
     }



}; // End of TermCamshiftBase



} // End namespace phaistos

#endif
