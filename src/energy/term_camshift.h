// Copyright (C) 2011 by Anders Steen Christensen, Wouter Boomsma
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

#ifndef TERM_CAMSHIFT
#define TERM_CAMSHIFT

#include "term_camshift_base.h"

namespace phaistos {

class TermCamshift: public TermCamshiftBase<TermCamshift> {

public:

     // For convenience, define local EnergyTermCommon
     typedef phaistos::TermCamshiftBase<TermCamshift> TermCamshiftBase;

     // Use settings from base class
     typedef TermCamshiftBase::Settings Settings;

     //! Runs the CamShift predictor on a protein structure (uncached version)
     //! \param chain The chain object
     //! \return A Matrix containing six chemical shifts for each residue
     std::vector< std::vector<double> > predict(phaistos::ChainFB& chain) {

          using namespace phaistos;
          using namespace definitions;

          std::vector<ResidueHBondData> h_bond_network = get_hydrogen_bonding_network(chain);
          std::vector<AromaticRing> ring_current_info = get_ring_current_info(chain);
          std::vector< std::vector<double> > protein_predicted_cs;

          for (ResidueIterator<ChainFB> res1(chain); !(res1).end(); ++res1) {

              //Don't calculate chemical shifts of first and last residue.
              if ((res1->terminal_status == NTERM) || (res1->terminal_status == CTERM)) {
                        protein_predicted_cs.push_back(vector_utils::make_vector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
                        continue;
              }

              std::vector<double> random_coil_contribution        = random_coil_contribution_all_residues[res1->index];
              std::vector<double> dihedral_angle_contribution     = get_dihedral_angle_contribution(*res1);
              std::vector<double> backbone_distance_contribution  = get_backbone_distance_contribution(*res1);
              std::vector<double> current_side_chain_contribution = get_current_side_chain_contribution(*res1);
              std::vector<double> extra_distances_contribution    = get_extra_distances_contribution((*res1), chain);
              std::vector<double> five_angs_contribution          = get_five_angs_contribution((*res1), spn_hybridization_all_atoms, chain);
              std::vector<double> ring_current_contribution       = get_ring_current_contribution((*res1), ring_current_info);
              std::vector<double> hydrogen_bonding_contribution   = get_hydrogen_bonding_contribution((*res1), h_bond_network);
              // std::vector<double> ss_bond_contribution            = get_ss_bond_contribution((*res1), chain);
              std::vector<double> residue_predicted_cs            = vector_utils::make_vector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

              for (unsigned int i = 0; i < 6; i++) {
                   residue_predicted_cs[i] += backbone_distance_contribution[i]
                                       + random_coil_contribution[i]
                                       + current_side_chain_contribution[i]
                                       + dihedral_angle_contribution[i]
                                       + extra_distances_contribution[i]
                                       + five_angs_contribution[i]
                                       + hydrogen_bonding_contribution[i]
                                       + ring_current_contribution[i];
                                       // + ss_bond_contribution[i];


               }

               protein_predicted_cs.push_back(residue_predicted_cs);
          }

          return protein_predicted_cs;
     }

     //Constructor
     TermCamshift(ChainFB *chain, 
                  const Settings &settings=Settings(),
                  RandomNumberEngine *random_number_engine = &random_global)
          : TermCamshiftBase(chain, "camshift", settings, random_number_engine) {

     }

     //Copy constructor
     TermCamshift(const TermCamshift &other, 
                  RandomNumberEngine *random_number_engine,
                  int thread_index, ChainFB *chain)
          : TermCamshiftBase(other, random_number_engine, thread_index, chain) {
     } 

     // Dummy function to avoid duplicate accept() code
     void accept_cache() {}

     // Dummy function to avoid duplicate reject() code
     void reject_cache() {}

};

//! Observable specialization for TermCamshift
template <>
class Observable<TermCamshift>: public TermCamshift, public ObservableBase {

public:

     //! Local settings class.
     const class Settings: public TermCamshift::Settings, public ObservableBase::Settings {
     public:

          //! Constructor. Defines default values for settings object.
          Settings(){}

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << static_cast<const typename TermCamshift::Settings>(settings);
               o << static_cast<const ObservableBase::Settings>(settings);
               return o;
          }
     } settings; //!< Local settings object·

     //! Constructor.
     //! \param energy_term TermCamshift energy term object
     //! \param settings Local Settings object
     //! \param reference_energy_function All observables have a pointer to a reference energy function which they can refer to.
     Observable(const TermCamshift &energy_term,
                const ObservableBase::Settings &settings=ObservableBase::Settings(),
                Energy<ChainFB> *reference_energy_function=NULL)
          : TermCamshift(energy_term),
            settings(dynamic_cast<const Settings&>(settings)) {
     }

     //! Copy Constructor.
     //! \param other Source object from which copy is made
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     Observable(const Observable &other, int thread_index, ChainFB *chain)
          : TermCamshift(other, random_number_engine, thread_index, chain),
            settings(other.settings) {
     }


     //! Clone: Corresponds to a virtual copy constructor
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermCamshift *clone(int thread_index=0, typename TermCamshift::ChainType *chain=NULL) {
          return new Observable<TermCamshift>(*this, thread_index, chain);
     }


     //! Chemical shift RMSD between two sets of chemical shifts
     //! \param cs1_this A matrix containing chemical shifts
     //! \param cs2_this A matrix containing chemical shifts
     //! \return A vector containing RMSDs
     std::vector<double> calc_rmsds(const std::vector< std::vector<double > > &cs1_this, 
                                    const std::vector< std::vector<double > > &cs2_this) {

          std::vector<double> chi_sq = vector_utils::make_vector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
          std::vector<int> n_bins = vector_utils::make_vector(0, 0, 0, 0, 0, 0);

          for (unsigned int i = 0; i <  std::min(cs1_this.size(), cs2_this.size()); i++) {
               for (unsigned int j = 0; j < 6; j++) {

                   if ((std::fabs(cs1_this[i][j]) > 0.0001)
                    && (std::fabs(cs2_this[i][j]) > 0.0001)
                    && !(std::isnan(cs1_this[i][j]))
                    && !(std::isnan(cs2_this[i][j]))) {
                        
                         const double diff  = cs1_this[i][j] - cs2_this[i][j];
                         chi_sq[j] += diff * diff;
                         n_bins[j] += 1;
                    }
               }
          }

          std::vector<double> rmsds = vector_utils::make_vector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

          for (unsigned int j = 0; j < 6; j++) {
               rmsds[j] = std::sqrt(chi_sq[j] / n_bins[j]);
          }

          return rmsds;

     }




     //! Make observation.
     virtual std::string observe(MoveInfo *move_info=NULL,
                                 PHAISTOS_LONG_LONG current_iteration=0,
                                 bool register_only=false) {

          // Energy to be returned
          double energy = this->evaluate();

          // Calculate new chemical shifts
          this->protein_predicted_cs = predict_base(*(this->chain));

          // Calculate RMSDs
          std::vector<double> rmsds = calc_rmsds(this->protein_predicted_cs,
                                                 this->chemshifts_from_file);

          // Output stream
          std::stringstream s;
          s << std::fixed << std::setprecision(5) << energy << "\t" << rmsds;

          return s.str();

     }

}; // End class Observable




} // End namespace phaistos
#endif
