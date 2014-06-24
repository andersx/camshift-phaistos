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

#ifndef CAMSHIFT_PREDICTOR
#define CAMSHIFT_PREDICTOR


#include <iostream>
#include <vector>
#include <math.h>

#include "camshift_data.h"



namespace phaistos {

class CamshiftBackendBase {

public:

    

     //! Define the cutoff defining an H-bond in angstroms (Must be 3 angstroms)
     double h_bond_cutoff; 

     //! A table containing the hybridiztion of all atoms in the protein
     std::vector< std::vector< std::vector<double> > > spn_hybridization_all_atoms;

     //! A table cotaining the hybridization of all atoms in the protein (cached version)
     std::vector< std::vector< std::vector< std::vector<double> > > > spn_hybridization_all_atoms_cached;

     //! A table containing the radom coil values for each atom
     std::vector< std::vector<double> > random_coil_contribution_all_residues;

     //! List that holds all hydrogen bond donors
     std::vector<phaistos::Atom*> donor_cache;

     //! List that holds all hydrogen bond acceptor oxygen atoms
     std::vector<phaistos::Atom*> acceptor_o_cache;

     //! List that holds all hydrogen bond acceptor nitrogen atoms
     std::vector<phaistos::Atom*> acceptor_n_cache;

     // //! The NMR STAR format file parser 
     // //! \param star_filename Name of NMR STAR datafile
     // //! \return value_matrix A list of chemical shifts ordered by residue number
     std::vector< std::vector<double> > get_chemshifts_from_file(const std::string star_filename) {
     
          std::vector<std::vector<double> > value_matrix;
          std::vector<std::string> lines = file_to_string_vector(star_filename);
          int residue_index_current = uninitialized<int>();
     
          for (unsigned int i=0; i<lines.size(); ++i) {
               std::string line = lines[i];
               // Remove white space at beginning and end
               boost::trim(line);
               if (line.size() == 0)
                    continue;
               std::vector<std::string> tokens;
               boost::split(tokens, line, boost::is_any_of(" \t"), boost::token_compress_on);
     
               // Detect whether this is a new residue
               int residue_index = boost::lexical_cast<int>(tokens[1]);
               if (residue_index != residue_index_current) {
                    int offset= (residue_index-residue_index_current);
                    if (!is_initialized(residue_index_current)) {
                         offset = 1;
                    }

                    value_matrix.resize(value_matrix.size()+offset, std::vector<double>(6, UNINITIALIZED));
                    residue_index_current = residue_index;
               }
               std::string atom_type = tokens[3];

               // Attempt to convert CS value to a double
               double value = uninitialized<double>();
               if (tokens.size() > 5) {
                    std::stringstream ss(tokens[5]);
                    if ((ss >> value).fail() || !(ss >> std::ws).eof()) {
                         value = uninitialized<double>();
                    }
               }

               if (atom_type == "CA") {
                     value_matrix.back()[1] = value;
               } else if (atom_type == "CB") {
                     value_matrix.back()[5] = value;
               } else if (atom_type == "C") {
                     value_matrix.back()[4] = value;
               } else if (atom_type == "N") {
                      value_matrix.back()[3] = value;
               } else if (atom_type == "HA") {
                      value_matrix.back()[0] = value;
               } else if (atom_type == "H") {
                      value_matrix.back()[2] = value;
               } else {
                     // std::cerr << "WARNING - NMR-STAR parser: Skipping unknown atom type " << atom_type << std::endl;
               }
          }
     
          return value_matrix;
     }

     //! Get the positions for the correct atoms in the given residue, the atoms types
     //! differ for GLY, PRO and other residues.
     //! \param residue The residue containing the atoms in question
     std::vector<Vector_3D> get_prediction_atom_positions(phaistos::ResidueFB& residue){

          using namespace definitions; 
          if (residue.residue_type == PRO) {
               return vector_utils::make_vector((residue)[HA]->position,
                                                (residue)[CA]->position,
                                                //We use CD instead of H for PRO
                                                (residue)[CD]->position,
                                                (residue)[N]->position,
                                                (residue)[C]->position,
                                                (residue)[CB]->position);
          } else if (residue.residue_type == GLY) {
                                                //We use HA3 instead of HA for GLY
               return vector_utils::make_vector((residue)[HA3]->position, //TODO: SWAP?
                                                (residue)[CA]->position,
                                                (residue)[H]->position,
                                                (residue)[N]->position,
                                                (residue)[C]->position,
                                                //We use HA2 instead of CB for GLY
                                                (residue)[HA2]->position); //TODO: SWAP?
          } else {
               return vector_utils::make_vector((residue)[HA]->position,
                                                (residue)[CA]->position,
                                                (residue)[H]->position,
                                                (residue)[N]->position,
                                                (residue)[C]->position,
                                                (residue)[CB]->position);
          }

     }



     //! A class holding the variable needed for a ring current contribution calculation
     //! in the point-dipole model
     class AromaticRing {
          public:
               //! \param Intensity Holds a vector of the ring-current intensities for different atom types
               std::vector< std::vector<double> > Intensity;
               //! \param ring_center The geometric center of the aromatic ring
               Vector_3D ring_center;
               //! \param ring_normal_vector A normal vector to the ring plane
               Vector_3D ring_normal_vector;
     };



     //! Get information about all aromatic rings in the protein, and store it in a vector.
     //! \param chain The current state of the chain
     //! \return all_rings A vector of AromaticRing objects, corresponding to each aromatic side chain in the chain
     std::vector<AromaticRing> get_ring_current_info(phaistos::ChainFB& chain) {

          using namespace phaistos;
          using namespace definitions;

           std::vector<AromaticRing> all_rings;

           for(ResidueIterator<ChainFB> res1(chain); !res1.end(); ++res1) {
                     if ((*res1).residue_type == PHE){
                             AromaticRing current_ring;
                             Vector_3D CE1sc, CE2sc, CGsc, CZsc;
                             CE1sc=(*res1)[CE1]->position;
                             CE2sc=(*res1)[CE2]->position;
                             CGsc=(*res1)[CG]->position;
                             CZsc=(*res1)[CZ]->position;
                             Vector_3D Phe_Centroid = (CGsc + CZsc)*0.5;
                             Vector_3D Phe_Normal = cross_product((CE1sc - CGsc), (CE2sc - CGsc)).normalize();
                             Phe_Normal = Phe_Normal + Phe_Centroid;
                             current_ring.ring_center = Phe_Centroid;
                             current_ring.ring_normal_vector = Phe_Normal;
                             current_ring.Intensity = vector_utils::make_vector(
                                                      vector_utils::make_vector(0.0197106186493400, 0.0107214108381650, 0.0248125664310836, -0.0121892669011571, 0.0379475059735857, 0.0228891866580828),
                                                      vector_utils::make_vector(0.0187096537702738, 0.0109848666800219, 0.0235627731481730, -0.0095191355123960, 0.0367113289223075, 0.0000000000000000),
                                                      vector_utils::make_vector(0.0186948051817643, 0.0109772155419481, 0.0000000000000000,  0.0000000000000000, 0.0369794702372783, 0.0228539667967397));
                             all_rings.push_back(current_ring);
                     } else if ((*res1).residue_type == TYR){
                             AromaticRing current_ring;
                             Vector_3D CE1sc, CE2sc, CGsc, CZsc;
                             CE1sc=(*res1)[CE1]->position;
                             CE2sc=(*res1)[CE2]->position;
                             CGsc=(*res1)[CG]->position;
                             CZsc=(*res1)[CZ]->position;
                             Vector_3D Tyr_Centroid = (CGsc + CZsc)*0.5;
                             Vector_3D Tyr_Normal = cross_product((CE1sc - CGsc), (CE2sc - CGsc)).normalize();
                             Tyr_Normal = Tyr_Normal + Tyr_Centroid;
                             current_ring.ring_center = Tyr_Centroid;
                             current_ring.ring_normal_vector = Tyr_Normal;
                             current_ring.Intensity = vector_utils::make_vector(
                                                      vector_utils::make_vector(0.0187376017237899, 0.0066462048321479, 0.0190884308526385, -0.0325441337476974, 0.0321684127282076, 0.0220080120663250),
                                                      vector_utils::make_vector(0.0190572333905550, 0.0068350431109663, 0.0152576959095551, -0.0345050214545893, 0.0291118543424899, 0.0000000000000000),
                                                      vector_utils::make_vector(0.0185362301819640, 0.0068703094892360, 0.0000000000000000,  0.0000000000000000, 0.0333030889700796, 0.0210886896115598));
                             all_rings.push_back(current_ring);
                     } else if ((*res1).residue_type == HIS){
                             AromaticRing current_ring;
                             Vector_3D CE1sc, CGsc, CD2sc;
                             CE1sc=(*res1)[CE1]->position;
                             CD2sc=(*res1)[CD2]->position;
                             CGsc=(*res1)[CG]->position;
                             Vector_3D His_Middle = (CGsc + CD2sc)*0.5;
                             Vector_3D His_Centroid = His_Middle*0.5527864045 + CE1sc*0.447213955;
                             Vector_3D His_Normal = cross_product((CGsc - CE1sc),(CD2sc - CE1sc));
                             His_Normal.normalize();
                             His_Normal = His_Normal + His_Centroid;
                             current_ring.ring_center = His_Centroid;
                             current_ring.ring_normal_vector = His_Normal;
                             current_ring.Intensity = vector_utils::make_vector(
                                                      vector_utils::make_vector(0.0119130052996101, 0.0036751726504010, 0.0130313394708430, -0.0110219144681428, 0.0176256168220394, 0.0131119278040671),
                                                      vector_utils::make_vector(0.0116936288982352, 0.0028029776470729, 0.0128327339461995, -0.0118152881953194, 0.0169791834536970, 0.0000000000000000),
                                                      vector_utils::make_vector(0.0121214853108319, 0.0038377495012746, 0.0000000000000000,  0.0000000000000000, 0.0193006129047789, 0.0144458048289313));
                             all_rings.push_back(current_ring);
                     } else if ((*res1).residue_type == TRP){
                             AromaticRing current_ring1;
                             Vector_3D NE1sc, CE2sc, CGsc, CD2sc;
                             NE1sc=(*res1)[NE1]->position;
                             CE2sc=(*res1)[CE2]->position;
                             CD2sc=(*res1)[CD2]->position;
                             CGsc=(*res1)[CG]->position;
                             Vector_3D Trp5_Middle = (CGsc + CD2sc)*0.5;
                             Vector_3D Trp5_Centroid = Trp5_Middle*0.5527864045 + NE1sc*0.447213955;
                             Vector_3D Trp5_Normal = cross_product((CE2sc - CGsc),(NE1sc - CGsc));
                             Trp5_Normal.normalize();
                             Trp5_Normal = Trp5_Normal + Trp5_Centroid;
                             current_ring1.ring_center = Trp5_Centroid;
                             current_ring1.ring_normal_vector = Trp5_Normal;
                             current_ring1.Intensity = vector_utils::make_vector(
                                                       vector_utils::make_vector(0.0202471655040728, 0.0060401060685761, 0.0255547010749335, -0.0134924111011466, 0.0344568285650518, 0.0192569022832064),
                                                       vector_utils::make_vector(0.0198814513646570, 0.0079561244360142, 0.0213769453566236, -0.0133617574701984, 0.0330345615296924, 0.0000000000000000),
                                                       vector_utils::make_vector(0.0194632244095829, 0.0120661699378942, 0.0000000000000000,  0.0000000000000000, 0.0365632161068835, 0.0191179003674346));
                             all_rings.push_back(current_ring1);
                             AromaticRing current_ring2;
                             Vector_3D CE3sc, CZ2sc, CH2sc;
                             CE3sc=(*res1)[CE3]->position;
                             CZ2sc=(*res1)[CZ2]->position;
                             CH2sc=(*res1)[CH2]->position;
                             Vector_3D Trp6_Centroid = (CE3sc + CZ2sc)*0.5;
                             Vector_3D Trp6_Normal = cross_product((CH2sc - CE3sc), (CE2sc - CE3sc)).normalize();
                             Trp6_Normal = Trp6_Normal + Trp6_Centroid;
                             current_ring2.ring_center = Trp6_Centroid;
                             current_ring2.ring_normal_vector = Trp6_Normal;
                             current_ring2.Intensity = vector_utils::make_vector(
                                                       vector_utils::make_vector(0.0112872405047719, 0.0250234697642502, 0.0199792013870296, -0.0196035137369542, 0.0024794007624462, 0.0273835953537524),
                                                       vector_utils::make_vector(0.0116558200308564, 0.0225182549511629, 0.0224576802064048, -0.0205959718464885, 0.0041195063438243, 0.0000000000000000),
                                                       vector_utils::make_vector(0.0119456590380501, 0.0208549364318374, 0.0000000000000000,  0.0000000000000000, 0.0008892397581068, 0.0260999762717652));
                             all_rings.push_back(current_ring2);
                     }
             }

          return all_rings;
     }



     //! Get the ring current contribution from all aromatic rings in the chain acting upon one residue
     //! \param res1 The residue in question
     //! \param all_rings A list of all AromaticRing objects in the chain
     //! \return ring_current_contribution a vector of contributions for the residue in question
     std::vector<double> get_ring_current_contribution(phaistos::ResidueFB& res1, std::vector<AromaticRing>& all_rings) {

                using namespace phaistos;
                using namespace definitions;

               std::vector<Vector_3D> backbone_atoms_list = get_prediction_atom_positions((res1));
               std::vector<double> ring_current_contribution = vector_utils::make_vector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

               unsigned int residue_type_index;
               if (res1.residue_type == GLY) residue_type_index = 1;
               if (res1.residue_type == PRO) residue_type_index = 2;
               else residue_type_index = 0;

               for (unsigned int i = 0; i < backbone_atoms_list.size(); i++) {
                if ((res1.residue_type == PRO) && (i == 2)) {
                     //Skip H for PRO
                     ring_current_contribution[i] = 0.0;
                     continue;
                } else if ((res1.residue_type == PRO) && (i == 3)) {
                     //Skip N for PRO
                     ring_current_contribution[i] = 0.0;
                     continue;
                } else if ((res1.residue_type == GLY) && (i == 5)) {
                     //Skip CB for GLY
                     ring_current_contribution[i] = 0.0;
                     continue;
                }
                    for (unsigned int j = 0; j < all_rings.size(); j++) {
                         Vector_3D r_vector;
                         r_vector = backbone_atoms_list[i] - all_rings[j].ring_center;
                         double theta = calc_angle(backbone_atoms_list[i], all_rings[j].ring_center, all_rings[j].ring_normal_vector);
                         // A factor of 1000 is added, since the parameters listed in the paper are 1000x less than they're supposed to
                         ring_current_contribution[i] += (1.0 - 3.0*std::pow(cos(theta), 2))/std::pow(r_vector.norm(), 3)*(all_rings[j].Intensity[residue_type_index][i]) * 1000.0;
                    }
               }

          return ring_current_contribution;
     }

     //! Returns the number of atoms in each atom type group for which experimental and predicted chemical shifts
     //! are available.
     std::vector<int> get_number_of_atoms_per_atom_type(std::vector< std::vector<double > > &experimental_cs, phaistos::ChainFB& chain) {

          using namespace phaistos;
          using namespace definitions;

          std::vector<int> bin_counts = vector_utils::make_vector(0, 0, 0, 0, 0, 0);
          for(ResidueIterator<ChainFB> res1(chain); !(res1).end(); ++res1) {
               int residue_index = res1->index;
               if ((res1->terminal_status == NTERM) || (res1->terminal_status == CTERM)) {
                    // Add nothing
                    continue;
               } else if (res1->residue_type == PRO) {
                    for (unsigned int i = 0; i < 6; i++) {
                         // Skip N and H for PRO
                         if ((i == 2) || (i == 3)) {
                              // Add nothing
                         } else {
                           if ((experimental_cs[residue_index][i] != 0.0)
                           && !(std::isnan(experimental_cs[residue_index][i])))
                                bin_counts[i] += 1;
                         }
                    }
               } else if (res1->residue_type == GLY) {
                    for (unsigned int i = 0; i < 6; i++) {
                         // Skip CB for GLY
                         if (i == 5) {
                              // Add nothing
                         } else {
                           if ((experimental_cs[residue_index][i] != 0.0)
                           && !(std::isnan(experimental_cs[residue_index][i])))
                                bin_counts[i] += 1;
                         }
                    }
               } else {
                    for (unsigned int i = 0; i < 6; i++) {
                           if ((experimental_cs[residue_index][i] != 0.0)
                           && !(std::isnan(experimental_cs[residue_index][i])))
                                bin_counts[i] += 1;
                    }
               }
          }
          return bin_counts;
     }


     //! Get the hybridization of an atom, by atom-typing via its mass and counting the number of neighboring atoms
     //! \param atom The atom for with a hybridization calulation is desired
     //! \return n where n is the hybridization in an sp(n) scheme, eg. sp2, sp3, etc.
     int get_spn_hybridization(phaistos::Atom& atom) {

          using namespace phaistos;
          using namespace definitions;


          if (atom.mass == atom_h_weight) {
               return 3;
          } else if (atom.mass == atom_c_weight) {
               // covalent_neighbours contains an extra atom, representing a
               // pseudo side chain for CA. This is an explicit work around,
               // which *should* always work, regardless of future changes.
               if ((atom.atom_type == CA)) {
                    return 3;
               } else if (((atom.covalent_neighbours).size() == 3)) {
                    return 2;
               } else {
                    return 3;
               }

          } else if (atom.mass == atom_n_weight) {
              if (atom.atom_type == N) { // Amide nitrogen is, in fact sp3 hybridized.
                   return 3;
              } else if ((atom.residue->residue_type == TRP)
                      || (atom.residue->residue_type == HIS)
                      || (atom.residue->residue_type == ARG)) {
                   return 2;
              } else {
                   return 3;
              }

          } else if (atom.mass == atom_o_weight) {
               if(((atom.covalent_neighbours).size() == 2)) {
                    return 3;
               } else {
                    return 2;
               }
          } else if (atom.mass == atom_s_weight) {
               return 3;
          } else {
               std::cerr << "Unknown atom hybridization: mass = " << atom.mass << "   n_bonds = " << (atom.covalent_neighbours).size() << "   " << atom.atom_type <<std::endl;
               exit(1);
          }
     }



     //! Create a list of the coefficients pertaining to the specific atom type (sp2, sp3, etc) for use in the five_ang function.
     //! \param chain The protein chain
     //! \return spn_hybridization_all_atoms A list containg the hybridization specific coefficients for all atoms in the chain for use in the five_angs function
     std::vector< std::vector< std::vector<double> > > get_spn_hybridization_all_atoms(phaistos::ChainFB& chain) {

          using namespace phaistos;
          using namespace definitions;

          std::vector< std::vector<double> > power_1_coefficients_list_STD;
          std::vector< std::vector<double> > power_1_coefficients_list_GLY;
          std::vector< std::vector<double> > power_1_coefficients_list_PRO;
          std::vector< std::vector<double> > power_minus_3_coefficients_list_STD;
          std::vector< std::vector<double> > power_minus_3_coefficients_list_GLY;
          std::vector< std::vector<double> > power_minus_3_coefficients_list_PRO;

          for (AtomIterator<ChainFB, definitions::ALL> atom(chain); !atom.end(); ++atom) {

              std::vector<double> power_1_coefficients_STD;
              std::vector<double> power_1_coefficients_GLY;
              std::vector<double> power_1_coefficients_PRO;
              std::vector<double> power_minus_3_coefficients_STD;
              std::vector<double> power_minus_3_coefficients_GLY;
              std::vector<double> power_minus_3_coefficients_PRO;

              int hybridization = get_spn_hybridization(*atom);
              if (hybridization == 2) {
                   unsigned int atom_type_index; //May give compiler warning
                   if ((*atom).mass == atom_c_weight) {
                        atom_type_index = 0;
                   } else if ((*atom).mass == atom_n_weight) {
                        atom_type_index = 1;
                   } else if ((*atom).mass == atom_o_weight) {
                        atom_type_index = 2;
                   } else {
                        //Throw an error, if the atom does not have any correct atom typing by Camshift
                        std::cerr << "Unknown sp2 hybridization: m = " << (*atom).mass << "   N = " << ((*atom).covalent_neighbours).size() << std::endl;
                        exit(1);
                   }
                   power_1_coefficients_STD = camshift_constants::power_1_sp2_coefficients_STD[atom_type_index];
                   power_minus_3_coefficients_STD = camshift_constants::power_minus_3_sp2_coefficients_STD[atom_type_index];

                   power_1_coefficients_GLY = camshift_constants::power_1_sp2_coefficients_GLY[atom_type_index];
                   power_minus_3_coefficients_GLY = camshift_constants::power_minus_3_sp2_coefficients_GLY[atom_type_index];

                   power_1_coefficients_PRO = camshift_constants::power_1_sp2_coefficients_PRO[atom_type_index];
                   power_minus_3_coefficients_PRO = camshift_constants::power_minus_3_sp2_coefficients_PRO[atom_type_index];

              } else {
                   unsigned int atom_type_index; //May give compiler warning
                   if ((*atom).mass == atom_c_weight) {
                        atom_type_index = 0;
                   } else if ((*atom).mass == atom_h_weight) {
                        atom_type_index = 1;
                   } else if ((*atom).mass == atom_n_weight) {
                        atom_type_index = 2;
                   } else if ((*atom).mass == atom_o_weight) {
                        atom_type_index = 3;
                   } else if ((*atom).mass == atom_s_weight) {
                        atom_type_index = 4;
                   } else {
                        //Throw an error, if the atom does not have any correct atom typing by Camshift
                        std::cerr << "Unknown sp3 hybridization: m = " << (*atom).mass << "   Number of neighbours: " << ((*atom).covalent_neighbours).size() << std::endl;
                        exit(1);
                   }
                   power_1_coefficients_STD = camshift_constants::power_1_sp3_coefficients_STD[atom_type_index];
                   power_minus_3_coefficients_STD = camshift_constants::power_minus_3_sp3_coefficients_STD[atom_type_index];

                   power_1_coefficients_GLY = camshift_constants::power_1_sp3_coefficients_GLY[atom_type_index];
                   power_minus_3_coefficients_GLY = camshift_constants::power_minus_3_sp3_coefficients_GLY[atom_type_index];

                   power_1_coefficients_PRO = camshift_constants::power_1_sp3_coefficients_PRO[atom_type_index];
                   power_minus_3_coefficients_PRO = camshift_constants::power_minus_3_sp3_coefficients_PRO[atom_type_index];

              }

              power_1_coefficients_list_STD.push_back(power_1_coefficients_STD);
              power_1_coefficients_list_GLY.push_back(power_1_coefficients_GLY);
              power_1_coefficients_list_PRO.push_back(power_1_coefficients_PRO);
              power_minus_3_coefficients_list_STD.push_back(power_minus_3_coefficients_STD);
              power_minus_3_coefficients_list_GLY.push_back(power_minus_3_coefficients_GLY);
              power_minus_3_coefficients_list_PRO.push_back(power_minus_3_coefficients_PRO);

          }

          return vector_utils::make_vector(power_1_coefficients_list_STD, power_minus_3_coefficients_list_STD,
                                           power_1_coefficients_list_GLY, power_minus_3_coefficients_list_GLY,
                                           power_1_coefficients_list_PRO, power_minus_3_coefficients_list_PRO);

     }



     //! Create a list of the coefficients pertaining to the specific atom type (sp2, sp3, etc) for use in the cached five_ang function.
     //! \param chain The protein chain
     //! \return spn_hybridization_all_atoms A list containg the hybridization specific coefficients for all atoms in the chain for use in the five_angs function
     std::vector< std::vector< std::vector< std::vector<double> > > > get_spn_hybridization_all_atoms_cached(phaistos::ChainFB& chain) {

          using namespace phaistos;
          using namespace definitions;

          std::vector< std::vector< std::vector< std::vector<double> > > > chain_spn_info;

          for(ResidueIterator<ChainFB> res(chain); !res.end(); ++res) {

               std::vector< std::vector< std::vector<double> > > residue_spn_info;

               for (AtomIterator<ChainFB, definitions::ALL> atom(*res); !atom.end(); ++atom) {

                    std::vector< std::vector<double> > atom_spn_info;

                   int hybridization = get_spn_hybridization(*atom);
                   if (hybridization == 2) {
                        unsigned int atom_type_index; //May give compiler warning
                        if ((*atom).mass == atom_c_weight) {
                             atom_type_index = 0;
                        } else if ((*atom).mass == atom_n_weight) {
                             atom_type_index = 1;
                        } else if ((*atom).mass == atom_o_weight) {
                             atom_type_index = 2;
                        } else {
                             std::cerr << "Unknown sp2 hybridization: m = " << (*atom).mass << "   N = " << ((*atom).covalent_neighbours).size() << std::endl;
                             exit(1);
                        }

                        atom_spn_info = vector_utils::make_vector(camshift_constants::power_1_sp2_coefficients_STD[atom_type_index],
                                                                  camshift_constants::power_1_sp2_coefficients_GLY[atom_type_index],
                                                                  camshift_constants::power_1_sp2_coefficients_PRO[atom_type_index],
                                                                  camshift_constants::power_minus_3_sp2_coefficients_STD[atom_type_index],
                                                                  camshift_constants::power_minus_3_sp2_coefficients_GLY[atom_type_index],
                                                                  camshift_constants::power_minus_3_sp2_coefficients_PRO[atom_type_index]);

                   } else {
                        unsigned int atom_type_index; //May give compiler warning
                        if ((*atom).mass == atom_c_weight) {
                             atom_type_index = 0;
                        } else if ((*atom).mass == atom_h_weight) {
                             atom_type_index = 1;
                        } else if ((*atom).mass == atom_n_weight) {
                             atom_type_index = 2;
                        } else if ((*atom).mass == atom_o_weight) {
                             atom_type_index = 3;
                        } else if ((*atom).mass == atom_s_weight) {
                             atom_type_index = 4;
                        } else {
                             std::cerr << "Unknown sp3 hybridization: m = " << (*atom).mass << "   N = " << ((*atom).covalent_neighbours).size() << std::endl;
                             exit(1);
                        }

                        atom_spn_info = vector_utils::make_vector(camshift_constants::power_1_sp3_coefficients_STD[atom_type_index],
                                                                  camshift_constants::power_1_sp3_coefficients_GLY[atom_type_index],
                                                                  camshift_constants::power_1_sp3_coefficients_PRO[atom_type_index],
                                                                  camshift_constants::power_minus_3_sp3_coefficients_STD[atom_type_index],
                                                                  camshift_constants::power_minus_3_sp3_coefficients_GLY[atom_type_index],
                                                                  camshift_constants::power_minus_3_sp3_coefficients_PRO[atom_type_index]);
                   }
                   residue_spn_info.push_back(atom_spn_info);
              }
              chain_spn_info.push_back(residue_spn_info);
         }

         return chain_spn_info;
     }



     //! Get the contribution from all extra atoms within a 5 angstrom sphere from each atom in a residue
     //! \param current_residue The current residue
     //! \param spn_hybridization_all_atoms A list of the hybridization of all atoms in the chain
     //! \param chain The current state of the protein chain
     //! \return five_angs_contribution The contribution from all extra atoms within a five angstrom sphere from each atom in the residue
     std::vector<double> get_five_angs_contribution(phaistos::ResidueFB& current_residue, std::vector< std::vector< std::vector<double> > >& spn_hybridization_all_atoms, phaistos::ChainFB& chain) {

          using namespace phaistos;
          using namespace definitions;

          int atom_type_index;
          if (current_residue.residue_type == GLY) {
               atom_type_index = 2;
          } else if (current_residue.residue_type == PRO) {
               atom_type_index = 4;
          } else {
               atom_type_index = 0;
          }
          std::vector< std::vector<double> > power_1_coefficients_list = spn_hybridization_all_atoms[atom_type_index];
          std::vector< std::vector<double> > power_minus_3_coefficients_list = spn_hybridization_all_atoms[atom_type_index+1];

          std::vector<Vector_3D> current_residue_query_atom_positions = get_prediction_atom_positions(current_residue);
          std::vector<double> five_angs_contribution = vector_utils::make_vector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

          unsigned int iterator_atom_index = 0;
          unsigned int current_residue_index = current_residue.index;

          //This double loop is about 1/6 faster when iterations are carried out in this order.
          for (AtomIterator<ChainFB, definitions::ALL> atom(chain); !atom.end(); ++atom) {
               Vector_3D atom_position = (*atom).position;
               for (unsigned int i = 0; i < 6; i++) {
                    if ((current_residue.residue_type == PRO) && (i == 2)) {        //Skip H for Proline
                    } else if ((current_residue.residue_type == PRO) && (i == 3)) { //Skip N for Proline
                    } else if ((current_residue.residue_type == GLY) && (i == 5)) { //Skip CB for Glycine
                    } else if ((abs(((*atom).residue)->index - current_residue_index)  < 2)) { //Skip atoms on same or neighbour residues
                    } else {
                         double pair_distance = (current_residue_query_atom_positions[i] - atom_position).norm();
                         if (pair_distance < 5.0) {

                             five_angs_contribution[i] += ( power_1_coefficients_list[iterator_atom_index][i]*pair_distance
                                                          + power_minus_3_coefficients_list[iterator_atom_index][i]*std::pow(pair_distance, -3.0) );
                         }
                    }
               }
               iterator_atom_index++;
          }
          return five_angs_contribution;
     }


     //! Get contributions from a set of nearby atoms, as defined by Camshift
     //! \param current_residue The current residue
     //! \param chain The current state of the protein chain
     //! \return extra_distances_contribution A 6-vector of chemical shift contributions
     std::vector<double> get_extra_distances_contribution(phaistos::ResidueFB& current_residue, phaistos::ChainFB& chain) {

          using namespace phaistos;
          using namespace definitions;

          std::vector<double> extra_distances_contribution = vector_utils::make_vector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
          std::vector< std::vector<double> > extra_distances_coefficients;

          if (current_residue.residue_type == PRO) {
               extra_distances_coefficients = camshift_constants::extra_distances_coefficients_PRO;
          } else if (current_residue.residue_type == GLY) {
               extra_distances_coefficients = camshift_constants::extra_distances_coefficients_GLY;
          } else {
               extra_distances_coefficients = camshift_constants::extra_distances_coefficients_STD;
          }

          int current_residue_index = current_residue.index;

          for (unsigned int i = 0; i < extra_distances_coefficients.size(); i++) {
               if ((((chain)[current_residue_index + camshift_constants::extra_distances_index_a[i]]).has_atom(camshift_constants::extra_distances_atoms_a[i]))
                && (((chain)[current_residue_index + camshift_constants::extra_distances_index_b[i]]).has_atom(camshift_constants::extra_distances_atoms_b[i]))) {

                    double extra_distance = (((chain)[current_residue_index + camshift_constants::extra_distances_index_a[i]])[camshift_constants::extra_distances_atoms_a[i]]->position
                                           - ((chain)[current_residue_index + camshift_constants::extra_distances_index_b[i]])[camshift_constants::extra_distances_atoms_b[i]]->position).norm();

                    for (unsigned int j = 0; j < 6; j++) {
                         extra_distances_contribution[j] += extra_distance * extra_distances_coefficients[i][j];
                    }
               //TODO: Clean this up ... ugly code to take care of some special cases.
               //If we're asking about the CG atom, and residue only has CG1
               } else if ( (camshift_constants::extra_distances_atoms_a[i] == CG)
                && (((chain)[current_residue_index + camshift_constants::extra_distances_index_a[i]]).has_atom(CG1))
                && (((chain)[current_residue_index + camshift_constants::extra_distances_index_b[i]]).has_atom(camshift_constants::extra_distances_atoms_b[i]))) {
                     double extra_distance = (((chain)[current_residue_index + camshift_constants::extra_distances_index_a[i]])[CG1]->position
                                           - ((chain)[current_residue_index + camshift_constants::extra_distances_index_b[i]])[camshift_constants::extra_distances_atoms_b[i]]->position).norm();

                    for (unsigned int j = 0; j < 6; j++) {
                         extra_distances_contribution[j] += extra_distance * extra_distances_coefficients[i][j];
                    }
               } else if ( (camshift_constants::extra_distances_atoms_a[i] == CG)
                && (((chain)[current_residue_index + camshift_constants::extra_distances_index_a[i]]).has_atom(OG))
                && (((chain)[current_residue_index + camshift_constants::extra_distances_index_b[i]]).has_atom(camshift_constants::extra_distances_atoms_b[i]))) {
                     double extra_distance = (((chain)[current_residue_index + camshift_constants::extra_distances_index_a[i]])[OG]->position
                                           - ((chain)[current_residue_index + camshift_constants::extra_distances_index_b[i]])[camshift_constants::extra_distances_atoms_b[i]]->position).norm();

                    for (unsigned int j = 0; j < 6; j++) {
                         extra_distances_contribution[j] += extra_distance * extra_distances_coefficients[i][j];
                    }
               } else if ( (camshift_constants::extra_distances_atoms_a[i] == CG)
                && (((chain)[current_residue_index + camshift_constants::extra_distances_index_a[i]]).has_atom(OG1))
                && (((chain)[current_residue_index + camshift_constants::extra_distances_index_b[i]]).has_atom(camshift_constants::extra_distances_atoms_b[i]))) {
                     double extra_distance = (((chain)[current_residue_index + camshift_constants::extra_distances_index_a[i]])[OG1]->position
                                           - ((chain)[current_residue_index + camshift_constants::extra_distances_index_b[i]])[camshift_constants::extra_distances_atoms_b[i]]->position).norm();

                    for (unsigned int j = 0; j < 6; j++) {
                         extra_distances_contribution[j] += extra_distance * extra_distances_coefficients[i][j];
                    }
               } else if ( (camshift_constants::extra_distances_atoms_a[i] == CG)
                && (((chain)[current_residue_index + camshift_constants::extra_distances_index_a[i]]).has_atom(SG))
                && (((chain)[current_residue_index + camshift_constants::extra_distances_index_b[i]]).has_atom(camshift_constants::extra_distances_atoms_b[i]))) {
                     double extra_distance = (((chain)[current_residue_index + camshift_constants::extra_distances_index_a[i]])[SG]->position
                                           - ((chain)[current_residue_index + camshift_constants::extra_distances_index_b[i]])[camshift_constants::extra_distances_atoms_b[i]]->position).norm();

                    for (unsigned int j = 0; j < 6; j++) {
                         extra_distances_contribution[j] += extra_distance * extra_distances_coefficients[i][j];
                    }
                } else if ( (camshift_constants::extra_distances_atoms_b[i] == HA)
                && (((chain)[current_residue_index + camshift_constants::extra_distances_index_a[i]]).has_atom(camshift_constants::extra_distances_atoms_a[i]))
                && (((chain)[current_residue_index + camshift_constants::extra_distances_index_b[i]]).has_atom(HA3))) {
                     double extra_distance = (((chain)[current_residue_index + camshift_constants::extra_distances_index_a[i]])[camshift_constants::extra_distances_atoms_a[i]]->position
                                            - ((chain)[current_residue_index + camshift_constants::extra_distances_index_b[i]])[HA3]->position).norm();

                    for (unsigned int j = 0; j < 6; j++) {
                         extra_distances_contribution[j] += extra_distance * extra_distances_coefficients[i][j];
                    }
              }
          }

          return extra_distances_contribution;
     }



     //! The formula to which contributions due to backbone dihedral angles are fitted
     //! \param angle A dihedral backbone angle
     //! \param parameters A set of vectors corresponding to the dihedral and residuetype
     //! \return value The value of the dihedral angle contribution
     double dihedral_angle_function(const double angle, const std::vector<double> parameters) {

          double p1 = parameters[0];
          double p2 = parameters[1];
          double p3 = parameters[2];
          double p4 = parameters[3];
          double p5 = parameters[4];
          //NOTE: The formula below is different in the paper, but this is an error in the paper. The
          //formula below is correct (Source: personal communication with Kai Kolhoff), and consistent
          //with the formula implemented in CamShift 1.35.
          double value = p1 * cos(3.0*angle + p4) + p2* cos(angle + p5) + p3;

          return value;
     }

     //! Get contributions from a set of nearby atoms, as defined by Camshift
     //! \param current_residue The current residue
     //! \return dihedral_angle_contribution A 6-vector of chemical shift contributions
     std::vector<double> get_dihedral_angle_contribution(phaistos::ResidueFB& current_residue){

          using namespace phaistos;
          using namespace definitions;

          std::vector< std::vector<double> > dihedral_angle_coefficients;
          std::vector<double> backbone_angles;

          if (current_residue.residue_type == PRO) {
               dihedral_angle_coefficients = camshift_constants::dihedral_angle_coefficients_PRO;
               backbone_angles = vector_utils::make_vector(current_residue.get_phi(),
                                                           current_residue.get_psi(),
                                                           current_residue.get_sidechain_dof_values(SIDECHAIN_ATOMS)[0]); //CHI1 angle
          } else if (current_residue.residue_type == GLY) {
               dihedral_angle_coefficients = camshift_constants::dihedral_angle_coefficients_GLY;
               backbone_angles = vector_utils::make_vector(current_residue.get_phi(),
                                                           current_residue.get_psi(),
                                                           0.0); //CHI1 angle -- Set to zero for Glycine.
          } else {
               dihedral_angle_coefficients = camshift_constants::dihedral_angle_coefficients_STD;
               backbone_angles = vector_utils::make_vector(current_residue.get_phi(),
                                                           current_residue.get_psi(),
                                                           current_residue.get_sidechain_dof_values(SIDECHAIN_ATOMS)[0]); //CHI1 angle
          }

          std::vector<double> dihedral_angle_contribution = vector_utils::make_vector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);


          for (unsigned int i = 0; i < camshift_constants::dihedral_angle_data.size(); i++) {
               for (unsigned int j = 0; j < backbone_angles.size(); j++) {
                    //Skip dihedral contribution for CHI1 angle (j==2), if residue type is GLY or ALA
                      if (!((current_residue.residue_type == GLY) && j == 2)
                       && !((current_residue.residue_type == ALA) && j == 2)) {
                         dihedral_angle_contribution[i] += dihedral_angle_function(backbone_angles[j], camshift_constants::dihedral_angle_data[i][j])*dihedral_angle_coefficients[j][i];
                    }
               }
          }
          return dihedral_angle_contribution;

     }

     //! A class that holds all information about hydrogen bonding at the amide proton and carbonyl oxygen for a residue
     struct ResidueHBondData {

          bool has_O_bond;
          bool has_H_bond;

          double g1_O_degrees;
          double g2_O_degrees;
          double dist_O_angs;

          double g1_H_degrees;
          double g2_H_degrees;
          double dist_H_angs;

     };

     //! Lennard-Jones like formula, to which chemical shift perturbation due to hydrogen bonding is fitted
     //! \param x The hydrogen bonding distance (Formula S1, Camshift supplementary material)
     //! \param parameters The set of parameters corresponding to the particular hydrogenbond
     double h_bond_function(const double x, const std::vector<double> parameters ) {

          double p1 = parameters[0];
          double p2 = parameters[1];
          double p3 = parameters[2];
          double p4 = parameters[3];
          double p5 = parameters[4];
          double r  = parameters[5];
          double s  = parameters[6];

	  double value = p1 * (std::pow(p2/x + p4, r) - std::pow(p3/x + p5, s));

          return value;
     }


     //! Identify all hydrogen bonds to backbone atoms in the protein and store them
     //! \param chain The current protein
     //! \return h_bond_network A vector of struct ResidueHBondData items, one for each residue in the protein.
     std::vector<ResidueHBondData> get_hydrogen_bonding_network(const phaistos::ChainFB& chain) {
        
          using namespace phaistos;
          using namespace definitions;

          std::vector<ResidueHBondData> h_bond_network;

          double radians_to_degrees = 180.0 / 3.141592654;

          for(ResidueIterator<phaistos::ChainFB> res1(chain); !(res1).end(); ++res1) {

               ResidueHBondData res1_h_bond_data;

               res1_h_bond_data.has_H_bond = false;
               res1_h_bond_data.has_O_bond = false;

               double h_bond_min_distance_H2 = h_bond_cutoff*h_bond_cutoff;
               double h_bond_min_distance_O2 = h_bond_cutoff*h_bond_cutoff;

               bool res_has_H_and_N = ((res1->has_atom(H)) && (res1->has_atom(N)));
               bool res_has_O_and_C = ((res1->has_atom(O)) && (res1->has_atom(C)));


               if (res_has_O_and_C) {
                    for (std::vector<Atom*>::iterator atom1 = donor_cache.begin(); atom1 < donor_cache.end(); ++atom1) {

                         if (res1->index == (*atom1)->residue->index) continue;

                         // Calculate the distance to the atom, which must be less than 5 angstrom
                         double o_h_distance2 = ((*res1)[O]->position - (*atom1)->position).norm_squared();

                         if (o_h_distance2 < h_bond_min_distance_O2) {    // Minumum bond distance

                              // Get the atom, to which the hydrogen bonding H-atom is covalently bonded.
                              // There may be an easier way, but it works. :-/
                              AtomEnum h_covalent_atom = (*atom1)->covalent_neighbours[0].first;
                              Residue *res2 = (*atom1)->residue;
                              Atom *atom2 = (*res2)[h_covalent_atom];

                              double angle_coh = calc_angle((*res1)[C]->position, (*res1)[O]->position, (*atom1)->position) * radians_to_degrees;
                              double angle_ohc = calc_angle((*res1)[O]->position, (*atom1)->position, atom2->position) * radians_to_degrees;

                              if ((angle_coh > 90.0) && (angle_ohc > 90.0)) {

                                   res1_h_bond_data.has_O_bond = true;

                                   h_bond_min_distance_O2 = o_h_distance2;
                                   res1_h_bond_data.g1_O_degrees = angle_coh;
                                   res1_h_bond_data.g2_O_degrees = angle_ohc;
                                   res1_h_bond_data.dist_O_angs = sqrt(o_h_distance2);

                              }
                         }
                    }
               }

              if (res_has_H_and_N) {

                    // Search for bonds to oxygen atoms
                    for (std::vector<Atom*>::iterator atom1 = acceptor_o_cache.begin(); atom1 < acceptor_o_cache.end(); ++atom1) {

                         if (res1->index == (*atom1)->residue->index) continue;

                         double h_o_distance2 = ((*res1)[H]->position - (*atom1)->position).norm_squared();

                         if (h_o_distance2 < h_bond_min_distance_H2) {  // Minimum bond distance

                              AtomEnum h_covalent_atom = (*atom1)->covalent_neighbours[0].first; // Luckily the C-atom to which the O-atom is bonded is ALWAYS first
                              Residue *res2 = (*atom1)->residue;                                 // and on same residue. See phaistos/src/atom.cpp if in doubt.
                              Atom *atom2 = (*res2)[h_covalent_atom];

                              double angle_nho = calc_angle((*res1)[N]->position, (*res1)[H]->position, (*atom1)->position) * radians_to_degrees;
                              double angle_hoc = calc_angle((*res1)[H]->position, (*atom1)->position, atom2->position) * radians_to_degrees;

                              if ((angle_nho > 90.0) && (angle_hoc > 90.0)) {


                                   h_bond_min_distance_H2 = h_o_distance2;
                                   res1_h_bond_data.has_H_bond = true;
                                   res1_h_bond_data.g1_H_degrees = angle_nho;
                                   res1_h_bond_data.g2_H_degrees = angle_hoc;
                                   res1_h_bond_data.dist_H_angs = sqrt(h_o_distance2);
                              }
                         }
                    }


                    // Search for bonds to nitrogen atoms
                    for (std::vector<Atom*>::iterator atom1 = acceptor_n_cache.begin(); atom1 < acceptor_n_cache.end(); ++atom1) {

                         if (res1->index == (*atom1)->residue->index) continue;

                         double h_n_distance2 = ((*res1)[H]->position - (*atom1)->position).norm_squared();

                         if (h_n_distance2 < h_bond_min_distance_H2) {  // Minimum bond distance

                              Atom *atom2;

                              if (((*atom1)->atom_type == N) && ((*atom1)->residue->terminal_status != NTERM)) {

                                   atom2 = chain[(*atom1)->residue->index-1][C]; // Always exists if (atom1->residue->terminal_status != NTERM)

                              } else {

                                   AtomEnum h_covalent_atom = (*atom1)->covalent_neighbours[0].first;
                                   Residue *res2 = (*atom1)->residue;
                                   atom2 = (*res2)[h_covalent_atom];

                              }

                              double angle_nho = calc_angle((*res1)[N]->position, (*res1)[H]->position, (*atom1)->position) * radians_to_degrees;
                              double angle_hoc = calc_angle((*res1)[H]->position, (*atom1)->position, atom2->position) * radians_to_degrees;

                              if ((angle_nho > 90.0) && (angle_hoc > 90.0)) {

                                   h_bond_min_distance_H2 = h_n_distance2;
                                   res1_h_bond_data.has_H_bond = true;
                                   res1_h_bond_data.g1_H_degrees = angle_nho;
                                   res1_h_bond_data.g2_H_degrees = angle_hoc;
                                   res1_h_bond_data.dist_H_angs = sqrt(h_n_distance2);
                              }
                         }
                    }
               }

               h_bond_network.push_back(res1_h_bond_data);

          }

          return h_bond_network;
     }



     //! Get the contribution due to hydrogenbonding for a residue
     //! \param current_residue The residue for which hydrogen bonding contributions are to be calculated
     //! \param h_bond_network The vector holding all information about all hydrogen bonds in the protein
     //! \return hydrogen_bonding_contribution A 6-vector of chemical shift contributions
     std::vector<double> get_hydrogen_bonding_contribution(const phaistos::ResidueFB& current_residue, std::vector<ResidueHBondData>& h_bond_network) {

          using namespace phaistos;
          using namespace definitions;

          std::vector<double> hydrogen_bonding_contribution = vector_utils::make_vector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

          int current_residue_index = current_residue.index;
          std::vector< std::vector< std::vector<double> > > h_bond_coefficients;

          if (current_residue.residue_type == PRO) {
               h_bond_coefficients = camshift_constants::h_bond_coefficients_PRO;
          } else if (current_residue.residue_type == GLY) {
               h_bond_coefficients = camshift_constants::h_bond_coefficients_GLY;
          } else {
               h_bond_coefficients = camshift_constants::h_bond_coefficients_STD;
          }
          //contrib from n-1 residue (O)
          if (h_bond_network[current_residue_index-1].has_O_bond) {
               double d_contribution  = h_bond_function(h_bond_network[current_residue_index-1].dist_O_angs,  camshift_constants::h_bond_parameters[0]);
               double g1_contribution = h_bond_function(h_bond_network[current_residue_index-1].g2_O_degrees, camshift_constants::h_bond_parameters[1]); //note, that g2 and g1 are swapped
               double g2_contribution = h_bond_function(h_bond_network[current_residue_index-1].g1_O_degrees, camshift_constants::h_bond_parameters[2]);

               for (unsigned int i = 0; i< 6; i++) {
                    hydrogen_bonding_contribution[i] += d_contribution*h_bond_coefficients[0][0][i]
                                                     + g1_contribution*h_bond_coefficients[0][2][i]//Indexes are swapped in the CamShift code. This is also done here
                                                     + g2_contribution*h_bond_coefficients[0][1][i];
               }
          }

          //contrib fron n residue (H)
          if (h_bond_network[current_residue_index].has_H_bond) {
               double d_contribution  = h_bond_function(h_bond_network[current_residue_index].dist_H_angs,  camshift_constants::h_bond_parameters[0]);
               double g1_contribution = h_bond_function(h_bond_network[current_residue_index].g1_H_degrees, camshift_constants::h_bond_parameters[1]);
               double g2_contribution = h_bond_function(h_bond_network[current_residue_index].g2_H_degrees, camshift_constants::h_bond_parameters[2]);

               for (unsigned int i = 0; i< 6; i++) {
                    hydrogen_bonding_contribution[i] += d_contribution*h_bond_coefficients[1][0][i]
                                                     + g1_contribution*h_bond_coefficients[1][2][i]//Indexes are swapped in the CamShift code. This is also done here
                                                     + g2_contribution*h_bond_coefficients[1][1][i];
               }
          }

          //contrib from n residue (O)
          if (h_bond_network[current_residue_index].has_O_bond) {
               double d_contribution  = h_bond_function(h_bond_network[current_residue_index].dist_O_angs,  camshift_constants::h_bond_parameters[0]);
               double g1_contribution = h_bond_function(h_bond_network[current_residue_index].g2_O_degrees, camshift_constants::h_bond_parameters[1]); //note, that g2 and g1 are swapped
               double g2_contribution = h_bond_function(h_bond_network[current_residue_index].g1_O_degrees, camshift_constants::h_bond_parameters[2]);

               for (unsigned int i = 0; i< 6; i++) {
                    hydrogen_bonding_contribution[i] += d_contribution*h_bond_coefficients[2][0][i]
                                                     + g1_contribution*h_bond_coefficients[2][2][i]//Indexes are swapped in the CamShift code. This is also done here
                                                     + g2_contribution*h_bond_coefficients[2][1][i];
               }
          }

          //contrib from n+1 residue (H)
          if (h_bond_network[current_residue_index+1].has_H_bond) {
               double d_contribution  = h_bond_function(h_bond_network[current_residue_index+1].dist_H_angs,  camshift_constants::h_bond_parameters[0]);
               double g1_contribution = h_bond_function(h_bond_network[current_residue_index+1].g1_H_degrees, camshift_constants::h_bond_parameters[1]);
               double g2_contribution = h_bond_function(h_bond_network[current_residue_index+1].g2_H_degrees, camshift_constants::h_bond_parameters[2]);

               for (unsigned int i = 0; i< 6; i++) {
                    hydrogen_bonding_contribution[i] += d_contribution*h_bond_coefficients[3][0][i]
                                                     + g1_contribution*h_bond_coefficients[3][2][i]//Indexes are swapped in the CamShift code. This is also done here
                                                     + g2_contribution*h_bond_coefficients[3][1][i];
               }
          }

          return hydrogen_bonding_contribution;
     }


     //! Get certain side chain atoms fom a residue
     //! \param current_residue The residue in which the side chain is found
     //! \return side_chain_atoms The list of side chain atoms
     std::vector<phaistos::definitions::AtomEnum> get_side_chain_atoms(const phaistos::ResidueFB& current_residue){

          using namespace phaistos;
          using namespace definitions;

          std::vector<AtomEnum> side_chain_atoms;

          if      ((current_residue).residue_type == ALA) side_chain_atoms = camshift_constants::side_chain_ALA_atoms;
          else if ((current_residue).residue_type == ARG) side_chain_atoms = camshift_constants::side_chain_ARG_atoms;
          else if ((current_residue).residue_type == ASN) side_chain_atoms = camshift_constants::side_chain_ASN_atoms;
          else if ((current_residue).residue_type == ASP) side_chain_atoms = camshift_constants::side_chain_ASP_atoms;
          else if ((current_residue).residue_type == CYS) side_chain_atoms = camshift_constants::side_chain_CYS_atoms;
          else if ((current_residue).residue_type == GLN) side_chain_atoms = camshift_constants::side_chain_GLN_atoms;
          else if ((current_residue).residue_type == GLU) side_chain_atoms = camshift_constants::side_chain_GLU_atoms;
          else if ((current_residue).residue_type == GLY) side_chain_atoms = camshift_constants::side_chain_GLY_atoms;
          else if ((current_residue).residue_type == HIS) side_chain_atoms = camshift_constants::side_chain_HIS_atoms;
          else if ((current_residue).residue_type == ILE) side_chain_atoms = camshift_constants::side_chain_ILE_atoms;
          else if ((current_residue).residue_type == LEU) side_chain_atoms = camshift_constants::side_chain_LEU_atoms;
          else if ((current_residue).residue_type == LYS) side_chain_atoms = camshift_constants::side_chain_LYS_atoms;
          else if ((current_residue).residue_type == MET) side_chain_atoms = camshift_constants::side_chain_MET_atoms;
          else if ((current_residue).residue_type == PHE) side_chain_atoms = camshift_constants::side_chain_PHE_atoms;
          else if ((current_residue).residue_type == PRO) side_chain_atoms = camshift_constants::side_chain_PRO_atoms;
          else if ((current_residue).residue_type == SER) side_chain_atoms = camshift_constants::side_chain_SER_atoms;
          else if ((current_residue).residue_type == THR) side_chain_atoms = camshift_constants::side_chain_THR_atoms;
          else if ((current_residue).residue_type == TRP) side_chain_atoms = camshift_constants::side_chain_TRP_atoms;
          else if ((current_residue).residue_type == TYR) side_chain_atoms = camshift_constants::side_chain_TYR_atoms;
          else if ((current_residue).residue_type == VAL) side_chain_atoms = camshift_constants::side_chain_VAL_atoms;
          else {
               std::cerr << "Unknown residue type: " << current_residue.residue_type << " ... CRITICAL FAIL IN get_side_chain_atoms" << std::endl;
               exit(1);
          }

          return side_chain_atoms;
     }



     //! Get coefficients sets for side chain atoms in the side chain contribution term
     //! \param current_residue The residue for which the side chain atom coordinates are desired
     //! \return side_chain_data a vector containing the coefficients for the important atoms in the side chain
     std::vector< std::vector<double> > get_side_chain_data(const phaistos::ResidueFB& current_residue){

          using namespace phaistos;
          using namespace definitions;

          std::vector< std::vector<double> > side_chain_data;

          if      (current_residue.residue_type == ALA) side_chain_data = camshift_constants::side_chain_ALA_data;
          else if (current_residue.residue_type == ARG) side_chain_data = camshift_constants::side_chain_ARG_data;
          else if (current_residue.residue_type == ASN) side_chain_data = camshift_constants::side_chain_ASN_data;
          else if (current_residue.residue_type == ASP) side_chain_data = camshift_constants::side_chain_ASP_data;
          else if (current_residue.residue_type == CYS) side_chain_data = camshift_constants::side_chain_CYS_data;
          else if (current_residue.residue_type == GLN) side_chain_data = camshift_constants::side_chain_GLN_data;
          else if (current_residue.residue_type == GLU) side_chain_data = camshift_constants::side_chain_GLU_data;
          else if (current_residue.residue_type == GLY) side_chain_data = camshift_constants::side_chain_GLY_data;
          else if (current_residue.residue_type == HIS) side_chain_data = camshift_constants::side_chain_HIS_data;
          else if (current_residue.residue_type == ILE) side_chain_data = camshift_constants::side_chain_ILE_data;
          else if (current_residue.residue_type == LEU) side_chain_data = camshift_constants::side_chain_LEU_data;
          else if (current_residue.residue_type == LYS) side_chain_data = camshift_constants::side_chain_LYS_data;
          else if (current_residue.residue_type == MET) side_chain_data = camshift_constants::side_chain_MET_data;
          else if (current_residue.residue_type == PHE) side_chain_data = camshift_constants::side_chain_PHE_data;
          else if (current_residue.residue_type == PRO) side_chain_data = camshift_constants::side_chain_PRO_data;
          else if (current_residue.residue_type == SER) side_chain_data = camshift_constants::side_chain_SER_data;
          else if (current_residue.residue_type == THR) side_chain_data = camshift_constants::side_chain_THR_data;
          else if (current_residue.residue_type == TRP) side_chain_data = camshift_constants::side_chain_TRP_data;
          else if (current_residue.residue_type == TYR) side_chain_data = camshift_constants::side_chain_TYR_data;
          else if (current_residue.residue_type == VAL) side_chain_data = camshift_constants::side_chain_VAL_data;
          else {
               std::cerr << "Unknown residue type: " << current_residue.residue_type << " ... CRITICAL FAIL IN get_current_side_chain_contribution" << std::endl;
               exit(1);
          }

          return side_chain_data;
     }


     //! Get the chemical shift perturbation due to the atoms in the side chain of the current residue
     //! \param The residue for which the side chain chemical shift perturbation is to be calculated
     //! \return A 6-vector of chemical shift contributions
     std::vector<double> get_current_side_chain_contribution(phaistos::ResidueFB& current_residue){

          using namespace phaistos;
          using namespace definitions;

          std::vector<Vector_3D> current_residue_query_atom_positions = get_prediction_atom_positions(current_residue);
          std::vector<double> current_side_chain_contribution = vector_utils::make_vector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
          std::vector<AtomEnum> side_chain_atoms = get_side_chain_atoms(current_residue);
          std::vector< std::vector<double> > side_chain_data = get_side_chain_data(current_residue);

          for (unsigned int j = 0; j < 6; j++) {
                if ((current_residue.residue_type == PRO) && (j == 2)) {
                     //Skip H for PRO
                     current_side_chain_contribution[j] = 0.0;
                     continue;
                } else if ((current_residue.residue_type == PRO) && (j == 3)) {
                     //Skip N for PRO
                     current_side_chain_contribution[j] = 0.0;
                     continue;
                } else if ((current_residue.residue_type == GLY) && (j == 5)) {
                     //Skip CB for GLY
                     current_side_chain_contribution[j] = 0.0;
                     continue;
                }
                for (unsigned int i = 0; i < side_chain_atoms.size(); i++) {
                     if (current_residue.has_atom(side_chain_atoms[i])) {
                          current_side_chain_contribution[j] += (((current_residue)[side_chain_atoms[i]])->position - current_residue_query_atom_positions[j]).norm()*side_chain_data[i][j];
                     //TODO: The check below is redundant and should to be removed in the production
                     //      version, but it is here for now, so I can test that I don't ask for wrong
                     //      atom types.
                     } else {
                         std::cerr << "Unable to find atom:  " << side_chain_atoms[i] << "  in residue:" << current_residue.residue_type << " " << current_residue.index << std::endl;
                         exit(1);
                    }
               }
          }

          return current_side_chain_contribution;
     }



     //! Get positions of the backbone atoms of the n-1, n, n+1 residues, for the chemical shift perturbations
     //! due to inter-atomic distances to those distances. This function takes into account different atom names
     //! in GLY and PRO. The order of returned atoms corresponds to a coefficient vector in camshift_data.
     //! \param current_residue The n'th residue
     //! \return ordered_backbone_positions_list A vector of Vector_3D's
     std::vector<Vector_3D> get_ordered_backbone_positions_list(phaistos::ResidueFB& current_residue) {

          using namespace phaistos;
          using namespace definitions;

          Residue *previous_residue = current_residue.get_neighbour(-1);
          Residue *next_residue     = current_residue.get_neighbour(+1);
          //The order of the atoms in this list is important. It must follow the category 2 table
          //in the supplementary material.
          std::vector<Vector_3D> ordered_backbone_positions_list;

          //N, CA, HA, C, O from last residue (HA -> HA3 for GLY)
          ordered_backbone_positions_list.push_back(((*previous_residue)[N])->position);
          ordered_backbone_positions_list.push_back(((*previous_residue)[CA])->position);
          if (previous_residue->residue_type == GLY) {
               ordered_backbone_positions_list.push_back(((*previous_residue)[HA3])->position);
          } else {
               ordered_backbone_positions_list.push_back(((*previous_residue)[HA])->position);
          }
          ordered_backbone_positions_list.push_back(((*previous_residue)[C])->position);
          ordered_backbone_positions_list.push_back(((*previous_residue)[O])->position);

          //N, H, CA, HA, C, O from current residue (HA -> HA3 for GLY, H -> CD for PRO)
          ordered_backbone_positions_list.push_back(((current_residue)[N])->position);
          if (current_residue.residue_type == PRO) {
               ordered_backbone_positions_list.push_back(((current_residue)[CD])->position);
          } else {
               ordered_backbone_positions_list.push_back(((current_residue)[H])->position);
          }
          ordered_backbone_positions_list.push_back(((current_residue)[CA])->position);
          if (current_residue.residue_type == GLY) {
               ordered_backbone_positions_list.push_back(((current_residue)[HA3])->position);
          } else {
               ordered_backbone_positions_list.push_back(((current_residue)[HA])->position);
          }
          ordered_backbone_positions_list.push_back(((current_residue)[C])->position);
          ordered_backbone_positions_list.push_back(((current_residue)[O])->position);

          //N, H, CA, HA, C from next residue (HA -> HA3 for GLY, H -> CD for PRO)
          ordered_backbone_positions_list.push_back(((*next_residue)[N])->position);
          if (next_residue->residue_type == PRO) {
               ordered_backbone_positions_list.push_back(((*next_residue)[CD])->position);
          } else {
               ordered_backbone_positions_list.push_back(((*next_residue)[H])->position);
          }
          ordered_backbone_positions_list.push_back(((*next_residue)[CA])->position);
          if (next_residue->residue_type == GLY) {
               ordered_backbone_positions_list.push_back(((*next_residue)[HA3])->position);
          } else {
               ordered_backbone_positions_list.push_back(((*next_residue)[HA])->position);
          }
          ordered_backbone_positions_list.push_back(((*next_residue)[C])->position);

          return ordered_backbone_positions_list;
     }


     //! Get the chemical shift contributions form distances to n-1, n, n+1 backbone atoms
     //! \param
     //! \return
     std::vector<double> get_backbone_distance_contribution(phaistos::ResidueFB& current_residue){

          using namespace phaistos;
          using namespace definitions;

          std::vector<double> backbone_distance_contributions = vector_utils::make_vector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

          std::vector< std::vector<double> > backbone_distances_data;

          if (current_residue.residue_type == PRO) {
               backbone_distances_data = camshift_constants::backbone_distances_data_PRO;
          } else if (current_residue.residue_type == GLY) {
               backbone_distances_data = camshift_constants::backbone_distances_data_GLY;
          } else {
               backbone_distances_data = camshift_constants::backbone_distances_data_STD;
          }
          std::vector<Vector_3D> current_residue_query_atom_positions = get_prediction_atom_positions(current_residue);

          std::vector<Vector_3D> ordered_backbone_positions_list = get_ordered_backbone_positions_list((current_residue));

          for (unsigned int i = 0; i < 6; i++) {
                if ((current_residue.residue_type == PRO) && (i == 2)) { //Skip H for PRO
                } else if ((current_residue.residue_type == PRO) && (i == 3)) { //Skip N for PRO
                } else if ((current_residue.residue_type == GLY) && (i == 5)) { //Skip CB for GLY
                } else {
                    for (unsigned int j = 0; j < ordered_backbone_positions_list.size(); j++) {
                         backbone_distance_contributions[i] += (ordered_backbone_positions_list[j] - current_residue_query_atom_positions[i]).norm()*backbone_distances_data[i][j];
                    }
                }
          }
          return backbone_distance_contributions;
     }



     //! Get random coil chemical shifts for a given residue type from the Camshift dataset
     //! \param query_residue The residue for which random coil values are to be looked up
     //! \return random_coil_residue_data A 6-vector of random coil values
     std::vector<double> get_random_coil_data_vector(const phaistos::Residue *query_residue) {

          using namespace phaistos;
          using namespace definitions;

          std::vector<double> random_coil_residue_data;

          if      (query_residue->residue_type == ALA) random_coil_residue_data = camshift_constants::random_coil_data[0];
          else if (query_residue->residue_type == ARG) random_coil_residue_data = camshift_constants::random_coil_data[1];
          else if (query_residue->residue_type == ASN) random_coil_residue_data = camshift_constants::random_coil_data[2];
          else if (query_residue->residue_type == ASP) random_coil_residue_data = camshift_constants::random_coil_data[3];
          else if (query_residue->residue_type == CYS) random_coil_residue_data = camshift_constants::random_coil_data[4];
          else if (query_residue->residue_type == GLN) random_coil_residue_data = camshift_constants::random_coil_data[5];
          else if (query_residue->residue_type == GLU) random_coil_residue_data = camshift_constants::random_coil_data[6];
          else if (query_residue->residue_type == GLY) random_coil_residue_data = camshift_constants::random_coil_data[7];
          else if (query_residue->residue_type == HIS) random_coil_residue_data = camshift_constants::random_coil_data[8];
          else if (query_residue->residue_type == ILE) random_coil_residue_data = camshift_constants::random_coil_data[9];
          else if (query_residue->residue_type == LEU) random_coil_residue_data = camshift_constants::random_coil_data[10];
          else if (query_residue->residue_type == LYS) random_coil_residue_data = camshift_constants::random_coil_data[11];
          else if (query_residue->residue_type == MET) random_coil_residue_data = camshift_constants::random_coil_data[12];
          else if (query_residue->residue_type == PHE) random_coil_residue_data = camshift_constants::random_coil_data[13];
          else if (query_residue->residue_type == PRO) random_coil_residue_data = camshift_constants::random_coil_data[14];
          else if (query_residue->residue_type == SER) random_coil_residue_data = camshift_constants::random_coil_data[15];
          else if (query_residue->residue_type == THR) random_coil_residue_data = camshift_constants::random_coil_data[16];
          else if (query_residue->residue_type == TRP) random_coil_residue_data = camshift_constants::random_coil_data[17];
          else if (query_residue->residue_type == TYR) random_coil_residue_data = camshift_constants::random_coil_data[18];
          else if (query_residue->residue_type == VAL) random_coil_residue_data = camshift_constants::random_coil_data[19];
          else {
               std::cerr << "Unknown residue type: " << query_residue->residue_type << " ... CRITICAL FAIL IN get_random_coil_data_vector" << std::endl;
               exit(1);
          }

         return random_coil_residue_data;

     }



     //! Get the random_coil value for a give nucleus/atom
     //! \param query_atom An atom/nuclus
     //! \return random_coil_values The random coil value for a given atom
     double get_random_coil_value(phaistos::Atom *query_atom) {

          using namespace phaistos;
          using namespace definitions;

          double random_coil_value = 0.0;
          std::vector<double> random_coil_residue_data;
          random_coil_residue_data = get_random_coil_data_vector(query_atom->residue);

          if      (((*query_atom).atom_type) == HA) random_coil_value = random_coil_residue_data[0];
          else if (((*query_atom).atom_type) == CA) random_coil_value = random_coil_residue_data[1];
          else if (((*query_atom).atom_type) == H ) random_coil_value = random_coil_residue_data[2];
          //N has a correction, based on the previous residue_type.
          else if (((*query_atom).atom_type) == N ) random_coil_value = random_coil_residue_data[3]
                                                                      + get_random_coil_data_vector(query_atom->residue->get_neighbour(-1))[6];
          else if (((*query_atom).atom_type) == C ) random_coil_value = random_coil_residue_data[4];
          else if (((*query_atom).atom_type) == CB) random_coil_value = random_coil_residue_data[5];
          //Special option for GLY
          else if (((*query_atom).atom_type) == HA2) random_coil_value = random_coil_residue_data[0];
          else {
               std::cerr << "Unknown atom type: " << (*query_atom) << " ... CRITICAL FAIL IN get_random_coil_value" << std::endl;
               exit(1);
          }

          return random_coil_value;

     }



     //! Get random coil values for all six atoms in the residue
     //! \param res1 The residue for which random_coil values are to be looked up
     std::vector<double> get_random_coil_data_residue(phaistos::ResidueFB& res1){

          using namespace phaistos;
          using namespace definitions;

          if (res1.residue_type == PRO) {
               return vector_utils::make_vector(get_random_coil_value((res1)[HA]),
                                                get_random_coil_value((res1)[CA]),
                                                0.0,
                                                0.0,
                                                get_random_coil_value((res1)[C]),
                                                get_random_coil_value((res1)[CB]));
          } else if (res1.residue_type == GLY) {
               return vector_utils::make_vector(get_random_coil_value((res1)[HA2]),
                                                get_random_coil_value((res1)[CA]),
                                                get_random_coil_value((res1)[H]),
                                                get_random_coil_value((res1)[N]),
                                                get_random_coil_value((res1)[C]),
                                                0.0);
          } else {
               return vector_utils::make_vector(get_random_coil_value((res1)[HA]),
                                                get_random_coil_value((res1)[CA]),
                                                get_random_coil_value((res1)[H]),
                                                get_random_coil_value((res1)[N]),
                                                get_random_coil_value((res1)[C]),
                                                get_random_coil_value((res1)[CB]));
          }

     }


     //! Make the initial table of random coil values for all atoms in the prediction model
     //! \param chain The protein chain
     //! \return random_coil_data_all_residues A 6xN vector of random coil values
     std::vector< std::vector<double> > get_random_coil_contribution_all_residues(const phaistos::ChainFB& chain) {

          using namespace phaistos;
          using namespace definitions;

          std::vector< std::vector<double> > random_coil_data_all_residues;

          for(ResidueIterator<phaistos::ChainFB> res1(chain); !(res1).end(); ++res1) {

               if (res1->terminal_status == NTERM) {
                    random_coil_data_all_residues.push_back(vector_utils::make_vector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
               } else if (res1->terminal_status == CTERM) {
                    random_coil_data_all_residues.push_back(vector_utils::make_vector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
               } else if (res1->residue_type == PRO){
                    random_coil_data_all_residues.push_back(vector_utils::make_vector(get_random_coil_value((*res1)[HA]),
                                           get_random_coil_value((*res1)[CA]),
                                           0.0,
                                           0.0,
                                           get_random_coil_value((*res1)[C]),
                                           get_random_coil_value((*res1)[CB])));

               } else if (res1->residue_type == GLY)  {
                    random_coil_data_all_residues.push_back(vector_utils::make_vector(get_random_coil_value((*res1)[HA2]),
                                           get_random_coil_value((*res1)[CA]),
                                           get_random_coil_value((*res1)[H]),
                                           get_random_coil_value((*res1)[N]),
                                           get_random_coil_value((*res1)[C]),
                                           0.0));

               } else  {
                    random_coil_data_all_residues.push_back(vector_utils::make_vector(get_random_coil_value((*res1)[HA]),
                                           get_random_coil_value((*res1)[CA]),
                                           get_random_coil_value((*res1)[H]),
                                           get_random_coil_value((*res1)[N]),
                                           get_random_coil_value((*res1)[C]),
                                           get_random_coil_value((*res1)[CB])));

               }
          }
          return random_coil_data_all_residues;
     }



     //! Get the chemical shift contribution from CYS-bridges
     //! \param res1 the residue for which the presence of a CYS bridge is to be checked.
     //! \param chain The protein chain
     //! \return ss_bond_contribution A 6-vector of chemical shift contributions
     std::vector<double> get_ss_bond_contribution(phaistos::ResidueFB& res1, phaistos::ChainFB& chain) {

          using namespace phaistos;
          using namespace definitions;

          std::vector<double> ss_bond_contribution = vector_utils::make_vector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

          //Bail out early when query residue is of wrong type or does not have the appropriate SG atom.
          if ((res1.residue_type != CYS) || !(res1.has_atom(SG)))
               return vector_utils::make_vector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

          for(ResidueIterator<phaistos::ChainFB> res2(chain); !(res2).end(); ++res2) {
               if (( res2->residue_type == CYS)
                && (res2->has_atom(SG))) {
                    double ss_distance = (((*res2)[SG])->position - (res1[SG])->position).norm();
                    //Allow the SS_bond to deviate from the ideal bond length (2.05 angstrom) by
                    //a generous 0.3 angstrom. This is our "arbitrary" S-S bond criterium.
                    if ((std::fabs(ss_distance - 2.05) < 0.3) && (res1.index != res2->index)){
                         for (unsigned int i = 0; i < 6; i++) {
                             //A factor of 0.5 is missing in the published formula or parameters. We apply it here,
                             //and obtain results identical to those of CamShift 1.35.
                             ss_bond_contribution[i] += ss_distance * camshift_constants::ss_bond_coefficients[i] * 0.5;
                         }
                         return ss_bond_contribution;
                    }
               }
          }
          return ss_bond_contribution;
     }

     CamshiftBackendBase () {
     }

     //! Constructor for CamshiftBackendBase
     //! Camshift predictors inhereit from this class
       CamshiftBackendBase(phaistos::ChainFB *chain) {

           using namespace phaistos;
          using namespace definitions;

          for (int i = 0; i < chain->size(); i++) {
               for (int j = 0; j < (*chain)[i].size(); j++) {

                   Atom* atom1 = (*chain)[i][j];

                    if (atom1->mass == atom_h_weight) {
                         donor_cache.push_back(atom1);
                    } else if (atom1->mass == atom_o_weight) {
                         acceptor_o_cache.push_back(atom1);
                    } else if (atom1->mass == atom_n_weight) {
                         acceptor_n_cache.push_back(atom1);
                    }
               }
          }


          h_bond_cutoff = 3.0;

          spn_hybridization_all_atoms = get_spn_hybridization_all_atoms(*chain);
          spn_hybridization_all_atoms_cached = get_spn_hybridization_all_atoms_cached(*chain);
          random_coil_contribution_all_residues = get_random_coil_contribution_all_residues(*chain);
     }




     //! Runs the CamShift predictor on a protein structure.
     //! \param chain The chain object
     //! \return A Matrix containing six chemical shifts for each residue
     std::vector< std::vector<double> >predict_base(phaistos::ChainFB& chain) {

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


}; //End class CamshiftBackendBase

} // End namespace phaistos

#endif
