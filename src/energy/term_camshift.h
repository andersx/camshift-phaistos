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

} // End namespace phaistos
#endif
