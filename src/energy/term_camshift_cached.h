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

#ifndef TERM_CAMSHIFT_CACHED
#define TERM_CAMSHIFT_CACHED

#include "term_camshift_base.h"

namespace phaistos {

//! Cached version of the Camshift energy term
class TermCamshiftCached: public TermCamshiftBase<TermCamshiftCached> {

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::TermCamshiftBase<TermCamshiftCached> TermCamshiftBase;

     //! Define NodeType locally for ease of reference
     typedef ChainFB::ChainTree::NodeType NodeType;


     //! The type of value registered to the cached iterator in 
     //! each iteration step
     class CamshiftContributionType {

          //! Internal data
          double data[6];

     public:

          //! Constructor
          CamshiftContributionType() {
               data[0] = 0;
               data[1] = 0;
               data[2] = 0;
               data[3] = 0;
               data[4] = 0;
               data[5] = 0;
          }
          
          //! Overload [] indexing operator (const)
          double operator[](const int index) const {
               return data[index];
          }

          //! Overload indexing operator (non-const)
          double& operator[](const int index) {
               return data[index];
          }

          //! Overload += operator
          const CamshiftContributionType& operator+=(const CamshiftContributionType& c) {
               data[0] += c.data[0];
               data[1] += c.data[1];
               data[2] += c.data[2];
               data[3] += c.data[3];
               data[4] += c.data[4];
               data[5] += c.data[5];
               return *this;
          }          

          //! Overload output operator
          friend std::ostream& operator<<(std::ostream& out, const CamshiftContributionType &c) {
               out << c.data[0] << " " 
                   << c.data[1] << " " 
                   << c.data[2] << " " 
                   << c.data[3] << " " 
                   << c.data[4] << " " 
                   << c.data[5];
               return out;
          }

     };

     
     //! Datastructure used as a RETURN_TYPE in the cached iterator
     template <typename CONTRIBUTION_TYPE>
     class CamshiftReturnType: public VectorReturnType<ChainFB, double> {
          
     public:
          
          //! Define how data vector should be resized given a node
          void initialize(NodeType *node) {
               this->data.resize(6, 0.0);
          }
          
          //! Specifies how a contribution is stored in this type - given
          //! the current node in the chaintree and the entity type
          //! potential_duplicate is true in the case where 
          //! entity_self=entity_other (identical nodes), which
          //! is registered twice - the second with potential_duplicate=true, so
          //! it can be filtered out
          void register_contribution(const NodeType *node, 
                                     const NodeType *entity_self, 
                                     const NodeType *entity_other, 
                                     const CONTRIBUTION_TYPE &contribution,
                                     bool potential_duplicate) {
               if (!potential_duplicate) {
                    this->data[0] += contribution[0];
                    this->data[1] += contribution[1];
                    this->data[2] += contribution[2];
                    this->data[3] += contribution[3];
                    this->data[4] += contribution[4];
                    this->data[5] += contribution[5];
               }
          }
     };

     //! Cached iterator
     CachedIterator<chaintree::PairIterator<ChainFB,NodeType,NodeType>,
                    std::pair<CamshiftContributionType, CamshiftContributionType >,
                    std::vector<CamshiftReturnType<CamshiftContributionType > > > cached_it;

     //! Cached iterator settings
     chaintree::PairIterator<ChainFB,
                             NodeType,
                             NodeType>::Settings it_settings;


     //! Local lookup table from atom_type to local index (HA,  CA,  H,   N,   C,   CB)
     std::vector<int> atom_type_to_local_index;

public:

     //! Local settings class
     const class Settings: public TermCamshiftBase::Settings {
     } settings;

     //! Constructor
     //!
     //! \param chain Molecule chain object
     //! \param settings Local settings object
     TermCamshiftCached(ChainFB *chain, const Settings &settings=Settings(),
                        RandomNumberEngine *random_number_engine = &random_global)
          : TermCamshiftBase(chain, "camshift-cached", settings, random_number_engine),
            cached_it(*chain),
            atom_type_to_local_index(definitions::ATOM_ENUM_SIZE, -1),
            settings(settings) {

          using namespace phaistos;
          using namespace definitions;

          const bool only_modified_pairs = true;
          const int minimum_residue_distance = 0;
          const double cut_off = 5.0;

          it_settings = chaintree::PairIterator<ChainFB, NodeType,
                                                         NodeType>::Settings(cut_off,
                                                                             only_modified_pairs,
                                                                             minimum_residue_distance);
          
          atom_type_to_local_index[HA]  = 0;
          atom_type_to_local_index[HA3] = 0; //GLY
          atom_type_to_local_index[CA]  = 1;
          atom_type_to_local_index[H]   = 2;
          atom_type_to_local_index[N]   = 3;
          atom_type_to_local_index[C]   = 4;
          atom_type_to_local_index[CB]  = 5;
     }

     //! Copy constructor
     //!
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermCamshiftCached(const TermCamshiftCached &other, 
                        RandomNumberEngine *random_number_engine,
                        int thread_index, ChainFB *chain)
          : TermCamshiftBase(other, random_number_engine, thread_index, chain),
            cached_it(*chain),
            it_settings(other.it_settings),
            atom_type_to_local_index(other.atom_type_to_local_index),
            settings(other.settings) {


     } 
 
     //! Runs the CamShift predictor on a protein structure.
     //! \param chain The chain object
     //! \return A Matrix containing six chemical shifts for each residue
     std::vector< std::vector<double> >predict(phaistos::ChainFB& chain) {

          using namespace phaistos;
          using namespace definitions;

          std::vector<ResidueHBondData> h_bond_network = get_hydrogen_bonding_network(*(this->chain));
          std::vector<AromaticRing> ring_current_info = get_ring_current_info(*(this->chain));
          std::vector< std::vector<double> > temp_protein_predicted_cs;


          for(ResidueIterator<ChainFB> res1(*(this->chain)); !(res1).end(); ++res1) {

               //Don't calculate chemical shifts of first and last residue.
               if ((res1->terminal_status == NTERM) || (res1->terminal_status == CTERM)){
                    temp_protein_predicted_cs.push_back(vector_utils::make_vector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
                    continue;
               }

               //Constant
               std::vector<double> random_coil_contribution        = random_coil_contribution_all_residues[res1->index];
               //Local
               std::vector<double> dihedral_angle_contribution     = get_dihedral_angle_contribution(*res1);
               //Local
               std::vector<double> backbone_distance_contribution  = get_backbone_distance_contribution(*res1);
               //Local
               std::vector<double> current_side_chain_contribution = get_current_side_chain_contribution(*res1);
               //Local, residue n-1, n and n+1
               std::vector<double> extra_distances_contribution    = get_extra_distances_contribution((*res1), (*(this->chain)));

               std::vector<double> residue_predicted_cs = vector_utils::make_vector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

               for (unsigned int i = 0; i < 6; i++) {
                    double predicted_cs = backbone_distance_contribution[i]  
                                        + random_coil_contribution[i]        
                                        + current_side_chain_contribution[i] 
                                        + dihedral_angle_contribution[i]     
                                        + extra_distances_contribution[i];    

                    residue_predicted_cs[i] += predicted_cs;

               }

               temp_protein_predicted_cs.push_back(residue_predicted_cs);
         }


         for(ResidueIterator<ChainFB> res1(*(this->chain)); !(res1).end(); ++res1) {
               //Don't calculate chemical shifts of first and last residue.
               if ((res1->terminal_status == NTERM) || (res1->terminal_status == CTERM)) continue;

               //Global, no cut-off
               std::vector<double> ring_current_contribution     = get_ring_current_contribution((*res1), ring_current_info);
               //Global, no cut-off
               std::vector<double> hydrogen_bonding_contribution = get_hydrogen_bonding_contribution((*res1), h_bond_network);
               //Global, no cut-off
               std::vector<double> ss_bond_contribution          = get_ss_bond_contribution((*res1), (*(this->chain)));

               for (unsigned int i = 0; i < 6; i++) {
                    temp_protein_predicted_cs[res1->index][i] += hydrogen_bonding_contribution[i]
                                                          + ring_current_contribution[i]    
                                                          + ss_bond_contribution[i];        
               }
          }

          // Cut-off distance is 5 angstrom -- the squared cut_off is 25.
          const double cut_off2 = 25.0;

          // Iterate over modified pairs
          for (cached_it(*(this->chain), it_settings); !cached_it.end(); ++cached_it) {
               // Pairs are accessible as cached_it->first and cached_it->second
     
               // Create contribution objects (wrapped double[6] objects)
               CamshiftContributionType contribution1;
               CamshiftContributionType contribution2;

               for (unsigned int i=0; i<cached_it->first->size(); ++i) {
                    Atom *atom1 = (*cached_it->first)[i];
                    int local_index1 = atom_type_to_local_index[atom1->atom_type];

                    for (unsigned int j=0; j<cached_it->second->size(); ++j) {

                         Atom *atom2 = (*cached_it->second)[j];
                         int local_index2 = atom_type_to_local_index[atom2->atom_type];

                         if ((abs(atom1->residue->index - atom2->residue->index) > 1)  // CHANGED FOR SPEED
                            && ((local_index1 > -1) || (local_index2 > -1))) {         // CHANGED FOR SPEED


                              double pair_distance2 = (atom2->position - atom1->position).norm_squared();

                              // No need to evaluate the norm, when we can use norm squared and skip a
                              // sqrt() call.
                              if (pair_distance2 < cut_off2) {

                                   int residue1_type = 0;
                                   int residue2_type = 0;

                                   double pair_distance = sqrt(pair_distance2); 
                                   double inv_pair_distance = 1.0/pair_distance;
                                   double inv3_pair_distance = inv_pair_distance*inv_pair_distance*inv_pair_distance;

                                   if (atom1->residue->residue_type == GLY) {
                                        residue1_type = 1;
                                   } else if (atom1->residue->residue_type == PRO) {
                                        residue1_type = 2;
                                   }
                                   if (atom2->residue->residue_type == GLY) {
                                        residue2_type = 1;
                                   } else if (atom2->residue->residue_type == PRO) {
                                        residue2_type = 2;
                                   }

                                   if (local_index1 != -1) {


                                       contribution1[local_index1]+= spn_hybridization_all_atoms_cached[atom2->residue->index][atom2->index][residue1_type][local_index1]*pair_distance
                                                                     + spn_hybridization_all_atoms_cached[atom2->residue->index][atom2->index][residue1_type+3][local_index1]*inv3_pair_distance;
                                   }
                                   if (local_index2 != -1) {
                                       contribution2[local_index2]+= spn_hybridization_all_atoms_cached[atom1->residue->index][atom1->index][residue2_type][local_index2]*pair_distance
                                                                     + spn_hybridization_all_atoms_cached[atom1->residue->index][atom1->index][residue2_type+3][local_index2]*inv3_pair_distance;

                                   }
                              }
                         }                         
                    }
               }

               // Register contribution for this node pair to the cache
               cached_it.register_contribution(std::make_pair(contribution1, contribution2));

          }


          // Compute total value from cache
          std::vector<CamshiftReturnType<CamshiftContributionType> > &chemical_shifts = cached_it.compute_total();

          std::vector<std::vector<double> > five_angs_contribution_cached(this->chain->size(), vector_utils::make_vector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0));


          // Calculate chemical shift energy
          ChainFB::ChainTree *chaintree = this->chain->get_chain_tree();
          for (unsigned int i=0; i<chemical_shifts.size(); ++i) {
               NodeType *node = chaintree->nodes[i];

               std::vector<AtomEnum> chemical_shift_atom_types = vector_utils::make_vector(HA,CA,H,N,C,CB);
               for (unsigned int j=0; j<chemical_shift_atom_types.size(); ++j) {
                    AtomEnum atom_type = chemical_shift_atom_types[j];
                    if (node->has_atom(atom_type)) {
                         Residue *res = (*node)[atom_type]->residue;
                          if (res->residue_type != PRO) {
                              five_angs_contribution_cached[res->index][j] += chemical_shifts[i][j];
                         } else {
                         //Skip H and N for PRO
                              if ((j != 2) && (j != 3)) {
                                   five_angs_contribution_cached[res->index][j] += chemical_shifts[i][j];
                              }
                         }
                    }
               }

               //Special case for GLY
               if (node->has_atom(HA3)) {
                    Residue *res = (*node)[HA3]->residue;
                    five_angs_contribution_cached[res->index][0] += chemical_shifts[i][0];
               }
          }

          //
          for(ResidueIterator<ChainFB> res1(*(this->chain)); !(res1).end(); ++res1) {
               //Don't calculate chemical shifts of first and last residue.
               if ((res1->terminal_status == NTERM) || (res1->terminal_status == CTERM)) continue;
               for (unsigned int i = 0; i < 6; i++) {
                    if ((i == 5) //Skip CB for GLY
                      && (res1->residue_type == GLY)) {
                         // Make sure chemical shift is zero, and do not add anything
                         temp_protein_predicted_cs[res1->index][i] = 0.0;
                    } else if ( ((i == 2) || (i==3)) //Skip H and N (indexes 2 and 3) for PRO
                                                     && (res1->residue_type == PRO)) {
                         // Make sure chemical shift is zero, and do not add anything
                         temp_protein_predicted_cs[res1->index][i] = 0.0;
                    } else {
                         temp_protein_predicted_cs[res1->index][i] += five_angs_contribution_cached[res1->index][i];
                   }
              }
         }
         return temp_protein_predicted_cs;
     }

     //! Function to accept the cache in the this->accept()
     void accept_cache() {
          cached_it.accept();
     }

     //! Function to reject the cache in the this->reject()
     void reject_cache() {
          cached_it.reject();
     }

};


} // End namespace phaistos
#endif
