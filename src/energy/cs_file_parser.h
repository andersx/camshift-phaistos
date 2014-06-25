#include <iostream>
#include <vector>

namespace camshift_parser {

     //! The NMR STAR format file parser 
     //! \param star_filename Name of NMR STAR datafile
     //! \return value_matrix A list of chemical shifts ordered by residue number
     std::vector< std::vector<double> > value_matrix_from_starfile(const std::string& star_filename) {

          using namespace phaistos;
          using namespace definitions;

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


     //! Detect if there are mismatches between input chemical shifts data and the actual sequence
     bool do_data_and_chain_match(const std::string& star_filename, phaistos::ChainFB& chain) {

          using namespace phaistos;
          using namespace definitions;

          std::vector<std::string> lines = file_to_string_vector(star_filename);
          bool contains_bugs = false;

          for (unsigned int i=0; i<lines.size(); ++i) {
               std::string line = lines[i];
               // Remove white space at beginning and end
               boost::trim(line);
               if (line.size() == 0)
                    continue;
               std::vector<std::string> tokens;
               boost::split(tokens, line, boost::is_any_of(" \t"), boost::token_compress_on);

               int residue_index = boost::lexical_cast<int>(tokens[1]);
               std::string res_type = tokens[2];

               for (ResidueIterator<phaistos::ChainFB> res1(chain); !(res1).end(); ++res1) {

                    if (residue_index == (*res1).index + 1) {
                         // std::cout << "Loading residue:" << (*res1).residue_type << (*res1).index << std::endl;
                         if (((*res1).residue_type == definitions::ALA) && (res_type != "ALA")) contains_bugs = true;
                         if (((*res1).residue_type == definitions::CYS) && (res_type != "CYS")) contains_bugs = true;
                         if (((*res1).residue_type == definitions::ASP) && (res_type != "ASP")) contains_bugs = true;
                         if (((*res1).residue_type == definitions::GLU) && (res_type != "GLU")) contains_bugs = true;
                         if (((*res1).residue_type == definitions::PHE) && (res_type != "PHE")) contains_bugs = true;
                         if (((*res1).residue_type == definitions::GLY) && (res_type != "GLY")) contains_bugs = true;
                         if (((*res1).residue_type == definitions::HIS) && (res_type != "HIS")) contains_bugs = true;
                         if (((*res1).residue_type == definitions::ILE) && (res_type != "ILE")) contains_bugs = true;
                         if (((*res1).residue_type == definitions::LYS) && (res_type != "LYS")) contains_bugs = true;
                         if (((*res1).residue_type == definitions::LEU) && (res_type != "LEU")) contains_bugs = true;
                         if (((*res1).residue_type == definitions::MET) && (res_type != "MET")) contains_bugs = true;
                         if (((*res1).residue_type == definitions::ASN) && (res_type != "ASN")) contains_bugs = true;
                         if (((*res1).residue_type == definitions::PRO) && (res_type != "PRO")) contains_bugs = true;
                         if (((*res1).residue_type == definitions::GLN) && (res_type != "GLN")) contains_bugs = true;
                         if (((*res1).residue_type == definitions::ARG) && (res_type != "ARG")) contains_bugs = true;
                         if (((*res1).residue_type == definitions::SER) && (res_type != "SER")) contains_bugs = true;
                         if (((*res1).residue_type == definitions::THR) && (res_type != "THR")) contains_bugs = true;
                         if (((*res1).residue_type == definitions::VAL) && (res_type != "VAL")) contains_bugs = true;
                         if (((*res1).residue_type == definitions::TRP) && (res_type != "TRP")) contains_bugs = true;
                         if (((*res1).residue_type == definitions::TYR) && (res_type != "TYR")) contains_bugs = true;
                    }
                    if (contains_bugs) {
                          std::cerr << "ERROR (chemshift): CS-data mismatch in " << (*res1).residue_type << (*res1).index + 1 << std::endl;
                          return false;
                    }
               }
          }

          
          return true;
     }



} //End namespace camshift_parser
