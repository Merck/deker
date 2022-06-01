/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//MIT License
//
//Copyright (c) 2021 Merck Sharp & Dohme Corp. a subsidiary of Merck & Co., Inc., Kenilworth, NJ, USA.
//Originally authored by Sean M.S. Hayes (MSD)
//
//Permission is hereby granted, free of charge, to any person obtaining a copy
//of this software and associated documentation files (the "Software"), to deal
//in the Software without restriction, including without limitation the rights
//to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//copies of the Software, and to permit persons to whom the Software is
//furnished to do so, subject to the following conditions:
//  
//The above copyright notice and this permission notice shall be included in all
//copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//SOFTWARE.
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef DEKER_IO_OUTPUT
#define DEKER_IO_OUTPUT
///////////////////////////////////////////
#include <deker/fit/output.h>
#include <deker/io/misc.h>
#include <iostream>
#include <fstream>
///////////////////////////////////////////
namespace deker{
  namespace fit{
    struct Output_deker_IO:Output_deker{
      //included from Output_deker:
      // double lambda_regularization;
      // double sigma_kernel_width;
      // double BIC;
      // int return_status;
      // Eigen::VectorXd v;
      // std::vector<unsigned> feature_index;
      // std::vector<unsigned> key_to_original_data;
      double null_BIC;
      std::vector<std::string> feature_names;
      unsigned which_response;
      std::string response_name;
      std::string data_filename;
      //
      Eigen::VectorXd lgbm_imp;
      double lgbm_l2_fit;
      /////////////////////////////////////
      //constructor using a base Output_deker
      Output_deker_IO():Output_deker(){}
      Output_deker_IO(const Output_deker& out_deker):Output_deker(out_deker){}
      /////////////////////////////////////
      Output_deker_IO(const Output_deker& out_deker,
                      const std::vector<std::string>& full_feature_names,
                      const std::string& data_filename,
                      const unsigned& which_response,
                      const double& null_BIC,
                      const Eigen::Ref<const Eigen::VectorXd> lgbm_imp_in,
                      const double& lgbm_l2_fit):
        Output_deker(out_deker),
        data_filename(data_filename),
        which_response(which_response),
        response_name(full_feature_names.at(which_response)),
        null_BIC(null_BIC),
        lgbm_l2_fit(lgbm_l2_fit)
        {
          lgbm_imp.resize(feature_index.size());
          for(unsigned i = 0; i<feature_index.size(); i++){
            feature_names.push_back(full_feature_names.at(key_to_original_data.at(feature_index.at(i))));
            lgbm_imp(i) = lgbm_imp_in(feature_index.at(i));
          }
        }
      ////////////////////////////////////////////////////
      void serialize(std::ifstream& infile)
      {
        Output_deker::serialize(infile);
        feature_names.clear();
        infile.read((char*) (&null_BIC),sizeof(null_BIC));
        size_t n_features;
        infile.read((char*) (&n_features),sizeof(n_features));
        unsigned temp_index;
        size_t temp_name_size;
        std::string temp_name;
        for(unsigned i = 0; i<n_features;i++){
          infile.read((char*) (&temp_name_size), sizeof(temp_name_size));
          temp_name.resize(temp_name_size);
          infile.read(&temp_name[0], temp_name_size);
          feature_names.push_back(temp_name);
        }
        infile.read((char*) (&temp_name_size), sizeof(temp_name_size));
        data_filename.resize(temp_name_size);
        infile.read(&data_filename[0],temp_name_size);
        infile.read((char*) (&which_response), sizeof(which_response));
        infile.read((char*) (&temp_name_size), sizeof(temp_name_size));
        response_name.resize(temp_name_size);
        infile.read(&response_name[0],temp_name_size);
        //
        infile.read((char*) (&lgbm_l2_fit), sizeof(lgbm_l2_fit));
        Eigen::MatrixXd::Index vec_size;
        infile.read((char*) (&vec_size),sizeof(Eigen::MatrixXd::Index));
        lgbm_imp.resize(vec_size);
        infile.read((char*) lgbm_imp.data(),vec_size*sizeof(Eigen::MatrixXd::Scalar));
      }
      ////////////////////////////////////////////////////
      void serialize(std::ofstream& outfile) const
      {
        Output_deker::serialize(outfile);
        outfile.write((char*) (&null_BIC), sizeof(null_BIC));
        size_t n_features = feature_names.size();
        outfile.write((char*) (&n_features), sizeof(n_features));
        for(unsigned i = 0;i<n_features;i++){
          size_t size=feature_names.at(i).size();
          outfile.write((char*) (&size),sizeof(size));
          outfile.write(feature_names.at(i).c_str(), size);
        }
        size_t data_filename_size = data_filename.size();
        outfile.write((char*) (&data_filename_size),sizeof(data_filename_size));
        outfile.write(data_filename.c_str(),data_filename_size);
        outfile.write((char*) (&which_response),sizeof(which_response));
        size_t response_name_size = response_name.size();
        outfile.write((char*) (&response_name_size),sizeof(response_name_size));
        outfile.write(response_name.c_str(),response_name_size);
        //
        outfile.write((char*) (&lgbm_l2_fit), sizeof(lgbm_l2_fit));
        Eigen::MatrixXd::Index vec_size=lgbm_imp.size();
        outfile.write((char*) (&vec_size), sizeof(Eigen::MatrixXd::Index));
        outfile.write((char*) lgbm_imp.data(), vec_size*sizeof(Eigen::MatrixXd::Scalar));
      }
      /////////////////////////////////////////////////////
      void write_to_csv_header(std::ofstream& outfile,
                               const std::vector<std::string>& extra_identifiers)
      {
        for(unsigned i = 0; i<extra_identifiers.size();i++)
          outfile<<extra_identifiers.at(i)<<",";
        outfile<<"data_filename,"<<
          "response,"<<
            "lambda,"<<
              "sigma,"<<
                "df,"<<
                  "logli,"<<
                    "BIC,"<<
                      "null_BIC,"<<
                        "return_status,"<<
                          "lgbm_l2_fit,"<<
                            "feature,"<<
                              "weight,"<<
                                "lgbm_imp\n";
      }
      /////////////////////////////////////////////////////
      void write_to_csv(std::ofstream& outfile,
                        const std::vector<std::string>& extra_identifiers)
      {
        for(unsigned i = 0; i<feature_names.size();i++){
          for(unsigned j = 0; j<extra_identifiers.size();j++)
            outfile<<extra_identifiers.at(j)<<",";
          outfile<<data_filename<<","<<
            response_name<<","<<
              lambda_regularization<<","<<
                sigma_kernel_width<<","<<
                  df<<","<<
                    logli<<","<<
                      BIC<<","<<
                        null_BIC<<","<<
                          return_status<<","<<
                            lgbm_l2_fit<<","<<
                              feature_names.at(i)<<","<<
                                v(i)<<","<<
                                  lgbm_imp(i)<<"\n";
        }
      }
      /////////////////////////////////////////////////////
      void write_to_csv(const std::string& file_location,
                        const std::vector<std::string>& extra_identifiers_header,
                        const std::vector<std::string>& extra_identifiers)
      {
        bool write_header = !io::file_exists_test(file_location);
        std::ofstream outfile(file_location,std::ios_base::app);
        if(write_header){
          write_to_csv_header(outfile,extra_identifiers_header);
        }
        write_to_csv(outfile,extra_identifiers);
        outfile.close();
      }
    };
  }
}
#endif