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
      // Eigen::VectorXd v;
      // std::vector<unsigned> feature_index;
      double null_BIC;
      unsigned return_status;
      std::vector<std::string> feature_names;
      unsigned which_response;
      std::string response_name;
      std::string data_filename;
      /////////////////////////////////////
      //constructor using a base Opt_Param
      Output_deker_IO():Output_deker(){}
      Output_deker_IO(const Output_deker& out_deker):Output_deker(out_deker){}
      /////////////////////////////////////
      Output_deker_IO(const Output_deker& out_deker,
                      const std::vector<std::string>& full_feature_names,
                      const std::string& data_filename,
                      const unsigned& which_response):
        Output_deker(out_deker),
        data_filename(data_filename),
        which_response(which_response),
        response_name(full_feature_names.at(which_response))
        {
          for(unsigned i = 0; i<feature_index.size(); i++)
            feature_names.push_back(full_feature_names.at(feature_index.at(i)));
        }
      /////////////////////////////////////
      Output_deker_IO(const std::string& input_file)
      {
        if(io::file_exists_test(input_file)&&!input_file.empty()){
          std::ifstream infile(input_file, std::ios::in | std::ios::binary);
          serialize(infile);
          infile.close();
        }
        else{
          throw std::invalid_argument("input file "+input_file+" not found");
        }
      }
      ////////////////////////////////////////////////////
      void serialize(std::ifstream& infile)
      {
        //read in from stream
        infile.read((char*) (&lambda_regularization),sizeof(lambda_regularization));
        infile.read((char*) (&sigma_kernel_width),sizeof(sigma_kernel_width));
        infile.read((char*) (&BIC),sizeof(BIC));
        infile.read((char*) (&return_status),sizeof(return_status));
        infile.read((char*) (&null_BIC),sizeof(null_BIC));
        size_t n_features;
        infile.read((char*) (&n_features),sizeof(n_features));
        v.resize(n_features);
        unsigned temp_index;
        size_t temp_name_size;
        std::string temp_name;
        double temp_v;
        for(unsigned i = 0; i<n_features;i++){
          infile.read((char*) (&temp_index), sizeof(unsigned));
          feature_index.push_back(temp_index);
          infile.read((char*) (&temp_name_size), sizeof(temp_name_size));
          temp_name.resize(temp_name_size);
          infile.read(&temp_name[0], temp_name_size);
          feature_names.push_back(temp_name);
          infile.read((char*) (&temp_v), sizeof(double));
          v(i) = temp_v;
        }
        infile.read((char*) (&temp_name_size), sizeof(temp_name_size));
        data_filename.resize(temp_name_size);
        infile.read(&data_filename[0],temp_name_size);
        infile.read((char*) (&which_response), sizeof(which_response));
        infile.read((char*) (&temp_name_size), sizeof(temp_name_size));
        response_name.resize(temp_name_size);
        infile.read(&response_name[0],temp_name_size);
      }
      ////////////////////////////////////////////////////
      void serialize(std::ofstream& outfile)
      {
        outfile.write((char*) (&lambda_regularization), sizeof(lambda_regularization));
        outfile.write((char*) (&sigma_kernel_width), sizeof(sigma_kernel_width));
        outfile.write((char*) (&BIC), sizeof(BIC));
        outfile.write((char*) (&return_status), sizeof(return_status));
        outfile.write((char*) (&null_BIC), sizeof(null_BIC));
        size_t n_features = v.size();
        outfile.write((char*) (&n_features), sizeof(n_features));
        for(unsigned i = 0;i<n_features;i++){
          outfile.write((char*) (&feature_index.at(i)),sizeof(unsigned));
          size_t size=feature_names.at(i).size();
          outfile.write((char*) (&size),sizeof(size));
          outfile.write(feature_names.at(i).c_str(), size);
          outfile.write((char*) (&v(i)),sizeof(double));
        }
        size_t data_filename_size = data_filename.size();
        outfile.write((char*) (&data_filename_size),sizeof(data_filename_size));
        outfile.write(data_filename.c_str(),data_filename_size);
        outfile.write((char*) (&which_response),sizeof(which_response));
        size_t response_name_size = response_name.size();
        outfile.write((char*) (&response_name_size),sizeof(response_name_size));
        outfile.write(response_name.c_str(),response_name_size);
      }
      /////////////////////////////////////////////////////
      void write_binary(const std::string& output_file)
      {
        std::ofstream outfile(output_file, std::ios::out | std::ios::binary | std::ios::trunc);
        serialize(outfile);
        outfile.close();
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
                "BIC,"<<
                  "null_BIC,"<<
                    "return_status,"<<
                      "feature,"<<
                        "weight\n";
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
                  BIC<<","<<
                    null_BIC<<","<<
                      return_status<<","<<
                        feature_names.at(i)<<","<<
                          v(i)<<"\n";
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