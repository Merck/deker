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
#ifndef DEKER_FIT_OUTPUT
#define DEKER_FIT_OUTPUT
///////////////////////////////////////////
#include <Eigen/Core>
#include <vector>
///////////////////////////////////////////
namespace deker{
  namespace fit{
    struct Output_deker{
      double lambda_regularization;
      double sigma_kernel_width;
      double df;
      double logli;
      double BIC;
      int return_status;
      Eigen::VectorXd v;
      std::vector<unsigned> feature_index;
      std::vector<unsigned> key_to_original_data;
      ////
      Output_deker()
      {
        return_status=-1;
        lambda_regularization=0;
        sigma_kernel_width=5;
      }
      ////
      void operator()(const Output_deker& new_output)
      {
        lambda_regularization = new_output.lambda_regularization;
        sigma_kernel_width = new_output.sigma_kernel_width;
        df = new_output.df;
        logli = new_output.logli;
        BIC = new_output.BIC;
        return_status = new_output.return_status;
        v = new_output.v;
        feature_index = new_output.feature_index;
        key_to_original_data = new_output.key_to_original_data;
      }
      ////////////////////////////////////////////////////
      void serialize(std::ifstream& infile)
      {
        feature_index.clear();
        key_to_original_data.clear();
        infile.read((char*) (&lambda_regularization),sizeof(lambda_regularization));
        infile.read((char*) (&sigma_kernel_width),sizeof(sigma_kernel_width));
        infile.read((char*) (&df),sizeof(df));
        infile.read((char*) (&logli),sizeof(logli));
        infile.read((char*) (&BIC),sizeof(BIC));
        infile.read((char*) (&return_status),sizeof(return_status));
        size_t n_features;
        infile.read((char*) (&n_features),sizeof(n_features));
        v.resize(n_features);
        unsigned tmp_index;
        double tmp_v;
        for(size_t i = 0; i<n_features;i++){
          infile.read((char*) (&tmp_index), sizeof(unsigned));
          feature_index.push_back(tmp_index);
          infile.read((char*) (&tmp_v), sizeof(double));
          v(i) = tmp_v;
        }
        size_t tmp_size;
        infile.read((char*) (&tmp_size),sizeof(tmp_size));
        for(size_t i = 0; i<tmp_size;i++){
          infile.read((char*) (&tmp_index), sizeof(tmp_index));
          key_to_original_data.push_back(tmp_index);
        }
      }
      ////////////////////////////////////////////////////
      void serialize(std::ofstream& outfile) const
      {
        outfile.write((char*) (&lambda_regularization), sizeof(lambda_regularization));
        outfile.write((char*) (&sigma_kernel_width), sizeof(sigma_kernel_width));
        outfile.write((char*) (&df), sizeof(df));
        outfile.write((char*) (&logli), sizeof(logli));
        outfile.write((char*) (&BIC), sizeof(BIC));
        outfile.write((char*) (&return_status), sizeof(return_status));
        size_t n_features = v.size();
        outfile.write((char*) (&n_features), sizeof(n_features));
        for(size_t i = 0;i<n_features;i++){
          outfile.write((char*) (&feature_index.at(i)),sizeof(unsigned));
          outfile.write((char*) (&v(i)),sizeof(double));
        }
        size_t tmp_size = key_to_original_data.size();
        outfile.write((char*) (&tmp_size), sizeof(tmp_size));
        for(size_t i = 0;i<tmp_size; i++){
          outfile.write((char*) (&key_to_original_data.at(i)),sizeof(unsigned));
        }
      }
    };
  }
}
#endif