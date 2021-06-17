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
      double BIC;
      int return_status;
      Eigen::VectorXd v;
      std::vector<unsigned> feature_index;
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
        BIC = new_output.BIC;
        return_status = new_output.return_status;
        v = new_output.v;
        feature_index = new_output.feature_index;
      }
    };
  }
}
#endif