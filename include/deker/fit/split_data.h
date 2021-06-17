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
#ifndef DEKER_FIT_SPLIT_DATA
#define DEKER_FIT_SPLIT_DATA
///////////////////////////////////////////
#include <Eigen/Core>
#include <vector>
///////////////////////////////////////////
namespace deker{
  namespace fit{
  /////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////
  //stolen from stackoverflow, credit Lukasz Wiklendt
  std::vector<unsigned> sort_indexes(const std::vector<double>& v) {
    // initialize original index locations
    std::vector<unsigned> idx(v.size());
    for(unsigned i = 0; i<idx.size();i++) idx.at(i) = i;
    // sort indexes based on comparing values in v
    sort(idx.begin(), idx.end(), [&v](unsigned i1, unsigned i2) {return v.at(i1) < v.at(i2);});
    return idx;
  }
  /////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////
  //stolen from stackoverflow, credit Lukasz Wiklendt
  std::vector<unsigned> sort_indexes(const Eigen::Ref<const Eigen::VectorXd> v) {
    // initialize original index locations
    std::vector<unsigned> idx(v.size());
    for(unsigned i = 0; i<idx.size();i++) idx.at(i) = i;
    // sort indexes based on comparing values in v
    sort(idx.begin(), idx.end(), [&v](unsigned i1, unsigned i2) {return v(i1) < v(i2);});
    return idx;
  }
  /////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////
  struct Split_Data{
    const Eigen::Ref<const Eigen::MatrixXd> predictors;
    Eigen::VectorXd response;
    std::vector<unsigned> response_sort_index;
    std::vector<unsigned> predictors_index;
    ///////////////////////////////////////////
    Split_Data(const Eigen::Ref<const Eigen::MatrixXd> input_data,
               const std::vector<unsigned>& predictors_use_index,
               const unsigned& which_response):
      predictors(input_data),
      predictors_index(predictors_use_index)
    {
      //extract response 
      response = predictors.col(which_response);
      for(unsigned i = 0; i<predictors_index.size();i++)
        if(which_response==predictors_index.at(i)){
          predictors_index.at(i) = predictors_index.back();
          predictors_index.pop_back();
          break;
        }
        response_sort_index = sort_indexes(response);
    }
  };
  /////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////
  }
}
#endif