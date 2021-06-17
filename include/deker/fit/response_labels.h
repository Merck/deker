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
#ifndef DEKER_FIT_RESPONSE_LABELS
#define DEKER_FIT_RESPONSE_LABELS
///////////////////////////////////////////
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <Eigen/LU>
#include <vector>
///////////////////////////////////////////
namespace deker{
  namespace fit{
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    struct Response_Labels{
      std::vector<unsigned> threshold_inds;
      std::vector<unsigned> threshold_mult;
      std::vector<unsigned> threshold_sorty_map;
      Eigen::MatrixXd y_labels; //samples x samples
      Eigen::MatrixXd y_labels_pseudoinv;
      Eigen::MatrixXd y_labels_pseudoident;
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
    public:
      Response_Labels(){}
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      void build(const Eigen::Ref<const Eigen::VectorXd> response,
                 const std::vector<unsigned>& response_sort_index)
      {
        //find number of classifiers to build (number of unique values in response)
        //ignore the lowest point - will not have two classes to split on
        //threshold_vals.push_back(response(response_sort_index.at(0)));
        threshold_inds.push_back(0);
        threshold_mult.push_back(1);
        threshold_sorty_map.push_back(0); 
        double last_threshold_val = response(response_sort_index.at(0));
        for(unsigned i = 1; i<response_sort_index.size();i++){
          if(response(response_sort_index.at(i))>last_threshold_val){
            last_threshold_val = response(response_sort_index.at(i));
            threshold_inds.push_back(i);
            threshold_mult.push_back(1);
          }else{
            threshold_mult.back()++;
          }
          threshold_sorty_map.push_back(threshold_inds.size()-1); //minimum 0
        }
        //thresholds start with lowest, which will not have a valid split (all >=)
        //so skip lowest split in y_labels
        y_labels = Eigen::MatrixXd::Zero(response_sort_index.size(),threshold_inds.size()-1);
        for(unsigned i = 1;i<threshold_inds.size();i++){
          for(unsigned j = threshold_inds.at(i);j<response_sort_index.size();j++){
            y_labels(j,i-1) = 1.0;
          }
        }
        y_labels_pseudoinv = (y_labels.transpose()*y_labels).colPivHouseholderQr().solve(Eigen::MatrixXd::Identity(y_labels.cols(),y_labels.cols()))*y_labels.transpose();
        y_labels_pseudoident = y_labels*y_labels_pseudoinv;
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      Response_Labels(const Eigen::Ref<const Eigen::VectorXd> response,
                      const std::vector<unsigned>& response_sort_index)
      {
        build(response,response_sort_index);
      }
    };
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
  }
}
#endif