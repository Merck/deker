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
#include <deker/fit.h>
#include <deker/opt_lambda.h>
#include <deker/io/misc.h>
#include <gtest/gtest.h>
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
class HoldoutFitSmallTest : public ::testing::Test{
protected:
  HoldoutFitSmallTest(){
    data_matrix.resize(20,3);
    data_matrix << -1,7,8,
                   -5,5,10,
                   5,6,1,
                   3,5,2,
                   -1,9,10,
                   9,10,1,
                   0,3,3,
                   0,2,2,
                   -4,3,7,
                   0,1,1,
                   5,9,4,
                   -1,3,4,
                   1,5,4,
                   4,9,5,
                   0,7,7,
                   3,4,1,
                   2,4,2,
                   1,3,2,
                   -5,1,6,
                   2,4,2;
    predictors_use_index.push_back(1);
    predictors_use_index.push_back(2);
    which_response = 0;
  }
  Eigen::MatrixXd data_matrix;
  std::vector<unsigned> predictors_use_index;
  unsigned which_response;
  deker::fit::Opt_Param control;
};
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
TEST_F(HoldoutFitSmallTest, DropSampleTest){
  deker::io::Eigen_Matrix_IO holdout_test;
  std::vector<std::string> dummy_names;
  for(unsigned i = 0; i<data_matrix.rows(); i++){
    holdout_test.build(data_matrix,dummy_names);
    std::vector<unsigned> drop_inds(1,i);
    holdout_test.drop_rows(drop_inds);
    for(unsigned j = 0; j<holdout_test.rows; j++){
      if(j == i){
        ASSERT_EQ((data_matrix.row(data_matrix.rows()-1)-holdout_test.matrix.row(j)).sum(),0)<<"holdout matrix is not as expected after dropping row "<<i<<" at new row "<<j<<"\n";
      }else{
        ASSERT_EQ((data_matrix.row(j)-holdout_test.matrix.row(j)).sum(),0)<<"holdout matrix is not as expected after dropping row "<<i<<" at new row "<<j<<"\n";
      }
    }
  }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
TEST_F(HoldoutFitSmallTest, HoldoutPredictionsTest){
  deker::io::Eigen_Matrix_IO holdout_test;
  std::vector<std::string> dummy_names;
  //
  Eigen::VectorXd v(2);
  v<<0.639349, 0.584701;
  std::string write_test_out = "write_holdout_test.csv";
  std::vector<std::string> header;
  header.push_back("hold_out_index");
  std::vector<std::string> header_ind(1);
  std::vector<unsigned> test_samples(1);
  for(unsigned i = 0; i<data_matrix.rows();i++){
    holdout_test.build(data_matrix,dummy_names);
    std::vector<unsigned> drop_inds(1,i);
    header_ind.at(0) = std::to_string(i);
    holdout_test.drop_rows(drop_inds);
    deker::fit::Solve_deker deker_test(holdout_test.matrix,predictors_use_index,which_response,control);
    std::vector<unsigned> feature_index = deker_test.get_input_data().predictors_index;
    Eigen::MatrixXd test_predictors(1,data_matrix.cols());
    test_predictors.row(0) = data_matrix.row(i);
    
    Eigen::VectorXd test_response(1);
    test_response(0) = data_matrix(i,which_response);
    //
    deker::fit::Predictions_deker_holdout predict_out(deker_test.get_input_data(),
                                                      deker_test.get_response_labels(),
                                                      6.56519,
                                                      v,
                                                      feature_index,
                                                      test_predictors,
                                                      test_response);
    predict_out.write_to_csv(write_test_out,header,header_ind);
  }
  remove(write_test_out.c_str());
}