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
#include <gtest/gtest.h>
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
class FitMathSmallTest : public ::testing::Test{
protected:
  FitMathSmallTest(){
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
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
TEST_F(FitMathSmallTest, CheckSplitData){
  deker::fit::Split_Data split_data_test(data_matrix,predictors_use_index,which_response);
  //
  double response_diff_check = (split_data_test.response-data_matrix.col(0)).array().abs().sum();
  ASSERT_NEAR(response_diff_check,0,1e-10)<<"response does not match imputed data";
  //
  for(unsigned i = 0; i<split_data_test.predictors_index.size();i++){
    double predictor_diff_check = (split_data_test.predictors.col(split_data_test.predictors_index.at(i))-data_matrix.col(split_data_test.predictors_index.at(i))).array().abs().sum();
    ASSERT_NEAR(predictor_diff_check,0,1e-10)<<"predictor "<<i<<" does not match imputed data";
  }
}
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
TEST_F(FitMathSmallTest, CheckResponseLabelsCalcs){
  deker::fit::Split_Data split_data_test(data_matrix,predictors_use_index,which_response);
  deker::fit::Response_Labels response_labels_test(split_data_test.response,
                                                   split_data_test.response_sort_index);
  
  Eigen::MatrixXd y_lab_true(20,9);
  y_lab_true << 0,0,0,0,0,0,0,0,0,
                0,0,0,0,0,0,0,0,0,
                1,0,0,0,0,0,0,0,0,
                1,1,0,0,0,0,0,0,0,
                1,1,0,0,0,0,0,0,0,
                1,1,0,0,0,0,0,0,0,
                1,1,1,0,0,0,0,0,0,
                1,1,1,0,0,0,0,0,0,
                1,1,1,0,0,0,0,0,0,
                1,1,1,0,0,0,0,0,0,
                1,1,1,1,0,0,0,0,0,
                1,1,1,1,0,0,0,0,0,
                1,1,1,1,1,0,0,0,0,
                1,1,1,1,1,0,0,0,0,
                1,1,1,1,1,1,0,0,0,
                1,1,1,1,1,1,0,0,0,
                1,1,1,1,1,1,1,0,0,
                1,1,1,1,1,1,1,1,0,
                1,1,1,1,1,1,1,1,0,
                1,1,1,1,1,1,1,1,1;
  double y_lab_check = (y_lab_true-response_labels_test.y_labels).array().abs().sum();
  ASSERT_NEAR(y_lab_check,0,1e-10)<<"y_labels does not match hand-calculated matrix";
  Eigen::MatrixXd y_pseudo_true(9,20);
  y_pseudo_true << 0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                   0,0,-1,.333,.333,.333,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                   0,0,0,-.333,-.333,-.333,.25,.25,.25,.25,0,0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,-.25,-.25,-.25,-.25,.5,.5,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,0,0,-.5,-.5,.5,.5,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,0,0,0,0,-.5,-.5,.5,.5,0,0,0,0,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,-.5,-.5,1,0,0,0,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,.5,.5,0,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-.5,-.5,1;
  double y_pseudo_check = (y_pseudo_true-response_labels_test.y_labels_pseudoinv).array().abs().sum();
  ASSERT_NEAR(y_pseudo_check,0,1e-2)<<"y_pseudoinv does not match hand-calculated matrix";
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
TEST_F(FitMathSmallTest, OptSigmaTest){
  deker::fit::Solve_deker deker_test(data_matrix,predictors_use_index,which_response,control);
  //
  Eigen::VectorXd v = Eigen::VectorXd::Constant(2,1);
  std::vector<unsigned> feature_index = deker_test.get_input_data().predictors_index;
  //
  double sigma = deker_test.opt_sigma(v,feature_index);
  double df,logli,A,B;
  std::tie(df,logli,A,B) = deker_test.calc_model_fit(sigma,v,feature_index);
  double BIC = 2*df-2*logli;
  std::cerr<<"sigma "<<sigma<<"\n";
  std::cerr<<"df "<<df<<"\n";
  std::cerr<<"logli "<<logli<<"\n";
  std::cerr<<"A "<<A<<"\n";
  std::cerr<<"B "<<B<<"\n";
  //ASSERT_NEAR(sigma,1.54834,1e-1)<<"sigma does not match hand-calculated value";
  //ASSERT_NEAR(AIC,23.0512542,1e-1)<<"BIC does not match hand-calculated value";
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
TEST_F(FitMathSmallTest, OptVTest){
  deker::fit::Solve_deker deker_test(data_matrix,predictors_use_index,which_response,control);
  //
  Eigen::VectorXd v = Eigen::VectorXd::Constant(2,1);
  std::vector<unsigned> feature_index = deker_test.get_input_data().predictors_index;
  //
  double current_sigma = 1.54834;
  double lambda = .01;
  Eigen::VectorXd current_v = deker_test.opt_v(lambda,current_sigma,v,feature_index);
  //
  std::cerr<<"v "<<current_v.transpose()<<"\n";
  double df,logli,A,B;
  std::tie(df,logli,A,B) = deker_test.calc_model_fit(current_sigma,current_v,feature_index);
  std::cerr<<"df "<<df<<"\n";
  std::cerr<<"logli "<<logli<<"\n";
  std::cerr<<"A "<<A<<"\n";
  std::cerr<<"B "<<B<<"\n";
  //ASSERT_NEAR(current_sigma,.885,5e-3);
  //ASSERT_NEAR(AIC,24,1e-1);
  //ASSERT_NEAR(current_v(0),.853,5e-3);
  //ASSERT_NEAR(current_v(1),.635,5e-3);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
TEST_F(FitMathSmallTest, SteffIterVTest){
  deker::fit::Solve_deker deker_test(data_matrix,predictors_use_index,which_response,control);
  //
  deker::fit::Output_deker sol_test;
  sol_test.lambda_regularization = .01;
  sol_test.sigma_kernel_width = 1.54834;
  sol_test.v = Eigen::VectorXd::Constant(2,1);
  sol_test.feature_index = deker_test.get_input_data().predictors_index;
  //
  deker_test.steff_iter_v(sol_test.lambda_regularization,sol_test);
  //
  std::cerr<<"sigma "<<sol_test.sigma_kernel_width<<"\n";
  std::cerr<<"v "<<sol_test.v.transpose()<<"\n";
  std::cerr<<"return "<<sol_test.return_status<<"\n";
  //
  double df,logli,A,B;
  std::tie(df,logli,A,B) = deker_test.calc_model_fit(sol_test.sigma_kernel_width,
                                                     sol_test.v,
                                                     sol_test.feature_index);
  std::cerr<<"df "<<df<<"\n";
  std::cerr<<"logli "<<logli<<"\n";
  std::cerr<<"A "<<A<<"\n";
  std::cerr<<"B "<<B<<"\n";
  //ASSERT_NEAR(current_sigma,.885,5e-3);
  //ASSERT_NEAR(AIC,24,1e-1);
  //ASSERT_NEAR(current_v(0),.853,5e-3);
  //ASSERT_NEAR(current_v(1),.635,5e-3);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
TEST_F(FitMathSmallTest, InnerFullOptTest){
  deker::fit::Solve_deker deker_test(data_matrix,predictors_use_index,which_response,control);
  //
  deker::fit::Output_deker sol_test;
  double lambda = .01;
  //
  deker_test.inner_full_opt(lambda,sol_test);
  //
  std::cerr<<"sigma "<<sol_test.sigma_kernel_width<<"\n";
  std::cerr<<"v "<<sol_test.v.transpose()<<"\n";
  std::cerr<<"return "<<sol_test.return_status<<"\n";
  std::cerr<<"BIC "<<sol_test.BIC<<"\n";
  //
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
TEST_F(FitMathSmallTest, OptLambdaTest){
  deker::fit::Solve_deker deker_test(data_matrix,predictors_use_index,which_response,control);
  //
  deker::fit::Output_deker sol_test;
  //
  deker_test.opt_lambda(sol_test);
  //
  std::cerr<<"lambda "<<sol_test.lambda_regularization<<"\n";
  std::cerr<<"sigma "<<sol_test.sigma_kernel_width<<"\n";
  std::cerr<<"v "<<sol_test.v.transpose()<<"\n";
  std::cerr<<"return "<<sol_test.return_status<<"\n";
  std::cerr<<"BIC "<<sol_test.BIC<<"\n";
  //
}