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
#include <gtest/gtest.h>
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
class FitMathMedTest : public ::testing::Test{
protected:
  FitMathMedTest(){
    data_matrix.resize(25,6);
    data_matrix << 9,7,4,4,4,2,
                   8,7,5,5,3,1,
                   4,8,2,8,8,-4,
                   8,9,1,1,8,-1,
                   2,8,5,3,4,-6,
                   1,3,9,3,9,-2,
                   8,7,6,2,0,1,
                   1,5,6,6,1,-4,
                   5,6,1,7,3,-1,
                   0,7,7,3,2,-7,
                   2,4,0,6,0,-2,
                   6,0,8,0,7,6,
                   3,4,8,4,9,-1,
                   9,7,1,0,2,2,
                   3,7,2,4,3,-4,
                   7,3,8,9,9,4,
                   8,2,3,4,3,6,
                   6,2,7,1,2,4,
                   7,8,0,9,4,-1,
                   6,3,5,4,6,3,
                   3,1,0,2,9,2,
                   4,1,3,8,7,3,
                   3,8,0,3,4,-5,
                   3,5,3,7,8,-2,
                   9,0,4,4,2,9;
    predictors_use_index.push_back(0);
    predictors_use_index.push_back(1);
    predictors_use_index.push_back(2);
    predictors_use_index.push_back(3);
    predictors_use_index.push_back(4);
    which_response = 5;
  }
  Eigen::MatrixXd data_matrix;
  std::vector<unsigned> predictors_use_index;
  unsigned which_response;
  deker::fit::Opt_Param control;
};
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
TEST_F(FitMathMedTest, CheckSplitData){
  deker::fit::Split_Data split_data_test(data_matrix,which_response,predictors_use_index);
  //
  double response_diff_check = (split_data_test.response-data_matrix.col(5)).array().abs().sum();
  ASSERT_NEAR(response_diff_check,0,1e-10)<<"response does not match imputed data";
  //
  for(unsigned i = 0; i<split_data_test.predictors_index.size();i++){
    double predictor_diff_check = (split_data_test.predictors.col(split_data_test.predictors_index.at(i))-data_matrix.col(split_data_test.predictors_index.at(i))).array().abs().sum();
    ASSERT_NEAR(predictor_diff_check,0,1e-10)<<"predictor "<<i<<" does not match imputed data";
  }
}
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
TEST_F(FitMathMedTest, OptSigmaTest){
  std::cerr<<"Blorpglob\n";
  deker::fit::Solve_deker deker_test(data_matrix,predictors_use_index,which_response,control);
  std::cerr<<"Globblorp\n";
  //
  std::vector<unsigned> feature_index = deker_test.get_input_data().predictors_index;
  Eigen::VectorXd v = Eigen::VectorXd::Constant(feature_index.size(),1);
  //
  std::cerr<<"AA\n";
  double sigma = deker_test.opt_sigma(v,feature_index);
  double df,logli;
  std::cerr<<"checkechk\n";
  std::tie(df,logli) = deker_test.calc_model_fit(sigma,v,feature_index);
  std::cerr<<"sigma "<<sigma<<"\n";
  std::cerr<<"df "<<df<<"\n";
  std::cerr<<"logli "<<logli<<"\n";
  //ASSERT_NEAR(sigma,1.54834,1e-1)<<"sigma does not match hand-calculated value";
  //ASSERT_NEAR(AIC,23.0512542,1e-1)<<"BIC does not match hand-calculated value";
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
TEST_F(FitMathMedTest, OptVTest){
  deker::fit::Solve_deker deker_test(data_matrix,predictors_use_index,which_response,control);
  //
  std::vector<unsigned> feature_index = deker_test.get_input_data().predictors_index;
  Eigen::VectorXd v = Eigen::VectorXd::Constant(feature_index.size(),1);
  //
  double current_sigma = 1.54834;
  double lambda = .01;
  Eigen::VectorXd current_v = deker_test.opt_v(lambda,current_sigma,v,feature_index);
  //
  std::cerr<<"v "<<current_v.transpose()<<"\n";
  double df,logli;
  std::tie(df,logli) = deker_test.calc_model_fit(current_sigma,current_v,feature_index);
  std::cerr<<"df "<<df<<"\n";
  std::cerr<<"logli "<<logli<<"\n";
  //ASSERT_NEAR(current_sigma,.885,5e-3);
  //ASSERT_NEAR(AIC,24,1e-1);
  //ASSERT_NEAR(current_v(0),.853,5e-3);
  //ASSERT_NEAR(current_v(1),.635,5e-3);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
TEST_F(FitMathMedTest, SteffIterVTest){
  deker::fit::Solve_deker deker_test(data_matrix,predictors_use_index,which_response,control);
  //
  deker::fit::Output_deker sol_test;
  sol_test.lambda_regularization = .01;
  sol_test.sigma_kernel_width = 1.54834;
  sol_test.feature_index = deker_test.get_input_data().predictors_index;
  sol_test.v = Eigen::VectorXd::Constant(sol_test.feature_index.size(),1);
  //
  deker_test.steff_iter_v(sol_test.lambda_regularization,sol_test);
  //
  std::cerr<<"sigma "<<sol_test.sigma_kernel_width<<"\n";
  std::cerr<<"v "<<sol_test.v.transpose()<<"\n";
  std::cerr<<"return "<<sol_test.return_status<<"\n";
  //
  double df,logli;
  std::tie(df,logli) = deker_test.calc_model_fit(sol_test.sigma_kernel_width,
           sol_test.v,
           sol_test.feature_index);
  std::cerr<<"df "<<df<<"\n";
  std::cerr<<"logli "<<logli<<"\n";
  //ASSERT_NEAR(current_sigma,.885,5e-3);
  //ASSERT_NEAR(AIC,24,1e-1);
  //ASSERT_NEAR(current_v(0),.853,5e-3);
  //ASSERT_NEAR(current_v(1),.635,5e-3);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
TEST_F(FitMathMedTest, InnerFullOptTest){
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
TEST_F(FitMathMedTest, InnerFullOptTestV2){
  deker::fit::Solve_deker deker_test(data_matrix,predictors_use_index,which_response,control);
  //
  deker::fit::Output_deker sol_test;
  double lambda = .5;
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
TEST_F(FitMathMedTest, OptLambdaTest){
  control.lambda_init_max=1;
  control.lambda_bo_maxiter = 50;
  deker::fit::Solve_deker deker_test(data_matrix,predictors_use_index,which_response,control);
  //
  using Lambda_Opt_t = deker::fit::deker_Lambda_Optimizer<deker::fit::Solve_deker,&deker::fit::Solve_deker::inner_full_opt>;
  Lambda_Opt_t lambda_opt;
  lambda_opt.build(control);
  //
  while(lambda_opt.opt_iteration(deker_test)){};
  //
  std::cerr<<"lambda "<<lambda_opt.get_best_sol().lambda_regularization<<"\n";
  std::cerr<<"sigma "<<lambda_opt.get_best_sol().sigma_kernel_width<<"\n";
  std::cerr<<"v "<<lambda_opt.get_best_sol().v.transpose()<<"\n";
  std::cerr<<"return "<<lambda_opt.get_best_sol().return_status<<"\n";
  std::cerr<<"BIC "<<lambda_opt.get_best_sol().BIC<<"\n";
  //
}