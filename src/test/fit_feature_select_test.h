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
#include <deker/fit/feature_select.h>
#include <gtest/gtest.h>
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
class FitFeatureSelectTest : public ::testing::Test{
protected:
  FitFeatureSelectTest(){
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    input_matrix.resize(20,7);
    input_matrix << -1,7,8,9,4,8,3,
                    -5,5,10,8,6,6,10,
                    5,6,1,9,3,2,8,
                    3,5,2,9,2,5,4,
                    -1,9,10,6,10,5,2,
                    9,10,1,9,3,9,6,
                    0,3,3,8,5,5,10,
                    0,2,2,5,7,5,9,
                    -4,3,7,10,8,2,10,
                    0,1,1,9,8,10,6,
                    5,9,4,5,2,5,9,
                    -1,3,4,10,2,7,9,
                    1,5,4,1,2,2,9,
                    4,9,5,2,5,4,9,
                    0,7,7,8,9,10,8,
                    3,4,1,9,10,6,3,
                    2,4,2,3,2,2,4,
                    1,3,2,1,1,2,8,
                    -5,1,6,5,8,9,3,
                    2,4,2,9,3,9,5;
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    use_predictors.push_back(0);
    use_predictors.push_back(1);
    use_predictors.push_back(2);
    use_predictors.push_back(3);
    use_predictors.push_back(4);
    use_predictors.push_back(5);
    use_predictors.push_back(6);
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    which_response = 0;
  }
  double sigma = 0;
  double lambda = .01;
  Eigen::MatrixXd input_matrix;
  std::vector<unsigned> use_predictors;
  unsigned which_response;
  deker::fit::Opt_Param control;
};
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
TEST_F(FitFeatureSelectTest, DekerTest){
  //deker::fit::Split_Data input_data(input_matrix,use_predictors,which_response);
  deker::fit::Deker deker_test(input_matrix,use_predictors,which_response,sigma,lambda,control,false);
  //
  ASSERT_NEAR(deker_test.get_sigma(),.885,5e-3);
  Eigen::VectorXd check_w = Eigen::VectorXd::Zero(7);
  for(unsigned i = 0; i<deker_test.get_feature_index().size();i++){
    check_w(deker_test.get_feature_index().at(i)) = deker_test.get_w()(i);
  }
  Eigen::VectorXd true_w(7);
  true_w<<0,.853,.635,0,0,0,0;
  true_w = true_w.array().square();
  ASSERT_NEAR((check_w-true_w).array().abs().maxCoeff(),0,5e-3);
}