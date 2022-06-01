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
#include <deker/fit/predictions.h>
#include <deker/fit.h>
#include <gtest/gtest.h>
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
class FitPredictionsTest : public ::testing::Test{
protected:
    FitPredictionsTest(){
        data_matrix_vsm.resize(10,2);
        data_matrix_vsm << 0,1,
                           0,2,
                           0,3,
                           0,4,
                           0,5,
                           0,6,
                           1,7,
                           1,8,
                           1,9,
                           1,10;
        predictors_use_index_vsm.push_back(1);
        which_response_vsm = 0;

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
        //
        data_matrix_reg.resize(12,6);
        data_matrix_reg << 9,7,4,4,4,2,
                           8,7,5,5,3,1,
                           4,8,2,8,8,-4,
                           8,9,1,1,8,-1,
                           2,8,5,3,4,-6,
                           1,3,9,3,9,-2,
                           //8,7,6,2,0,1,
                           //1,5,6,6,1,-4,
                           ///5,6,1,7,3,-1,
                           0,7,7,3,2,-7,
                           //2,4,0,6,0,-2,
                           6,0,8,0,7,6,
                           //3,4,8,4,9,-1,
                           //9,7,1,0,2,2,
                           //3,7,2,4,3,-4,
                           7,3,8,9,9,4,
                           //8,2,3,4,3,6,
                           //6,2,7,1,2,4,
                           //7,8,0,9,4,-1,
                           6,3,5,4,6,3,
                           //3,1,0,2,9,2,
                           //4,1,3,8,7,3,
                           3,8,0,3,4,-5,
                           //3,5,3,7,8,-2,
                           9,0,4,4,2,9;
        
        predictors_use_index.push_back(0);
        predictors_use_index.push_back(1);
        predictors_use_index.push_back(2);
        predictors_use_index.push_back(3);
        predictors_use_index.push_back(4);
        which_response = 5;
        //
        data_matrix_sm.resize(20,3);
        data_matrix_sm << -1,7,8,
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
        predictors_use_index_sm.push_back(1);
        predictors_use_index_sm.push_back(2);
        which_response_sm = 0;
    }
    
    deker::fit::Opt_Param control;
    //
    Eigen::MatrixXd data_matrix_vsm;
    std::vector<unsigned> predictors_use_index_vsm;
    unsigned which_response_vsm;
    //
    Eigen::MatrixXd data_matrix;
    Eigen::MatrixXd data_matrix_reg;
    std::vector<unsigned> predictors_use_index;
    unsigned which_response;
    //
    Eigen::MatrixXd data_matrix_sm;
    std::vector<unsigned> predictors_use_index_sm;
    unsigned which_response_sm;
};
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
TEST_F(FitPredictionsTest, PredictTest){
    /////////////////////////////////
    std::cerr<<"#################\n";
    std::cerr<<"#################\n";
    deker::fit::Solve_deker deker_test_vsm(data_matrix_vsm,predictors_use_index_vsm,which_response_vsm,control);
    std::vector<unsigned> feature_index_vsm = deker_test_vsm.get_input_data().predictors_index;
    Eigen::VectorXd v_vsm = Eigen::VectorXd::Constant(feature_index_vsm.size(),1);
    std::cerr<<"#################\n";
    std::cerr<<"#################\n";
    deker::fit::Predictions_deker predict_1_vsm(deker_test_vsm.get_input_data(),
                                                deker_test_vsm.get_response_labels(),
                                                -5,
                                                v_vsm,
                                                feature_index_vsm);
    std::cerr<<"#################\n";
    deker::fit::Predictions_deker predict_2_vsm(deker_test_vsm.get_input_data(),
                                                deker_test_vsm.get_response_labels(),
                                                -2.5,
                                                v_vsm,
                                                feature_index_vsm);
    std::cerr<<"#################\n";
    deker::fit::Predictions_deker predict_3_vsm(deker_test_vsm.get_input_data(),
                                                deker_test_vsm.get_response_labels(),
                                                -1,
                                                v_vsm,
                                                feature_index_vsm);
    std::cerr<<"#################\n";
    deker::fit::Predictions_deker predict_4_vsm(deker_test_vsm.get_input_data(),
                                                deker_test_vsm.get_response_labels(),
                                                -.1,
                                                v_vsm,
                                                feature_index_vsm);
    std::cerr<<"#################\n";
    deker::fit::Predictions_deker predict_5_vsm(deker_test_vsm.get_input_data(),
                                                deker_test_vsm.get_response_labels(),
                                                0,
                                                v_vsm,
                                                feature_index_vsm);
    std::cerr<<"#################\n";
    deker::fit::Predictions_deker predict_6_vsm(deker_test_vsm.get_input_data(),
                                                deker_test_vsm.get_response_labels(),
                                                .1,
                                                v_vsm,
                                                feature_index_vsm);
    std::cerr<<"#################\n";
    deker::fit::Predictions_deker predict_7_vsm(deker_test_vsm.get_input_data(),
                                                deker_test_vsm.get_response_labels(),
                                                1,
                                                v_vsm,
                                                feature_index_vsm);
    std::cerr<<"#################\n";
    deker::fit::Predictions_deker predict_8_vsm(deker_test_vsm.get_input_data(),
                                                deker_test_vsm.get_response_labels(),
                                                2.5,
                                                v_vsm,
                                                feature_index_vsm);
    std::cerr<<"#################\n";
    deker::fit::Predictions_deker predict_9_vsm(deker_test_vsm.get_input_data(),
                                                deker_test_vsm.get_response_labels(),
                                                5,
                                                v_vsm,
                                                feature_index_vsm);
    std::cerr<<"#################\n";
    /////////////////////////////////
    std::cerr<<"#################\n";
    deker::fit::Solve_deker deker_test_sm(data_matrix_sm,predictors_use_index_sm,which_response_sm,control);
    //
    std::vector<unsigned> feature_index_sm = deker_test_sm.get_input_data().predictors_index;
    Eigen::VectorXd v_sm = Eigen::VectorXd::Constant(feature_index_sm.size(),1);

    deker::fit::Predictions_deker predict_1_sm(deker_test_sm.get_input_data(),
                                          deker_test_sm.get_response_labels(),
                                          -1,
                                          v_sm,
                                          feature_index_sm);
    std::cerr<<"#################\n";
    deker::fit::Predictions_deker predict_2_sm(deker_test_sm.get_input_data(),
                                          deker_test_sm.get_response_labels(),
                                          0,
                                          v_sm,
                                          feature_index_sm);
    std::cerr<<"#################\n";
    deker::fit::Predictions_deker predict_3_sm(deker_test_sm.get_input_data(),
                                          deker_test_sm.get_response_labels(),
                                          1,
                                          v_sm,
                                          feature_index_sm);
    std::cerr<<"#################\n";
    deker::fit::Predictions_deker predict_4_sm(deker_test_sm.get_input_data(),
                                          deker_test_sm.get_response_labels(),
                                          10,
                                          v_sm,
                                          feature_index_sm);
    std::cerr<<"#################\n";
    std::cerr<<"#################\n";
    deker::fit::Solve_deker deker_test_med(data_matrix,predictors_use_index,which_response,control);
    std::cerr<<"#################\n";
    std::cerr<<"#################\n";
    std::vector<unsigned> feature_index_med = deker_test_med.get_input_data().predictors_index;
    Eigen::VectorXd v_med = Eigen::VectorXd::Constant(feature_index_med.size(),1);

    deker::fit::Predictions_deker predict_1_med(deker_test_med.get_input_data(),
                                            deker_test_med.get_response_labels(),
                                            -1,
                                            v_med,
                                            feature_index_med);
    std::cerr<<"#################\n";
    deker::fit::Predictions_deker predict_2_med(deker_test_med.get_input_data(),
                                            deker_test_med.get_response_labels(),
                                            0,
                                            v_med,
                                            feature_index_med);
    std::cerr<<"#################\n";
    deker::fit::Predictions_deker predict_3_med(deker_test_med.get_input_data(),
                                            deker_test_med.get_response_labels(),
                                            1,
                                            v_med,
                                            feature_index_med);
    std::cerr<<"#################\n";
    deker::fit::Predictions_deker predict_4_med(deker_test_med.get_input_data(),
                                            deker_test_med.get_response_labels(),
                                            10,
                                            v_med,
                                            feature_index_med);
    std::cerr<<"#################\n";
    std::cerr<<"#################\n";
    std::cerr<<"#################\n";
    deker::fit::Solve_deker deker_test_reg(data_matrix_reg,predictors_use_index,which_response,control);
    std::cerr<<"#################\n";
    std::cerr<<"#################\n";
    std::vector<unsigned> feature_index_reg = deker_test_reg.get_input_data().predictors_index;
    Eigen::VectorXd v_reg = Eigen::VectorXd::Constant(feature_index_reg.size(),1);

    deker::fit::Predictions_deker predict_1_reg(deker_test_reg.get_input_data(),
                                                deker_test_reg.get_response_labels(),
                                                -1,
                                                v_reg,
                                                feature_index_reg);
    std::cerr<<"#################\n";
    deker::fit::Predictions_deker predict_2_reg(deker_test_reg.get_input_data(),
                                                deker_test_reg.get_response_labels(),
                                                0,
                                                v_reg,
                                                feature_index_reg);
    std::cerr<<"#################\n";
    deker::fit::Predictions_deker predict_3_reg(deker_test_reg.get_input_data(),
                                                deker_test_reg.get_response_labels(),
                                                1,
                                                v_reg,
                                                feature_index_reg);
    std::cerr<<"#################\n";
    deker::fit::Predictions_deker predict_4_reg(deker_test_reg.get_input_data(),
                                                deker_test_reg.get_response_labels(),
                                                10,
                                                v_reg,
                                                feature_index_reg);
    std::cerr<<"#################\n";
}

