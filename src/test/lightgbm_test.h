#include <deker/lgbm_wrapper.h>
#include <deker/fit.h>
#include <deker/io/misc.h>
#include <gtest/gtest.h>
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
TEST (LightGBMTest, LightGBMTestOne){
    Eigen::MatrixXd data_matrix;
    data_matrix.resize(20,2);
    data_matrix << 7,8,
                   5,10,
                   6,1,
                   5,2,
                   9,10,
                   10,1,
                   3,3,
                   2,2,
                   3,7,
                   1,1,
                   9,4,
                   3,4,
                   5,4,
                   9,5,
                   7,7,
                   4,1,
                   4,2,
                   3,2,
                   1,6,
                   4,2;
    Eigen::VectorXd response;
    response.resize(20);
    response << -1,-5,5,3,-1,9,0,0,-4,0,5,-1,1,4,0,3,2,1,-5,2;
    ///
    //char *parameters = "num_threads=1 objective=regression min_data_in_leaf=1 force_row_wise=true min_gain_to_split=.01 metric=l1,l2,rmse";
    char *parameters = "num_threads=1 objective=regression min_data_in_leaf=5 force_col_wise=true min_gain_to_split=.01 boosting=rf bagging_fraction=.666 bagging_freq=1 verbosity=-1";
    deker::LGBM_wrapper booster(data_matrix,response,parameters);
    booster.build_booster();
    booster.run(10);
    std::cerr<<"imp split "<<booster.get_output().feature_imp_split.transpose()<<"\n";
    std::cerr<<"imp gain "<<booster.get_output().feature_imp_gain.transpose()<<"\n";
    std::cerr<<"fit_metric "<<booster.get_output().fit_metric<<"\n";
    //
    std::string write_test_loc = "lgbm_write_test";
    deker::io::write_to_binary<deker::LGBM_wrapper>(write_test_loc,booster);
    deker::LGBM_wrapper booster_2(data_matrix,response,parameters);
    deker::io::read_from_binary<deker::LGBM_wrapper>(write_test_loc,booster_2);
    std::cerr<<"imp split "<<booster_2.get_output().feature_imp_split.transpose()<<"\n";
    std::cerr<<"imp gain "<<booster_2.get_output().feature_imp_gain.transpose()<<"\n";
    std::cerr<<"fit_metric "<<booster_2.get_output().fit_metric<<"\n";
    remove(write_test_loc.c_str());
}
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
TEST (LightGBMTest, LightGBMTestTwo){
    Eigen::MatrixXd data_matrix;
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
    Eigen::VectorXd response;
    response.resize(20);
    response = data_matrix.col(0);
    ///
    //char *parameters = "num_threads=1 objective=regression min_data_in_leaf=1 force_row_wise=true min_gain_to_split=.01 metric=l1,l2,rmse";
    char *parameters = "num_threads=1 objective=regression min_data_in_leaf=5 force_col_wise=true min_gain_to_split=.01 boosting=rf bagging_fraction=.666 bagging_freq=1 verbosity=-1";
    deker::LGBM_wrapper booster(data_matrix.rightCols(5),response,parameters);
    booster.build_booster();
    booster.run(10);
    std::cerr<<"imp split "<<booster.get_output().feature_imp_split.transpose()<<"\n";
    std::cerr<<"imp gain "<<booster.get_output().feature_imp_gain.transpose()<<"\n";
    std::cerr<<"fit_metric "<<booster.get_output().fit_metric<<"\n";
}
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
// TEST (LightGBMTest, LightGBMTestThree){
//     
//     std::string data_file = "/SFS/user/ctc/hayesse/deker_test_data/dream4_data_1.csv.bcsv";
//     unsigned which_response = 1;
//     
//     deker::io::Eigen_Matrix_IO input_data;
//     deker::io::read_from_binary<deker::io::Eigen_Matrix_IO>(data_file,input_data);
//     
//     ////////////////////////////////////////////////
//     //build Solve_deker
//     deker::fit::Solve_deker sol_deker(input_data.matrix,which_response);
//     
//     //
//     std::vector<unsigned> use_predictors;
//     if(use_predictors.size()==0){
//         for(unsigned i = 0; i<input_data.cols; i++)
//             use_predictors.push_back(i);
//     }
//     
//     deker::fit::Opt_Param params;
//     
//     sol_deker.build(use_predictors,
//                     params);
//     
//     ////////////////////////////////////////////////
//     //run LGBM
//     char *parameters = "num_threads=1 objective=regression min_data_in_leaf=5 force_col_wise=true min_gain_to_split=.01 boosting=rf bagging_fraction=.666 bagging_freq=1 verbosity=-1";
//     //char *parameters = "num_threads=1 objective=rmse min_data_in_leaf=10 force_col_wise=true min_gain_to_split=.01 metric=rmse verbosity=-1";
//     Eigen::VectorXd dummy_response = sol_deker.get_input_data().response.array();
//     deker::LGBM_wrapper booster(sol_deker.get_input_data().predictors,dummy_response,parameters);
//     
//     booster.build_booster();
//     booster.run(1000);
//     
//     std::cerr<<"imp split "<<booster.get_feature_imp_split().transpose().array()/1000<<"\n";
//     std::cerr<<"imp gain "<<booster.get_feature_imp_gain().transpose().array().sqrt()/1000<<"\n";
//     std::cerr<<"booster_fit_metric "<<booster.get_fit_metric()<<"\n";
// }
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////