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
#include <deker/lgbm_wrapper.h>
#include <deker/opt_lambda.h>
#include <deker/io/misc.h>
#include <deker/io/digest.h>
#include <deker/io/output.h>
//
#include <experimental/filesystem>
#include <unistd.h>
#include <time.h>
////////////////////////////////////////////////////////////////////////////////////////////////
//main program control
int main(int argc,  char **argv){
  //setup for options
  bool all_out = false; //a
  bool binary_convert = false; //b
  bool do_collect_check = false; //c
  bool do_digest_check = false; //d
  bool make_edge_list = false; //e
  bool do_fit = false; //f
  bool do_holdout_predict = false; //h
  bool do_manual_fit = false; //m
  bool do_digest_count = false; //n
  bool predict_out = false; //p
  bool start_control = false; //s
  /////////////////////////////////////
  /////////////////////////////////////
  //parse options
  int opt;
  while((opt = getopt(argc, argv, "abcdefhmnps")) != -1){  
    switch(opt){
    case 'a':
      all_out = true; //for -e and -p, output data for all responses for a model
      break;
    case 'b': //switch to use deker_job to convert a csv file to binary
      binary_convert = true;
      break;
    case 'c': //check if fit for a response is done or errored
      do_collect_check = true;
      break;
    case 'd': //check digest contents
        do_digest_check = true;
        break;
    case 'e': //switch for creating edgelist output
      make_edge_list = true;
      break;
    case 'f': //switch for doing fit 
      do_fit = true;
      break;
    case 'h': //switch for doing hold out with predict-out from fit object
      do_holdout_predict = true;
      break;
    case 'm': //switch for doing manual fit 
      do_manual_fit = true;
      break;
    case 'n': //switch for doing digest count
      do_digest_count = true;
      break;
    case 'p': //switch for doing predict-out from fit object
      predict_out = true;
      break;
    case 's': //switch for starting opt control
      start_control = true;
      break;
    }  
  }
  /////////////////////////////////////
  /////////////////////////////////////
  if(binary_convert){
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //-b
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(argc-optind<2)
      throw std::invalid_argument("2 non-option arguments are required: data file location & location to write binary file\n");
    deker::io::Eigen_Matrix_IO to_convert(argv[optind],false);
    deker::io::write_to_binary<deker::io::Eigen_Matrix_IO>(argv[optind+1],to_convert);
    
    /////////////////////////////////////
  }
  else if(do_collect_check){
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //-c
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(argc-optind<3)
      throw std::invalid_argument("3 non-option arguments are required: which response to use, digest file, and location to write edge list file\n");
    deker::io::Control_Digest digest;
    unsigned input_response(deker::io::read_arg_to_val <unsigned> (argv[optind]));
    digest.read_single_record(argv[optind+1],input_response);
    std::experimental::filesystem::path digest_dir_path(digest.get_digest_dir());
    //
    deker::io::Eigen_Matrix_IO data_feature_names(digest.get_control().data_file,true,false,true);
    unsigned which_response = digest.get_digest_record(0).which_response;
    //
    std::string lgbm_file_loc = (digest_dir_path / std::string(digest.get_digest_record(0).lgbm_file)).string();
    std::string lopt_file_loc = (digest_dir_path / std::string(digest.get_digest_record(0).opt_file)).string();
    std::string run_file_loc = (digest_dir_path / std::string(digest.get_digest_record(0).run_file)).string();
    std::string output_file_loc = (digest_dir_path / std::string(digest.get_digest_record(0).output_file)).string();
    //
    bool run_error = deker::io::file_exists_test(run_file_loc);
    bool run_finished = deker::io::file_exists_test(output_file_loc);
    if(run_error||run_finished){
      //Setting up output to write
      std::vector<std::string> extra_id_header;
      extra_id_header.push_back("response_digest_key");
      extra_id_header.push_back("digest_filename");
      extra_id_header.push_back("edgelist_filename");
      std::vector<std::string> extra_id;
      extra_id.push_back(argv[optind]);
      extra_id.push_back(argv[optind+1]);
      extra_id.push_back(argv[optind+2]);
      if(run_error){
        deker::fit::Output_deker_IO to_write;
        //write error to csv
        to_write.response_name = data_feature_names.feature_names.at(which_response);
        to_write.feature_names.push_back("ERROR");
        to_write.v.resize(1);
        to_write.v(0) = 0;
        to_write.lgbm_imp.resize(1);
        to_write.lgbm_imp(0) = 0;
        to_write.write_to_csv(argv[optind+2],extra_id_header,extra_id);
        std::cout<<"1\n"; //feedback - response is done
      }else if(run_finished){
        //read in all the material for edgelist output
        //lambda_opt
        using Lambda_Opt_t = deker::fit::deker_Lambda_Optimizer<deker::fit::Solve_deker,&deker::fit::Solve_deker::inner_full_opt>;
        Lambda_Opt_t lambda_opt;
        deker::io::read_from_binary<Lambda_Opt_t>(lopt_file_loc,lambda_opt);
        //lgbm_output  
        deker::LGBM_out lgbm_output;
        deker::io::read_from_binary<deker::LGBM_out>(lgbm_file_loc,lgbm_output);
        //
        deker::fit::Output_deker_IO to_write(lambda_opt.get_best_sol(),
                                             data_feature_names.feature_names, //labels input data
                                             digest.get_control().data_file, //write location
                                             which_response, //input
                                             lambda_opt.get_null_BIC(),
                                             lgbm_output.feature_imp_gain, //from lgbm
                                             lgbm_output.fit_metric); //from lgbm
        to_write.write_to_csv(argv[optind+2],extra_id_header,extra_id);
        std::cout<<"1\n"; //feedback - response is done
      }
    }else{
      std::cout<<"0\n"; //feedback - response is not done
    }
    /////////////////////////////////////
  }
  else if(do_digest_check){
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //-d
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      if(argc-optind<1)
          throw std::invalid_argument("2 non-option arguments is required: digest file location and control file write location\n");
      deker::io::Control_Digest digest;
      digest.read_single_record(argv[optind],0);
      digest.get_control().write_control_file(argv[optind+1]);
      /////////////////////////////////////
  }
  else if(make_edge_list){
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //-e
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      if(argc-optind<3)
        throw std::invalid_argument("3 non-option arguments are required: which response to use, digest file, and location to write edge list file\n");
      //
      deker::io::Control_Digest digest;
      unsigned input_response(deker::io::read_arg_to_val <unsigned> (argv[optind]));
      digest.read_single_record(argv[optind+1],input_response);
      std::experimental::filesystem::path digest_dir_path(digest.get_digest_dir());
      //
      deker::io::Eigen_Matrix_IO data_feature_names(digest.get_control().data_file,true,false,true);
      unsigned which_response = digest.get_digest_record(0).which_response;
      //
      std::string lgbm_file_loc = (digest_dir_path / std::string(digest.get_digest_record(0).lgbm_file)).string();
      std::string lopt_file_loc = (digest_dir_path / std::string(digest.get_digest_record(0).opt_file)).string();
      //
      std::vector<std::string> extra_id_header;
      extra_id_header.push_back("response_digest_key");
      extra_id_header.push_back("model_key");
      extra_id_header.push_back("digest_filename");
      extra_id_header.push_back("edgelist_filename");
      std::vector<std::string> extra_id;
      extra_id.push_back(argv[optind]);
      extra_id.push_back("0");
      extra_id.push_back(argv[optind+1]);
      extra_id.push_back(argv[optind+2]);
      //read in all the material for edgelist output
      //lambda_opt
      using Lambda_Opt_t = deker::fit::deker_Lambda_Optimizer<deker::fit::Solve_deker,&deker::fit::Solve_deker::inner_full_opt>;
      Lambda_Opt_t lambda_opt;
      deker::io::read_from_binary<Lambda_Opt_t>(lopt_file_loc,lambda_opt);
      //lgbm_output  
      deker::LGBM_out lgbm_output;
      deker::io::read_from_binary<deker::LGBM_out>(lgbm_file_loc,lgbm_output);
      //
      if(all_out){
        for(unsigned i = 0; i<lambda_opt.get_sample_output().size(); i++){
          extra_id.at(1) = std::to_string(i);
          deker::fit::Output_deker_IO to_write(lambda_opt.get_sample_output().at(i),
                                               data_feature_names.feature_names, //labels input data
                                               digest.get_control().data_file, //write location
                                               which_response, //input
                                               lambda_opt.get_null_BIC(),
                                               lgbm_output.feature_imp_gain, //from lgbm
                                               lgbm_output.fit_metric); //from lgbm
          to_write.write_to_csv(argv[optind+2],extra_id_header,extra_id);
        }
      }else{
        extra_id.at(1) = std::to_string(lambda_opt.get_which_best());
        deker::fit::Output_deker_IO to_write(lambda_opt.get_best_sol(),
                                             data_feature_names.feature_names, //labels input data
                                             digest.get_control().data_file, //write location
                                             which_response, //input
                                             lambda_opt.get_null_BIC(),
                                             lgbm_output.feature_imp_gain, //from lgbm
                                             lgbm_output.fit_metric); //from lgbm
        to_write.write_to_csv(argv[optind+2],extra_id_header,extra_id);
      }
  }
  else if(do_fit){
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // -f
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    unsigned rf_n_trees = 1000;
    //
    time_t start_time;
    start_time = time(NULL);
    /////////////////////////////////////
    //values needed to run fit
    if(argc-optind<2)
      throw std::invalid_argument("2 non-option arguments are required: which response to use and digest file");
    deker::io::Control_Digest digest;
    unsigned input_response(deker::io::read_arg_to_val <unsigned> (argv[optind]));
    std::cerr<<argv[optind+1]<<" "<<input_response<<"\n";
    digest.read_single_record(argv[optind+1],input_response);
    //
    deker::io::Eigen_Matrix_IO input_data;
    deker::io::read_from_binary<deker::io::Eigen_Matrix_IO>(digest.get_control().data_file,input_data);
    //
    if(digest.get_control().exclude_sample.size()>0)
      input_data.drop_rows(digest.get_control().exclude_sample);
    //
    unsigned which_response = digest.get_digest_record(0).which_response;
    std::cerr<<"which_response "<<which_response<<"\n";
    std::experimental::filesystem::path digest_dir_path(digest.get_digest_dir());
    //
    std::string sol_file_loc = (digest_dir_path / std::string(digest.get_digest_record(0).solve_file)).string();
    std::string lgbm_file_loc = (digest_dir_path / std::string(digest.get_digest_record(0).lgbm_file)).string();
    std::string lopt_file_loc = (digest_dir_path / std::string(digest.get_digest_record(0).opt_file)).string();
    std::string run_file_loc = (digest_dir_path / std::string(digest.get_digest_record(0).run_file)).string();
    std::string output_file_loc = (digest_dir_path / std::string(digest.get_digest_record(0).output_file)).string();
    //
    //Check if run file exists - means a run was initiated and not completed
    if(deker::io::file_exists_test(run_file_loc)){
      throw std::invalid_argument("run-lock file found for specified response: previous iteration terminated before completion\n");
    }
    //Write run file - run initiated
    std::ofstream outfile(run_file_loc);
    outfile.close();
    ////////////////////////////////////////////////
    //build Solve_deker
    deker::fit::Solve_deker sol_deker(input_data.matrix,digest.get_digest_record(0).which_response);
    if(deker::io::file_exists_test(sol_file_loc)){
      std::cerr<<"read solve check\n";
      deker::io::read_from_binary<deker::fit::Solve_deker>(sol_file_loc,sol_deker);
      std::cerr<<"read solve clear\n";
    }else{
      std::vector<unsigned> use_predictors = deker::io::build_use_vectors(digest.get_control().ex_pred,digest.get_control().in_pred,(unsigned) input_data.cols);
      //if use_predictors is empty, make a vector of indicies for data columns - 0:(input_data.matrix.cols()-1)
      if(use_predictors.size()==0){
        for(unsigned i = 0; i<input_data.cols; i++)
          use_predictors.push_back(i);
      } 
      sol_deker.build(use_predictors,
                      digest.get_control().opt_param);
      deker::io::write_to_binary<deker::fit::Solve_deker>(sol_file_loc,sol_deker);
    }
    ////////////////////////////////////////////////
    //run LGBM
    bool continue_fit_from_lgbm = true;
    //deker::LGBM_out lgbm_output;
    //deker::io::read_from_binary<deker::LGBM_out>(lgbm_file_loc,lgbm_output);
    if(!deker::io::file_exists_test(lgbm_file_loc)){
      char *parameters = "num_threads=1 objective=regression min_data_in_leaf=5 force_col_wise=true min_gain_to_split=.01 boosting=rf bagging_fraction=.666 bagging_freq=1 verbosity=-1";
      deker::LGBM_wrapper booster(sol_deker.get_input_data().predictors,sol_deker.get_input_data().response,parameters);
      booster.build_booster();
      booster.run(rf_n_trees);
      //
      deker::LGBM_out lgbm_output = booster.get_output();
      lgbm_output.feature_imp_gain = lgbm_output.feature_imp_gain.array().sqrt()/((double)lgbm_output.iter);
      lgbm_output.fit_metric = sqrt(lgbm_output.fit_metric);
      //
      deker::io::write_to_binary<deker::LGBM_out>(lgbm_file_loc,lgbm_output);
      //
      // double lgbm_wt_cutoff = .065;
      // Eigen::VectorXd::Index which_max;
      // if(lgbm_output.feature_imp_split.maxCoeff(&which_max)<pow(lgbm_wt_cutoff*((double)lgbm_output.iter),2.0)){
      //  //break if no features above split cutoff
      //  //implement if/when a good cutoff can be defined
      // }
        time_t current_time;
        current_time = time(NULL);
        if(digest.get_control().max_job_time>0&&
           (current_time-start_time)/3600>digest.get_control().max_job_time){
          continue_fit_from_lgbm = false;
        }
    }
    if(continue_fit_from_lgbm){
      ////////////////////////////////////////////////
      //build Lambda_Optimizer
      using Lambda_Opt_t = deker::fit::deker_Lambda_Optimizer<deker::fit::Solve_deker,&deker::fit::Solve_deker::inner_full_opt>;
      Lambda_Opt_t lambda_opt;
      if(deker::io::file_exists_test(lopt_file_loc)){
        std::cerr<<"read lambda opt check\n";
        deker::io::read_from_binary<Lambda_Opt_t>(lopt_file_loc,lambda_opt);
        std::cerr<<"read lambda opt clear\n";
      }else{
        lambda_opt.build(digest.get_control().opt_param);
        lambda_opt.set_null_BIC(sol_deker.get_null_BIC());
      }
      ////////////////////////////////////////////////
      //run some Lambda_Optimizer iterations
      bool rangefind_complete, fill_complete, bo_complete;
      std::tie(rangefind_complete, fill_complete, bo_complete) = lambda_opt.get_status();
      bool still_opting = !(rangefind_complete&&fill_complete&&bo_complete);
      int lambda_iters = 0;
      std::cerr<<digest.get_control().lambda_iter_per_job<<"\n";
      while(still_opting){
        ///
        unsigned rf_it,fi_it,bo_it;
        std::tie(rf_it,fi_it,bo_it) = lambda_opt.get_iters();
        std::cerr<<"lambda total iters: "<<rf_it<<" "<<fi_it<<" "<<bo_it<<"\n";
        ///
        time_t init_time;
        init_time = time(NULL);
        still_opting = lambda_opt.opt_iteration(sol_deker);
        time_t finish_time;
        finish_time = time(NULL);
        lambda_iters++;
        //break if specified by max_time
        double predicted_elapsed_time = finish_time-start_time+finish_time-init_time;
        predicted_elapsed_time /= 3600;
        if(digest.get_control().max_job_time>0&&
           predicted_elapsed_time>digest.get_control().max_job_time){
          break;
        }
        //break if specified by iter_per_job
        if(digest.get_control().lambda_iter_per_job>0&&
           lambda_iters==digest.get_control().lambda_iter_per_job){
          break;
        }
      }
      ///
      unsigned rf_it,fi_it,bo_it;
      std::tie(rf_it,fi_it,bo_it) = lambda_opt.get_iters();
      std::cerr<<"lambda total iters: "<<rf_it<<" "<<fi_it<<" "<<bo_it<<"\n";
      ///
      deker::io::write_to_binary<Lambda_Opt_t>(lopt_file_loc,lambda_opt);
      ////////////////////////////////////////////////
      //write output
      if(!still_opting){
        //write output signal file
        std::ofstream outfile(output_file_loc);
        outfile.close();
      }
    }
    remove(run_file_loc.c_str());
    /////////////////////////////////////
  }
  else if(do_manual_fit){
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //-m
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //values needed to run fit
    if(argc-optind<4)
      throw std::invalid_argument("4 non-option arguments are required: output location, data file, which response to use, and lambda");
    deker::io::Eigen_Matrix_IO input_data;
    deker::io::read_from_binary<deker::io::Eigen_Matrix_IO>(argv[optind+1],input_data);
    unsigned which_response(deker::io::read_arg_to_val <unsigned> (argv[optind+2]));
    double lambda_regularization(deker::io::read_arg_to_val <double> (argv[optind+3]));
    std::cerr<<"manual fitting: "<<argv[optind+1]<<" feature "<<which_response<<" lambda "<<lambda_regularization<<"\n";
    //
    std::vector<unsigned> use_predictors;
    for(unsigned i = 0; i<input_data.cols; i++)
      use_predictors.push_back(i);
    //
    deker::fit::Opt_Param default_param;
    deker::fit::Solve_deker sol_deker(input_data.matrix,which_response);
    sol_deker.build(use_predictors,
                    default_param);
    //
    deker::fit::Output_deker sol;
    sol_deker.inner_full_opt(lambda_regularization,sol);
    //////////////////////////////////////////////////////////////
    double transformed_fit = 0;
    Eigen::VectorXd transformed_feature_importance = Eigen::VectorXd::Zero(sol_deker.get_input_data().predictors.cols());
    //write output to file
    deker::fit::Output_deker_IO to_write(sol,
                                         input_data.feature_names,
                                         argv[optind+1],
                                             which_response,
                                             sol_deker.get_null_BIC(),
                                             transformed_feature_importance,
                                             transformed_fit);
    //
    std::vector<std::string> extra_id_header;
    extra_id_header.push_back("edgelist_filename");
    extra_id_header.push_back("source_filename");
    std::vector<std::string> extra_id;
    extra_id.push_back(argv[optind]);
    extra_id.push_back(argv[optind+1]);
    to_write.write_to_csv(argv[optind],extra_id_header,extra_id);
    //
    deker::fit::Predictions_deker predict(sol_deker.get_input_data(),
                                          sol_deker.get_response_labels(),
                                          sol.sigma_kernel_width,
                                          sol.v,
                                          sol.feature_index);
    predict.write_to_csv(argv[optind]);
    /////////////////////////////////////
  }
  else if(do_digest_count){
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //-n
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(argc-optind<1)
      throw std::invalid_argument("1 non-option arguments is required: digest file");
    /////
    deker::io::Control_Digest digest;
    digest.read_digest_size(argv[optind]);
    std::cout<<digest.get_full_record_size()<<"\n";
  }
  else if(predict_out){
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //-p
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(argc-optind<2)
      throw std::invalid_argument("3 non-option arguments are required: which response to use, digest file, and where to write");
    /////
    deker::io::Control_Digest digest;
    unsigned input_response(deker::io::read_arg_to_val <unsigned> (argv[optind]));
    std::cerr<<argv[optind+1]<<" "<<input_response<<"\n";
    digest.read_single_record(argv[optind+1],input_response);
    deker::io::Eigen_Matrix_IO input_data;
    deker::io::read_from_binary<deker::io::Eigen_Matrix_IO>(digest.get_control().data_file,input_data);
    //
    if(do_holdout_predict&digest.get_control().exclude_sample.size()==0)
      throw std::invalid_argument("prediction of held-out samples requested, but no held-out samples found");
    //
    Eigen::MatrixXd test_predictors;
    Eigen::VectorXd test_response;
    if(do_holdout_predict){
      test_predictors.resize(digest.get_control().exclude_sample.size(),input_data.cols);
      test_response.resize(digest.get_control().exclude_sample.size());
      for(unsigned i = 0; i<digest.get_control().exclude_sample.size(); i++){
        test_predictors.row(i) = input_data.matrix.row(digest.get_control().exclude_sample.at(i));
        test_response(i) = input_data.matrix(digest.get_control().exclude_sample.at(i),digest.get_digest_record(0).which_response);
      }
    }
    //
    if(digest.get_control().exclude_sample.size()>0)
      input_data.drop_rows(digest.get_control().exclude_sample);
    //
    unsigned which_response = digest.get_digest_record(0).which_response;
    std::cerr<<"which_response "<<which_response<<"\n";
    std::experimental::filesystem::path digest_dir_path(digest.get_digest_dir());
    std::string sol_file_loc = (digest_dir_path / std::string(digest.get_digest_record(0).solve_file)).string();
    std::string lopt_file_loc = (digest_dir_path / std::string(digest.get_digest_record(0).opt_file)).string();
    /////
    if(!deker::io::file_exists_test(sol_file_loc)||!deker::io::file_exists_test(lopt_file_loc)){
      throw std::invalid_argument("sol and lopt files not found for specified response: cannot output prediction\n");
    }
    //
    std::vector<std::string> extra_id_header;
    extra_id_header.push_back("response_digest_key");
    extra_id_header.push_back("model_key");
    extra_id_header.push_back("digest_filename");
    extra_id_header.push_back("edgelist_filename");
    std::vector<std::string> extra_id;
    extra_id.push_back(argv[optind]);
    extra_id.push_back("0");
    extra_id.push_back(argv[optind+1]);
    extra_id.push_back(argv[optind+2]);
    /////
    deker::fit::Solve_deker sol_deker(input_data.matrix,digest.get_digest_record(0).which_response);
    deker::io::read_from_binary<deker::fit::Solve_deker>(sol_file_loc,sol_deker);
    //
    if(do_holdout_predict){ //align test_predictors with sol_deker.get_input_data()
      Eigen::MatrixXd tmp_test_predictors = test_predictors;
      test_predictors.resize(test_predictors.rows(),sol_deker.get_input_data().key_to_original_data.size());
      for(unsigned i = 0; i<sol_deker.get_input_data().key_to_original_data.size();i++){
        test_predictors.col(i) = tmp_test_predictors.col(sol_deker.get_input_data().key_to_original_data.at(i));
      }
    }
    //
    using Lambda_Opt_t = deker::fit::deker_Lambda_Optimizer<deker::fit::Solve_deker,&deker::fit::Solve_deker::inner_full_opt>;
    Lambda_Opt_t lambda_opt;
    deker::io::read_from_binary<Lambda_Opt_t>(lopt_file_loc,lambda_opt);
    //
    std::cerr<<"check1\n";
    if(all_out){
      for(unsigned i = 0; i<lambda_opt.get_sample_output().size();i++){
        extra_id.at(1) = std::to_string(i);
        if(lambda_opt.get_sample_output().at(i).return_status==1||lambda_opt.get_sample_output().at(i).return_status==3){
          if(do_holdout_predict){
            deker::fit::Predictions_deker_holdout predict(sol_deker.get_input_data(),
                                                          sol_deker.get_response_labels(),
                                                          lambda_opt.get_sample_output().at(i).sigma_kernel_width,
                                                          lambda_opt.get_sample_output().at(i).v,
                                                          lambda_opt.get_sample_output().at(i).feature_index,
                                                          test_predictors,
                                                          test_response);
            predict.write_to_csv(argv[optind+2],extra_id_header,extra_id);
          }else{
            deker::fit::Predictions_deker predict(sol_deker.get_input_data(),
                                                  sol_deker.get_response_labels(),
                                                  lambda_opt.get_sample_output().at(i).sigma_kernel_width,
                                                  lambda_opt.get_sample_output().at(i).v,
                                                  lambda_opt.get_sample_output().at(i).feature_index);
              std::cerr<<"bonehead\n";
            predict.write_to_csv(argv[optind+2],extra_id_header,extra_id);
          }
        }
      }
    }else{
      extra_id.at(1) = std::to_string(lambda_opt.get_which_best());
      if(do_holdout_predict){
        deker::fit::Predictions_deker_holdout predict(sol_deker.get_input_data(),
                                                      sol_deker.get_response_labels(),
                                                      lambda_opt.get_best_sol().sigma_kernel_width,
                                                      lambda_opt.get_best_sol().v,
                                                      lambda_opt.get_best_sol().feature_index,
                                                      test_predictors,
                                                      test_response);
        predict.write_to_csv(argv[optind+2],extra_id_header,extra_id);
      }else{
        deker::fit::Predictions_deker predict(sol_deker.get_input_data(),
                                              sol_deker.get_response_labels(),
                                              lambda_opt.get_best_sol().sigma_kernel_width,
                                              lambda_opt.get_best_sol().v,
                                              lambda_opt.get_best_sol().feature_index);
        predict.write_to_csv(argv[optind+2],extra_id_header,extra_id);
      }
    }
  }
  else if(start_control){
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //-s
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(argc-optind<2)
      throw std::invalid_argument("2 non-option argument are required: control file location & location to write digest\n");
    if(!deker::io::file_exists_test(argv[optind]))
      throw std::invalid_argument("Control file not found\n");
    if(deker::io::file_exists_test(argv[optind+1]))
      throw std::invalid_argument("File found at digest file location\n");
    deker::io::Control_Base control;
    control.read_from_txt_file(argv[optind]);
    deker::io::Eigen_Matrix_IO sizes_only_data(control.data_file,true,true);
    std::vector<unsigned> use_response = deker::io::build_use_vectors(control.ex_res,control.in_res,(unsigned) sizes_only_data.cols);
    if(use_response.size()==0)
      throw std::invalid_argument("No responses to run\n");
    //
    std::experimental::filesystem::path p(argv[optind+1]);
    std::string digest_dir = p.parent_path().string();
    std::experimental::filesystem::create_directory(digest_dir);
    deker::io::Control_Digest digest(digest_dir,control);
    for(unsigned i = 0; i<use_response.size();i++){
      digest.add_record(use_response.at(i));
    }
    deker::io::write_to_binary<deker::io::Control_Digest>(argv[optind+1],digest);
    //
    std::cout<<use_response.size()<<"\n";
    /////////////////////////////////////
  }
  else{
    /////////////////////////////////////
    std::cout<<"usage/help TODO\n";
    /////////////////////////////////////
  }
  /////////////////////////////////////
  /////////////////////////////////////
  return 0;
}