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
#ifndef DEKER_FIT
#define DEKER_FIT
///////////////////////////////////////////
#include <deker/fit/math.h>
#include <deker/fit/output.h>
#include <deker/fit/predictions.h>
#include <deker/fit/limbo_params.h>
#include <fstream>
///////////////////////////////////////////
namespace deker{
  namespace fit{
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    class Solve_deker{
    protected:
      Split_Data input_data;
      Response_Labels response_labels;
      Opt_Param control;
      double null_BIC;
      //
    private:
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      struct v_opt_prob{
        Eigen::MatrixXd z_cache;
        Eigen::VectorXd z_cache_col_multiple;
        const Solve_deker* parent;
        const double& lambda_regularization;
        ///////////
        v_opt_prob(const double& lambda_regularization,
                   const double& sigma_kernel_width,
                   const Eigen::Ref<const Eigen::VectorXd> v,
                   const std::vector<unsigned>& feature_index,
                   const Solve_deker* parent):
          lambda_regularization(lambda_regularization),parent(parent)
        {
          //std::vector<Eigen::MatrixXd> threshold_diff_mat(math::calc_threshold_diff_mat(parent->get_input_data(),parent->get_response_labels(),sigma_kernel_width,v,feature_index,true));
          std::tie(z_cache,z_cache_col_multiple) = math::calc_z_cache(parent->get_input_data(),parent->get_response_labels(),sigma_kernel_width,v,feature_index,true,true);
          //std::tie(z_cache,z_cache_col_multiple) = math::calc_z_cache(parent->get_input_data(),parent->get_response_labels(),sigma_kernel_width,v,feature_index,true); //rebuild z, excluding train from test
        }
        ///////////
        limbo::opt::eval_t operator() (const Eigen::Ref<const Eigen::VectorXd> opt_vec, bool eval_grad = false) const
        {
          //optimizing v
          double fx;
          Eigen::VectorXd gradient;
          //updates fx, gradient
          math::calculate_v_gradient(parent->control,opt_vec,z_cache,z_cache_col_multiple,lambda_regularization,fx,gradient);
          //because we're maximizing
          fx *= -1;
          gradient *= -1;
          return {fx, gradient};
        }
      };
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      struct sigma_opt_prob{
        const double sigma_min;
        const double sigma_max;
        const Eigen::Ref<const Eigen::VectorXd> v_current;
        const std::vector<unsigned>& feature_index_current;
        const Solve_deker* parent;
        ///////////
        sigma_opt_prob(const double& sigma_min,
                       const double& sigma_max,
                       const Eigen::Ref<const Eigen::VectorXd> v_current,
                       const std::vector<unsigned>& feature_index_current,
                       const Solve_deker* parent):
          sigma_min(sigma_min), sigma_max(sigma_max),v_current(v_current),feature_index_current(feature_index_current),parent(parent){}
        ///////////
        double convert_to_bounded(const double& unbounded_sigma) const 
        {
          return (unbounded_sigma-sigma_min)/(sigma_max-sigma_min);
        }
        ///////////
        double convert_from_bounded(const double& bounded_sigma) const
        {
          return bounded_sigma*(sigma_max-sigma_min)+sigma_min;
        }
        ///////////
        limbo::opt::eval_t operator() (const Eigen::Ref<const Eigen::VectorXd> opt_vec, bool eval_grad = false) const
        {
          double use_sigma = convert_from_bounded(opt_vec(0));
          double tmp_df, tmp_logli;
          std::tie(tmp_df,tmp_logli) = parent->calc_model_fit(use_sigma,v_current,feature_index_current);
          double fx = -log((double) parent->input_data.response.size())*tmp_df+2*tmp_logli;
          //std::cerr<<use_sigma<<" "<<tmp_df<<" "<<tmp_logli<<" "<<fx<<"\n";
          return limbo::opt::no_grad(fx);
        }
      };
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
    public:
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      Solve_deker(Eigen::MatrixXd& data_matrix,
                  const unsigned& which_response):
      input_data(data_matrix,which_response),
      response_labels(input_data.response)
      {}
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      Solve_deker(Eigen::MatrixXd& data_matrix,
                  const std::vector<unsigned>& predictors_use_index,
                  const unsigned& which_response,
                  const Opt_Param& input_control): 
      input_data(data_matrix,which_response,predictors_use_index),
      response_labels(input_data.response)
      {
        build(predictors_use_index,input_control);
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      void build(const std::vector<unsigned>& predictors_use_index,
                 const Opt_Param& input_control)
      {
        input_data.build(predictors_use_index);
        response_labels.build();
        null_BIC = math::calc_null_BIC(response_labels);
        control = input_control;
        //do some stuff to drop obviously extraneous features
        if(input_data.predictors_index.size()>2){
          std::vector<unsigned> initial_feature_index = input_data.predictors_index; //feature_index
          Eigen::VectorXd initial_v = Eigen::VectorXd::Constant(initial_feature_index.size(),1); //initial v
          double initial_sigma = opt_sigma(initial_v,initial_feature_index);
          double initial_lambda = 0;
          //builds z_cache etc.
          v_opt_prob v_opt_tmp(initial_lambda,initial_sigma,initial_v,initial_feature_index,this);
          //get column sums & censored (sum for classification problems where margin fails)
          Eigen::VectorXd z_col_sum = (v_opt_tmp.z_cache.array().transpose().colwise()*v_opt_tmp.z_cache_col_multiple.array()).colwise().sum()/(v_opt_tmp.z_cache_col_multiple.sum());
          Eigen::VectorXd cens_z_col_sum = Eigen::VectorXd::Zero(v_opt_tmp.z_cache.rows());
          double cens_cols = 0;
          for(unsigned i=0;i<v_opt_tmp.z_cache.cols();i++){
            if(v_opt_tmp.z_cache.col(i).sum()<=1.0){
              cens_z_col_sum+=v_opt_tmp.z_cache.col(i)*v_opt_tmp.z_cache_col_multiple(i);
              cens_cols+=v_opt_tmp.z_cache_col_multiple(i);
            }
          }
          if(cens_cols>0){
            cens_z_col_sum.array()/=cens_cols;
          }
          double elbow_sum = math::find_elbow(z_col_sum);
          double elbow_cens_sum = math::find_elbow(cens_z_col_sum);
          ///
          //debug only
          //std::cerr<<elbow_sum<<" "<<elbow_cens_sum<<"\n";
          //std::vector<unsigned> keep_features;
          //drop features that aren't contributing to total & censored - these will get zero'd out in the first run anyway
          for(int i=initial_feature_index.size()-1;i>=0;i--){
            if(z_col_sum(i)<elbow_sum&&cens_z_col_sum(i)<elbow_cens_sum){
              if(i < initial_feature_index.size()-1){
                initial_feature_index.at(i)=initial_feature_index.back();
              }
              initial_feature_index.pop_back();
            }
            //else{
            //debug only
            //  keep_features.push_back(i);
            //}
          }
          //out for debugging
          // for(unsigned i = 0;i<initial_feature_index.size();i++){
          //   std::cerr<<
          //     //initial_feature_index.at(i)<<" "<<
          //     input_data.predictors_index.at(keep_features.at(i))<<" "<<
          //       z_col_sum(keep_features.at(i))<<" "<<
          //         cens_z_col_sum(keep_features.at(i))<<"\n";
          // }
          input_data.subset_predictors(initial_feature_index);
        }
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      std::tuple<double,double> calc_model_fit(const double& sigma_kernel_width,
                                                             const Eigen::Ref<const Eigen::VectorXd> v,
                                                             const std::vector<unsigned>& feature_index) const
      {
        Predictions_deker predict(input_data,response_labels,sigma_kernel_width,v,feature_index);
        return predict.get_fit();
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      double opt_sigma(const Eigen::Ref<const Eigen::VectorXd> v,
                       const std::vector<unsigned>& feature_index) const
      {
        Eigen::MatrixXd temp_kernel = math::calc_kernel_matrix(input_data,0,v,feature_index);
        double max_el = -temp_kernel.minCoeff();
        if(max_el<=0){
          return 0;
        }
        //smallest_diff will be 10 or less, but greater than 0
        double smallest_diff = (temp_kernel.array() == -max_el).select(10.0,temp_kernel.array()+max_el).minCoeff();
        double min_sigma = .5*(log(smallest_diff/2.0)-log(1e5)); //sigma where the smallest difference is very large (only nearest neighbor)
        double max_sigma = .5*(log((max_el)/2.0)-log(1e-5)); //sigma where the largest difference is very small (global average)
        //find best sigma given v
        sigma_opt_prob sigma_opt_tmp(min_sigma,max_sigma,v,feature_index,this);
        Eigen::VectorXd sigma_init = Eigen::VectorXd::Constant(1,sigma_opt_tmp.convert_to_bounded((max_sigma+min_sigma)/2.0));
        limbo::opt::NLOptNoGrad<PF_Params,nlopt::LN_BOBYQA> limbo_nograd_solver;
        Eigen::VectorXd sigma_local_solution = limbo_nograd_solver(sigma_opt_tmp,sigma_init,true);
        return sigma_opt_tmp.convert_from_bounded(sigma_local_solution(0));
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      Eigen::VectorXd opt_v(const double& lambda_regularization,
                            const double& sigma_kernel_width,
                            const Eigen::Ref<const Eigen::VectorXd> v, 
                            const std::vector<unsigned>& feature_index) const
      {
        //find best v given sigma
        v_opt_prob v_opt_tmp(lambda_regularization,sigma_kernel_width,v,feature_index,this);
        limbo::opt::NLOptGrad<FS_Params> limbo_lbfgs_solver;
        Eigen::VectorXd v_new = limbo_lbfgs_solver(v_opt_tmp,v,false);
        return v_new.array().abs();
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      void steff_iter_v(const double& lambda_regularization,
                        Output_deker& sol) const
      {
        sol.return_status = 0; //nothing to report
        Eigen::MatrixXd v_i(sol.v.size()+1,4);
        v_i(0,0) = sol.sigma_kernel_width;
        v_i.block(1,0,sol.v.size(),1) = sol.v;
        //std::cerr<<" (v) "<<v_i.col(0).transpose()<<"\n";
        for(unsigned i = 1; i<v_i.cols(); i++){
          //find best sigma given v
          //std::cerr<<"***********************************\n";
          v_i(0,i) = opt_sigma(v_i.block(1,i-1,sol.v.size(),1),sol.feature_index);
          //find best v given sigma
          v_i.block(1,i,sol.v.size(),1) = opt_v(lambda_regularization,v_i(0,i),v_i.block(1,i-1,sol.v.size(),1),sol.feature_index);
          //std::cerr<<" (v) "<<v_i.col(i).transpose()<<"\n";
          //Early exit if fixed point is found in these iterations; no need to do steffensen, which may break
          double v_diff_test = math::diff_test(v_i.block(1,i-1,sol.v.size(),1),v_i.block(1,i,sol.v.size(),1));
          if(v_diff_test <= control.w_convergence_threshold && std::abs(v_i(0,i-1)-v_i(0,i))<1e-4){ 
            sol.sigma_kernel_width = v_i(0,i);
            sol.v = v_i.block(1,i,sol.v.size(),1);
            sol.return_status = 1; //w is a fixed point
            math::update_v_with_drop(control.feature_drop_threshold,
                                     sol.feature_index, //updates
                                     v_i.block(1,i-1,sol.v.size(),1),  
                                     sol.v); //updates
            break;
          }
        }
        if(sol.return_status==0){
          //use steffensen acceleration on the 4 observations to generate the next w in the fixed point iteration
          v_i.col(3) = math::steffensen_fixed_point_acceleration(v_i);
          sol.sigma_kernel_width = v_i(0,3);
          sol.v = v_i.block(1,3,sol.v.size(),1);
          double v_diff_test = math::diff_test(v_i.block(1,0,sol.v.size(),1),v_i.block(1,3,sol.v.size(),1));
          if(v_diff_test <= control.w_convergence_threshold && std::abs(v_i(0,0)-v_i(0,3))<1e-4){ 
            sol.return_status = 1;
          }
          math::update_v_with_drop(control.feature_drop_threshold,
                                   sol.feature_index, //updates
                                   v_i.block(1,0,sol.v.size(),1),  
                                   sol.v); //updates
        }
        if(sol.feature_index.size()==0){
          sol.return_status = 2;
        }  
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      void inner_full_opt(const double& lambda_regularization,
                          Output_deker& sol) const
      {
        std::cerr<<"***********\n";
        //set lambda in container
        sol.lambda_regularization = lambda_regularization;
        //initial values for v & sigma
        sol.feature_index = input_data.predictors_index;
        sol.key_to_original_data = input_data.key_to_original_data;
        sol.v = Eigen::VectorXd::Constant(sol.feature_index.size(),1);
        sol.sigma_kernel_width = 5;
        sol.return_status = 0;
        sol.df = 0;
        sol.logli = 0;
        sol.BIC = 0;
        ///
        int iter = 0;
        for(; iter<control.w_maxiter; iter++){
          steff_iter_v(lambda_regularization,sol);
          std::cerr<<"  sigma "<<sol.sigma_kernel_width<<" v "<<sol.v.transpose()<<"\n";
          if(sol.return_status != 0) break;
        }
        if(iter >= control.w_maxiter){
          sol.return_status = 3;
        }
        if(sol.return_status != 2){ //feature_index.size() > 0
          
          Predictions_deker predict(input_data,response_labels,sol.sigma_kernel_width,sol.v,sol.feature_index);
          //double df,lf_error;
          std::tie(sol.df,sol.logli) = predict.get_fit();
          sol.BIC = log((double)input_data.response.size())*sol.df-2*sol.logli;
          
          std::cerr<<"***********\n";
          std::cerr<<"full_opt df "<<sol.df<<" lf_error "<<sol.logli<<" BIC "<<sol.BIC<<"\n";
        }
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      const Split_Data& get_input_data() const {return input_data;}
      const Response_Labels& get_response_labels() const {return response_labels;}
      const Opt_Param& get_opt_param() const {return control;}
      double get_null_BIC() const {return null_BIC;}
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      void serialize(std::ifstream& infile)
      {
        input_data.serialize(infile);
        response_labels.serialize(infile);
        control.serialize(infile);
        infile.read((char*) (&null_BIC),sizeof(null_BIC));
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      void serialize(std::ofstream& outfile) const
      {
        input_data.serialize(outfile);
        response_labels.serialize(outfile);
        control.serialize(outfile);
        outfile.write((char*) (&null_BIC),sizeof(null_BIC));
      }
    };
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
  }
}
#endif