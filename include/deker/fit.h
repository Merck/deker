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
#include <deker/fit/platt_scaling.h>
#include <deker/fit/limbo_params.h>
///////////////////////////////////////////
namespace deker{
  namespace fit{
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    class Solve_deker{
    private:
      const Split_Data input_data;
      const Response_Labels response_labels;
      const Opt_Param control;
      //
      limbo::opt::NLOptNoGrad<PF_Params,nlopt::LN_BOBYQA> limbo_nograd_solver;
      limbo::opt::NLOptNoGrad<PF_Params,nlopt::GN_DIRECT_L> limbo_global_nograd_solver;
      limbo::opt::NLOptGrad<FS_Params> limbo_lbfgs_solver;
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      struct v_opt_prob{
        Eigen::MatrixXd z_cache;
        Eigen::VectorXd w_transform;
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
          std::tie(z_cache,w_transform,z_cache_col_multiple) = math::calc_z_cache(parent->get_input_data(),parent->get_response_labels(),sigma_kernel_width,v,feature_index,true); //rebuild z, excluding train from test
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
          double tmp_df, tmp_logli,tmp_A,tmp_B;
          std::tie(tmp_df,tmp_logli,tmp_A,tmp_B) = parent->calc_model_fit(use_sigma,v_current,feature_index_current);
          double fx = -log((double) parent->input_data.response.size())*tmp_df+2*tmp_logli;
          return limbo::opt::no_grad(fx);
        }
      };
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      std::tuple<double,std::vector<double>,std::vector<double>> lambda_rangefind(Output_deker& sol) const
      {
        double lambda_max = control.lambda_init_max;
        double stepsize = control.lambda_init_stepsize;
        int iter = 1;
        std::vector<double> sampled_lambda;
        std::vector<double> sampled_BIC;
        /////////////////////////////////////////////////////////////
        while(true){
          Output_deker new_sol;
          new_sol.lambda_regularization = stepsize * lambda_max;
          inner_full_opt(new_sol.lambda_regularization,new_sol);
          std::cerr<<"  lambda "<<new_sol.lambda_regularization<<" BIC "<<new_sol.BIC<<"\n";
          //
          if(new_sol.return_status==2){
            lambda_max = new_sol.lambda_regularization;
          }else if(new_sol.feature_index.size()>0){
            if(sampled_lambda.size()==0 || new_sol.BIC < sol.BIC){
              sol(new_sol);
            }
            stepsize += (1-stepsize)*.5; //reduce step size - gone too far
            sampled_lambda.push_back(new_sol.lambda_regularization);
            sampled_BIC.push_back(new_sol.BIC);
          }
          if(new_sol.lambda_regularization<control.lambda_init_min||stepsize>=.95){
            break;
          }
          if(iter==control.lambda_init_maxiter){
            break;
          }
          iter++;
        }
        return std::make_tuple(lambda_max,sampled_lambda,sampled_BIC);
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
    public: 
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      Solve_deker(const Eigen::Ref<const Eigen::MatrixXd> data_matrix,
                  const std::vector<unsigned>& predictors_use_index,
                  const unsigned& which_response,
                  const Opt_Param& control): 
      input_data(data_matrix,predictors_use_index,which_response),
      response_labels(input_data.response,input_data.response_sort_index),
      control(control){}
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      std::tuple<double,double,double,double> calc_model_fit(const double& sigma_kernel_width,
                                                             const Eigen::Ref<const Eigen::VectorXd> v,
                                                             const std::vector<unsigned>& feature_index) const
      {
        Eigen::MatrixXd tmp_z_cache = std::get<0>(math::calc_z_cache(input_data,response_labels,sigma_kernel_width,v,feature_index,false,false));
        //calculate properties = flat, samples in blocks with thresholds in blocks
        Eigen::VectorXd predicted = (tmp_z_cache.array().colwise()*(v.array().square()/v.array().square().sum())).colwise().sum();
        //do Platt scaling to get probabilities
        Eigen::Map<const Eigen::ArrayXd> lab_array(response_labels.y_labels.data(),response_labels.y_labels.size());
        Platt_Params platt_params;
        double A, B;
        std::tie(A,B) = platt_scaling(predicted,lab_array,platt_params);
        Eigen::VectorXd probs = rescale_deci_by_platt(predicted,A,B);
        //calculate degrees of freedom
        Eigen::Map<Eigen::MatrixXd> prob_mat(probs.data(),input_data.response_sort_index.size(),response_labels.y_labels.cols());
        Eigen::MatrixXd pseudohat = prob_mat*response_labels.y_labels_pseudoinv;
        //like taking the trace of the true hat matrix, but uses any value where the pseudoident (y*y^cross) is non-zero
        double df = (response_labels.y_labels_pseudoident.array() > 1e-10).select(pseudohat,0).sum();
        if(df < 0){df = std::numeric_limits<double>::max();}
        //calculate log likelihood
        double logli = (prob_mat.array()*response_labels.y_labels.array()+(1-prob_mat.array())*(1-response_labels.y_labels.array())).log().sum();
        logli /= (double) response_labels.y_labels.cols();
        return std::make_tuple(df,logli,A,B);
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      double opt_sigma(const Eigen::Ref<const Eigen::VectorXd> v,
                       const std::vector<unsigned>& feature_index) const
      {
        Eigen::MatrixXd temp_kernel = math::calc_kernel_matrix(input_data,0,v,feature_index);
        double max_el = -temp_kernel.minCoeff();
        //smallest_diff will be 10 or less, but greater than 0
        double smallest_diff = (temp_kernel.array() == -max_el).select(10.0,temp_kernel.array()+max_el).minCoeff();
        double min_sigma = .5*(log(smallest_diff/2.0)-log(1e5)); //sigma where the smallest difference is very large (only nearest neighbor)
        double max_sigma = .5*(log((max_el)/2.0)-log(1e-5)); //sigma where the largest difference is very small (global average)
        //find best sigma given v
        sigma_opt_prob sigma_opt_tmp(min_sigma,max_sigma,v,feature_index,this);
        Eigen::VectorXd sigma_init = Eigen::VectorXd::Constant(1,sigma_opt_tmp.convert_to_bounded((max_sigma+min_sigma)/2.0));
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
        Eigen::VectorXd transform_v = (v.array().square()/v_opt_tmp.w_transform.array()).sqrt();
        Eigen::VectorXd v_max_temp = limbo_lbfgs_solver(v_opt_tmp,transform_v,false);
        //calculate the next v by transforming with w_transform
        v_max_temp = (v_max_temp.array().square()*v_opt_tmp.w_transform.array()).sqrt();
        return v_max_temp;
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
        for(unsigned i = 1; i<v_i.cols(); i++){
          //find best sigma given v
          v_i(0,i) = opt_sigma(v_i.block(1,i-1,sol.v.size(),1),sol.feature_index);
          //find best v given sigma
          v_i.block(1,i,sol.v.size(),1) = opt_v(lambda_regularization,v_i(0,i),v_i.block(1,i-1,sol.v.size(),1),sol.feature_index);
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
        //set lambda in container
        sol.lambda_regularization = lambda_regularization;
        //initial values for v & sigma
        sol.feature_index = input_data.predictors_index;
        sol.v = Eigen::VectorXd::Constant(sol.feature_index.size(),1);
        sol.sigma_kernel_width = 5;
        sol.return_status = 0;
        sol.BIC = 0;
        ///
        int iter = 0;
        for(; iter<control.w_maxiter; iter++){
          steff_iter_v(lambda_regularization,sol);
          if(sol.return_status != 0) break;
        }
        if(iter >= control.w_maxiter){
          sol.return_status = 3;
        }
        if(sol.return_status != 2){ //feature_index.size() > 0
          double df,lf_error,A,B;
          std::tie(df,lf_error,A,B) = calc_model_fit(sol.sigma_kernel_width,
                                                     sol.v,
                                                     sol.feature_index);
          sol.BIC = log((double)input_data.response.size())*df-2*lf_error;
        }
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      void opt_lambda(Output_deker& sol) const
      {
        double lambda_max;
        std::vector<double> sampled_lambda,sampled_BIC;
        std::tie(lambda_max,sampled_lambda,sampled_BIC) = lambda_rangefind(sol);
        std::cerr<<"  rangefind done\n";
        //typedefs for BO machinery
        using kernel_t = limbo::kernel::Exp<Limbo_Params>;
        using mean_t = limbo::mean::Data<Limbo_Params>;
        using hopt_t = limbo::model::gp::KernelLFOpt<Limbo_Params,limbo::opt::Rprop<Limbo_Params>>;
        using model_t = limbo::model::GP<Limbo_Params, kernel_t, mean_t, hopt_t>;
        using acqui_optimizer_t = limbo::opt::NLOptNoGrad<Limbo_Params>;
        using acqui_function_t = limbo::acqui::GP_UCB<Limbo_Params, model_t>;
        //BO machinery
        model_t bo_model(1,1);
        //
        FirstElem afun;
        
        acqui_optimizer_t acqui_optimizer;
        //put lambda into a vector, divide by max (so it's in [0,1])
        Eigen::VectorXd lambda_vec(1);
        Eigen::VectorXd BIC_vec(1);
        //
        for(unsigned i = 0; i<sampled_lambda.size();i++){
          lambda_vec(0) = sampled_lambda.at(i)/lambda_max;
          BIC_vec(0) = -sampled_BIC.at(i);
          bo_model.add_sample(lambda_vec,BIC_vec);
        } 
        bo_model.optimize_hyperparams();
        ///////////////////////////////////////////////////////////////////////////////////////
        for(unsigned i = 0; i<control.lambda_bo_maxiter; i++){
          acqui_function_t acqui(bo_model, i+sampled_lambda.size());
          auto acqui_optimization = [&](const Eigen::VectorXd& x, bool g) { return acqui(x, afun, g); };
          Eigen::VectorXd starting_point = limbo::tools::random_vector(1, true);
          Eigen::VectorXd new_sample = acqui_optimizer(acqui_optimization, starting_point, true);
          Output_deker new_sol;
          new_sol.lambda_regularization = new_sample(0)*lambda_max;
          //
          inner_full_opt(new_sol.lambda_regularization,new_sol);
          std::cerr<<"  lambda "<<new_sol.lambda_regularization<<" BIC "<<new_sol.BIC<<"\n";
          //set objective to lowest possible if all features are removed
          if(new_sol.return_status == 2){
            new_sol.BIC = std::numeric_limits<double>::max();
          }else{
            if(new_sol.BIC<sol.BIC){
              sol(new_sol);
            }
          }
          //build for next iteration
          //divide by max so it's in [0,1]
          lambda_vec(0) = new_sol.lambda_regularization/lambda_max;
          BIC_vec(0) = -new_sol.BIC;
          //
          bo_model.add_sample(lambda_vec,BIC_vec);
          bo_model.optimize_hyperparams();
        }
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      const Response_Labels& get_response_labels() const {return response_labels;}
      const Split_Data& get_input_data() const {return input_data;}
    };
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
  }
}
#endif