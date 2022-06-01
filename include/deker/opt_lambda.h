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
#ifndef DEKER_FIT_OPT_LAMBDA
#define DEKER_FIT_OPT_LAMBDA
///////////////////////////////////////////
#include <deker/fit/output.h>
#include <deker/fit/opt_param.h>
#include <deker/fit/limbo_params.h>
#include <fstream>
///////////////////////////////////////////
namespace deker{
  namespace fit{
    template <typename ClassT, void(ClassT::*ClassFunc)(const double&,Output_deker&) const>
    class deker_Lambda_Optimizer{
    protected:
      //control parameters
      double lambda_max; //will be updated while rangefinding
      double lambda_min;
      double stepsize;
      int rangefind_maxiter;
      int fill_maxiter;
      double fill_mindist;
      int bo_maxiter;
      //status checks
      int rangefind_iter;
      bool rangefind_complete;
      int fill_iter;
      bool fill_complete;
      int bo_iter;
      bool bo_complete;
      //sampling history
      std::vector<double> sampled_lambda;
      std::vector<double> sampled_BIC;
      //output history
      std::vector<Output_deker> sample_sols;
      //best solution
      unsigned which_best;
      //Output_deker best_sol;
      double null_BIC;
      //BO machinery
      //model typedefs
      using kernel_t = limbo::kernel::Exp<Limbo_Params>;
      using mean_t = limbo::mean::Data<Limbo_Params>;
      using hopt_t = limbo::model::gp::KernelLFOpt<Limbo_Params,limbo::opt::Rprop<Limbo_Params>>;
      using model_t = limbo::model::GP<Limbo_Params, kernel_t, mean_t, hopt_t>;
      //aquisition functions typedefs 
      using acqui_optimizer_t = limbo::opt::NLOptNoGrad<Limbo_Params>;
      using acqui_function_t = limbo::acqui::GP_UCB<Limbo_Params, model_t>;
      //
      model_t bo_model;
      acqui_optimizer_t acqui_optimizer;
      FirstElem afun;
    public:
      deker_Lambda_Optimizer(){}
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      deker_Lambda_Optimizer(const Opt_Param& control)
      {
        build(control);
      }
      void build(const Opt_Param& control){
        lambda_max = control.lambda_init_max;
        lambda_min = control.lambda_init_min;
        stepsize = control.lambda_init_stepsize;
        rangefind_maxiter = control.lambda_init_maxiter;
        fill_mindist = control.lambda_fill_mindist;
        fill_maxiter = control.lambda_fill_maxiter;
        bo_maxiter = control.lambda_bo_maxiter;
        //reset bo_model
        model_t new_bo_model(1,1);
        bo_model = new_bo_model;
        ///
        rangefind_iter = 0;
        if(rangefind_maxiter>0){
          rangefind_complete = false;
        }else{
          rangefind_complete = true;
        }
        ///
        fill_iter = 0;
        if(fill_maxiter>0){
          fill_complete = false;
        }else{
          fill_complete = true;
        }
        ///
        bo_iter = 0;
        if(bo_maxiter>0){
          bo_complete = false;
        }else{
          bo_complete = true;
        }
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      bool opt_iteration(const ClassT& opt_class){
        Output_deker new_sol;
        if(!rangefind_complete){
          /////////////////////////
          //do a rangefind iteration
          new_sol.lambda_regularization = stepsize * lambda_max;
          //Run the model at the selected lambda
          std::cerr<<"check lambda: "<<new_sol.lambda_regularization<<"\n";
          (opt_class.*ClassFunc)(new_sol.lambda_regularization,new_sol);
          sample_sols.push_back(new_sol);
          std::cerr<<"(rf)  lambda "<<new_sol.lambda_regularization<<" sigma "<<new_sol.sigma_kernel_width<<" BIC "<<new_sol.BIC<<" "<<new_sol.return_status<<"\n";
          //Output cases
          if(new_sol.return_status==2){
            lambda_max = new_sol.lambda_regularization;
          }else if(new_sol.feature_index.size()>0){
            if(sampled_lambda.size()==0 || sample_sols.at(which_best).return_status==-1 || new_sol.BIC < sample_sols.at(which_best).BIC){
              which_best = sample_sols.size()-1;
            }
            stepsize += (1-stepsize)*.5; //reduce step size - gone too far
            sampled_lambda.push_back(new_sol.lambda_regularization);
            sampled_BIC.push_back(new_sol.BIC);
          }
          rangefind_iter++;
          if(new_sol.lambda_regularization<lambda_min||stepsize>=.95||rangefind_iter>=rangefind_maxiter){
            rangefind_complete=true;
          }
          /////////////////////////
        }else if(!fill_complete){
          //do a fill iteration
          double max_dist = 0;
          if(fill_iter==0){
            new_sol.lambda_regularization = 0;
            max_dist = 1;
          }else{
            ////////////////////////
            std::vector<double> sort_lambda = sampled_lambda;
            std::sort(sort_lambda.begin(),sort_lambda.end());
            sort_lambda.push_back(lambda_max);
            double max_dist_midpt=0;
            for(unsigned i=1; i<sort_lambda.size();i++){
              double tmp_diff=(sort_lambda.at(i)-sort_lambda.at(i-1))/lambda_max;
              if(tmp_diff>max_dist){
                max_dist=tmp_diff;
                max_dist_midpt=(sort_lambda.at(i)+sort_lambda.at(i-1))/2;
              }
            }
            new_sol.lambda_regularization = max_dist_midpt;
            
            ////////////////////////
            // double max_good_sample_lambda = *max_element(sampled_lambda.begin(), sampled_lambda.end());
            // max_dist = (lambda_max - max_good_sample_lambda)/lambda_max;
            // double max_dist_mid_pt = (lambda_max+max_good_sample_lambda)/2;
            // //
            // if(sampled_lambda.size()>1){
            //   for(unsigned i = 0; i<sampled_lambda.size();i++){
            //     bool init = true;
            //     double min_dist; 
            //     double min_dist_lambda;
            //     for(unsigned j = 0; j<sampled_lambda.size();j++){
            //       if(i!=j){
            //         double dist = std::abs(sampled_lambda.at(i)-sampled_lambda.at(j))/lambda_max;
            //         if(init||dist<min_dist){
            //           init = false;
            //           min_dist = dist;
            //           min_dist_lambda = sampled_lambda.at(j);
            //         }
            //       }
            //     }
            //     if(min_dist>max_dist){
            //       max_dist = min_dist;
            //       max_dist_mid_pt = (sampled_lambda.at(i)+min_dist_lambda)/2;
            //     }
            //   }
            // }
            // new_sol.lambda_regularization = max_dist_mid_pt;
          }
          //Run the model at the selected lambda
          (opt_class.*ClassFunc)(new_sol.lambda_regularization,new_sol);
          sample_sols.push_back(new_sol);
          std::cerr<<"(fi)  lambda "<<new_sol.lambda_regularization<<" sigma "<<new_sol.sigma_kernel_width<<" BIC "<<new_sol.BIC<<" "<<new_sol.return_status<<"\n";
          //Output cases
          if(new_sol.return_status == 2){
            new_sol.BIC = std::numeric_limits<double>::max();
          }else{
            if(sampled_lambda.size()==0 || sample_sols.at(which_best).return_status==-1 || new_sol.BIC < sample_sols.at(which_best).BIC){
              which_best = sample_sols.size()-1;
            }
          }
          sampled_lambda.push_back(new_sol.lambda_regularization);
          sampled_BIC.push_back(new_sol.BIC);
          fill_iter++;
          if(max_dist<=fill_mindist||fill_iter>=fill_maxiter){
            std::cerr<<max_dist<<"<="<<fill_mindist<<" or "<<fill_iter<<">="<<fill_maxiter<<"\n";
            fill_complete=true;
          }
        }else if(!bo_complete){
          /////////////////////////
          //do a bayesian optimization iteration
          Eigen::VectorXd lambda_vec(1);
          Eigen::VectorXd BIC_vec(1);
          if(bo_iter==0){
            /////////////////////////
            //first bo iteration - prepare bo_model
            for(unsigned i = 0; i<sampled_lambda.size();i++){
              lambda_vec(0) = sampled_lambda.at(i)/lambda_max;
              BIC_vec(0) = -sampled_BIC.at(i);
              bo_model.add_sample(lambda_vec,BIC_vec);
            } 
            bo_model.optimize_hyperparams();
            /////////////////////////
          }
          acqui_function_t acqui(bo_model, sampled_lambda.size());
          auto acqui_optimization = [&](const Eigen::VectorXd& x, bool g) { return acqui(x, afun, g); };
          Eigen::VectorXd starting_point = limbo::tools::random_vector(1, true);
          Eigen::VectorXd new_sample = acqui_optimizer(acqui_optimization, starting_point, true);
          new_sol.lambda_regularization = new_sample(0)*lambda_max;
          //Run the model at the selected lambda
          (opt_class.*ClassFunc)(new_sol.lambda_regularization,new_sol);
          sample_sols.push_back(new_sol);
          std::cerr<<"(bo)  lambda "<<new_sol.lambda_regularization<<" sigma "<<new_sol.sigma_kernel_width<<" BIC "<<new_sol.BIC<<" "<<new_sol.return_status<<"\n";
          //Output cases
          if(new_sol.return_status == 2){
            new_sol.BIC = std::numeric_limits<double>::max();
          }else{
            if(sampled_lambda.size()==0 || sample_sols.at(which_best).return_status==-1 || new_sol.BIC < sample_sols.at(which_best).BIC){
              //best_sol(new_sol);
              which_best = sample_sols.size()-1;
            }
          }
          sampled_lambda.push_back(new_sol.lambda_regularization);
          sampled_BIC.push_back(new_sol.BIC);
          //
          lambda_vec(0) = new_sol.lambda_regularization/lambda_max;
          BIC_vec(0) = -new_sol.BIC;
          bo_model.add_sample(lambda_vec,BIC_vec);
          bo_model.optimize_hyperparams();
          //
          bo_iter++;
          if(bo_iter>=bo_maxiter){
            bo_complete=true;
          }
          /////////////////////////
        }
        return !(rangefind_complete&&fill_complete&&bo_complete);
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      std::tuple<bool,bool,bool> get_status() const {return std::make_tuple(rangefind_complete,fill_complete,bo_complete);}
      std::tuple<unsigned,unsigned,unsigned> get_iters() const {return std::make_tuple(rangefind_iter,fill_iter,bo_iter);}
      const std::vector<double>& get_sampled_lambda() const {return sampled_lambda;}
      const std::vector<double>& get_sampled_BIC() const {return sampled_BIC;}
      double get_lambda_max() const {return lambda_max;}
      unsigned get_which_best() const {return which_best;}
      const Output_deker& get_best_sol() const {return sample_sols.at(which_best);}
      const std::vector<Output_deker>& get_sample_output() const {return sample_sols;}
      void set_null_BIC(const double& in_null_BIC){null_BIC=in_null_BIC;}
      double get_null_BIC(){return null_BIC;}
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      void serialize(std::ifstream& infile)
      {
        sampled_lambda.clear();
        sampled_BIC.clear();
        sample_sols.clear();
        model_t new_bo_model(1,1);
        bo_model = new_bo_model;
        ///
        infile.read((char*) (&lambda_max),sizeof(lambda_max));
        infile.read((char*) (&lambda_min),sizeof(lambda_min));
        infile.read((char*) (&stepsize),sizeof(stepsize));
        infile.read((char*) (&rangefind_maxiter),sizeof(rangefind_maxiter));
        infile.read((char*) (&fill_maxiter),sizeof(fill_maxiter));
        infile.read((char*) (&fill_mindist),sizeof(fill_mindist));
        infile.read((char*) (&bo_maxiter),sizeof(bo_maxiter));
        infile.read((char*) (&rangefind_iter),sizeof(rangefind_iter));
        infile.read((char*) (&rangefind_complete),sizeof(rangefind_complete));
        infile.read((char*) (&fill_iter),sizeof(fill_iter));
        infile.read((char*) (&fill_complete),sizeof(fill_complete));
        infile.read((char*) (&bo_iter),sizeof(bo_complete));
        infile.read((char*) (&bo_complete),sizeof(bo_complete));
        infile.read((char*) (&null_BIC),sizeof(null_BIC));
        ///
        size_t n_samples;
        double tmp;
        infile.read((char*) (&n_samples),sizeof(n_samples));
        for(size_t i = 0; i<n_samples;i++){
          infile.read((char*) (&tmp),sizeof(tmp));
          sampled_lambda.push_back(tmp);
          infile.read((char*) (&tmp),sizeof(tmp));
          sampled_BIC.push_back(tmp);
        }
        //
        infile.read((char*) (&n_samples),sizeof(n_samples));
        Output_deker temp;
        for(size_t i = 0; i<n_samples;i++){
          temp.serialize(infile);
          sample_sols.push_back(temp);
        }
        infile.read((char*) (&which_best),sizeof(which_best));
        //
        //rebuild bo_model if necessary
        if(rangefind_complete&&fill_complete&&!bo_complete&&bo_iter>0){
          Eigen::VectorXd lambda_vec(1);
          Eigen::VectorXd BIC_vec(1);
          for(unsigned i = 0; i<sampled_lambda.size();i++){
            lambda_vec(0) = sampled_lambda.at(i)/lambda_max;
            BIC_vec(0) = -sampled_BIC.at(i);
            bo_model.add_sample(lambda_vec,BIC_vec);
          } 
          bo_model.optimize_hyperparams();
        }
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      void serialize(std::ofstream& outfile) const
      {
        outfile.write((char*) (&lambda_max),sizeof(lambda_max));
        outfile.write((char*) (&lambda_min),sizeof(lambda_min));
        outfile.write((char*) (&stepsize),sizeof(stepsize));
        outfile.write((char*) (&rangefind_maxiter),sizeof(rangefind_maxiter));
        outfile.write((char*) (&fill_maxiter),sizeof(fill_maxiter));
        outfile.write((char*) (&fill_mindist),sizeof(fill_mindist));
        outfile.write((char*) (&bo_maxiter),sizeof(bo_maxiter));
        outfile.write((char*) (&rangefind_iter),sizeof(rangefind_iter));
        outfile.write((char*) (&rangefind_complete),sizeof(rangefind_complete));
        outfile.write((char*) (&fill_iter),sizeof(fill_iter));
        outfile.write((char*) (&fill_complete),sizeof(fill_complete));
        outfile.write((char*) (&bo_iter),sizeof(bo_complete));
        outfile.write((char*) (&bo_complete),sizeof(bo_complete));
        outfile.write((char*) (&null_BIC),sizeof(null_BIC));
        ///
        size_t n_samples = sampled_lambda.size();
        outfile.write((char*) (&n_samples),sizeof(n_samples));
        for(size_t i = 0; i<n_samples;i++){
          outfile.write((char*) (&sampled_lambda.at(i)),sizeof(double));
          outfile.write((char*) (&sampled_BIC.at(i)),sizeof(double));
        }
        ///
        n_samples = sample_sols.size();
        outfile.write((char*) (&n_samples),sizeof(n_samples));
        for(size_t i = 0; i<n_samples;i++){
          sample_sols.at(i).serialize(outfile);
        }
        outfile.write((char*) (&which_best),sizeof(which_best));
      }
    };
  }
}


#endif