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
#ifndef DEKER_FIT_LIMBO_PARAMS
#define DEKER_FIT_LIMBO_PARAMS
///////////////////////////////////////////
#define USE_NLOPT
#include <limbo/acqui/gp_ucb.hpp>
#include <limbo/bayes_opt/boptimizer.hpp>
#include <limbo/kernel/matern_five_halves.hpp>
#include <limbo/mean/data.hpp>
#include <limbo/model/gp.hpp>
#include <limbo/tools/macros.hpp>
#include <limbo/init/random_sampling.hpp>
#include <limbo/model/gp/kernel_lf_opt.hpp>
#include <limbo/opt/nlopt_no_grad.hpp>
#include <limbo/opt/nlopt_grad.hpp>
///////////////////////////////////////////
namespace deker{
  namespace fit{
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    struct PF_Params {
      struct opt_nloptnograd : public limbo::defaults::opt_nloptnograd {
        //BO_PARAM(int, iterations, 500);
        //BO_PARAM(double, fun_tolerance, -1);
        BO_PARAM(double, xrel_tolerance, 1e-2);
      };
    };
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    struct FS_Params {
      struct opt_nloptgrad : public limbo::defaults::opt_nloptgrad {
        //BO_PARAM(int, iterations, 1000);
        //BO_PARAM(double, eps_stop, 0.0);
      };
    };
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    struct Limbo_Params {
      struct opt_nloptnograd : public limbo::defaults::opt_nloptnograd {
      };
      // no noise
      struct kernel : public limbo::defaults::kernel {
        BO_PARAM(double, noise, 1e-10);
      };
      struct kernel_exp : public limbo::defaults::kernel_exp {
      };
      // 10 random samples to initialize the algorithm
      struct init_randomsampling {
        BO_PARAM(int, samples, 10);
      };
      // stop after 40 iterations
      struct stop_maxiterations {
        BO_PARAM(int, iterations, 40);
      };
      // default parameters for acqui_ucb
      struct acqui_gpucb : public limbo::defaults::acqui_gpucb {
        //BO_PARAM(double, delta, .35);
      };
      struct opt_rprop : public limbo::defaults::opt_rprop {
        //BO_PARAM(int, iterations, 1000);
        //BO_PARAM(double, eps_stop, 0.0);
      };
    };
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    struct FirstElem {
      using result_type = double;
      double operator()(const Eigen::VectorXd& x) const
      {
        return x(0);
      }
    };
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
  }
}
///////////////////////////////////////////
#endif
