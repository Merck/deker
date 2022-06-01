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
#ifndef DEKER_FIT_PLATT
#define DEKER_FIT_PLATT
///////////////////////////////////////////
#include <Eigen/Dense>
#include <iostream>
///////////////////////////////////////////
namespace deker{
  namespace fit{
    struct Platt_Params{
      unsigned maxiter; //maximum number of iterations
      double minstep; //minimum step taken in line search
      double sigma; //set to any value > 0
      Platt_Params(){
        //default values
        maxiter = 100; 
        minstep = 1e-10;
        sigma = 1e-12;
      }
    };
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    //implementation based on Lin et al. 2007 "A Note on Platt's Probabilistic Outputs for Support Vector Machines"
    //Appendix C Pseudocode
    std::tuple<double,double> platt_scaling(const Eigen::Ref<const Eigen::ArrayXd> deci,
                                            const Eigen::Ref<const Eigen::ArrayXd> label,
                                            const Platt_Params& params)
    {
      unsigned prior1 = label.sum();
      unsigned prior0 = label.size()-prior1;
      //make transformed target values 
      double hiTarget = (prior1+1.0)/(prior1+2.0);
      double loTarget = 1/(prior0+2.0);
      Eigen::ArrayXd t = label*hiTarget+(1-label)*loTarget;
      //
      double A = 0.0;
      double B = std::log(((double)prior0+1.0)/((double)prior1+1.0));
      double fval = 0.0;
      //departure from pseudocode, as we know the value of fApB for all values with A = 0
      double h11 = params.sigma;
      double h22 = params.sigma;
      double h21 = 0.0;
      double g1 = 0.0;
      double g2 = 0.0;
      //
      {
        double p, q;
        if(B>=0){
          fval = (t*B).sum()+std::log(1.0+exp(-B))*t.size();
          p = exp(-B)/(1.0+exp(-B));
          q = 1.0/(1.0+exp(-B));
        }else{
          fval = ((t-1)*B).sum()+std::log(1.0+exp(B))*t.size();
          p = 1.0/(1.0+exp(B));
          q = exp(B)/(1.0+exp(B));
        }
        double d2 = p*q;
        h11 = (deci.square()).sum()*d2;
        h22 = d2*deci.size();
        h21 = deci.sum()*d2;
        Eigen::ArrayXd d1 = t-p;
        g1 = (deci*d1).sum();
        g2 = d1.sum();
      }
      //
      if(std::isnan(g1)||std::isnan(g2)){
        throw std::logic_error("nans entered platt scaling calculations");
      }
      double stepsize;
      for(unsigned it = 0; it<params.maxiter;it++){
        //Compute modified Newton directions
        double det = h11*h22-h21*h21;
        double dA = -(h22*g1-h21*g2)/det;
        double dB = -(-h21*g1+h11*g2)/det;
        double gd = g1*dA+g2*dB;
        stepsize = 1;
        while(stepsize >= params.minstep){ //Line search
          double newA=A+stepsize*dA;
          double newB=B+stepsize*dB;
          double newf=0.0;
          for(unsigned i = 0; i<deci.size();i++){
            double fApB=deci(i)*newA+newB;
            if(fApB >= 0){
              newf+=t(i)*fApB+std::log(1.0+exp(-fApB));
            }else{
              newf+=(t(i)-1)*fApB+std::log(1.0+exp(fApB));
            }
          }
          if(newf<fval+.0001*stepsize*gd){
            A=newA;
            B=newB;
            fval=newf;
            break; //Sufficient decrease satisfied
          }else{
            stepsize /= 2.0;
          }
        }
        if(stepsize < params.minstep){
          if(std::isnan(g1)||std::isnan(g2)){
            throw std::logic_error("nans entered platt scaling calculations");
          }
          break;
        }
        //Update Gradient and Hessian (use H' = H + sigma I)
        h11 = params.sigma;
        h22 = params.sigma;
        h21 = 0.0;
        g1 = 0.0;
        g2 = 0.0;
        for(unsigned i=0; i<deci.size();i++){
          double fApB = deci(i)*A+B;
          double p,q;
          if(fApB >= 0){
            p = exp(-fApB)/(1.0+exp(-fApB));
            q = 1.0/(1.0+exp(-fApB));
          }else{
            p = 1.0/(1.0+exp(fApB));
            q = exp(fApB)/(1.0+exp(fApB));
          }
          double d2 = p*q;
          h11+=deci(i)*deci(i)*d2;
          h22+=d2;
          h21+=deci(i)*d2;
          double d1=t(i)-p;
          g1+=deci(i)*d1;
          g2+=d1;
        }
        //diverged slightly from publication here
        //reduced to 1.5e-5 from 1e-5; was getting a lot of failure errors where g1/g2 were stuck >1e-5 but <1.5e-5
        if(std::abs(g1)<1.5e-5&&std::abs(g2)<1.5e-5){
          break;
        }
        if(it==params.maxiter-1){
          std::cerr<<"platt_scaling reached maximum iterations.\n";
        }
      }
      return(std::make_tuple(A,B));
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    Eigen::VectorXd rescale_deci_by_platt(const Eigen::Ref<const Eigen::ArrayXd> deci, const double& A, const double& B)
    {
      return 1.0/(1.0+(A*deci+B).exp());
    }
  }
}
#endif
