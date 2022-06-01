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
#ifndef DEKER_FIT_OPT_PARAM
#define DEKER_FIT_OPT_PARAM
///////////////////////////////////////////
#include <fstream>
///////////////////////////////////////////
namespace deker{
  namespace fit{
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    struct Opt_Param{
      double lambda_init_max;
      double lambda_init_min;
      double lambda_init_stepsize;
      int lambda_init_maxiter;
      double lambda_fill_mindist;
      int lambda_fill_maxiter;
      int lambda_bo_maxiter;
      double h_huber_loss;
      double hinge_max;
      double w_convergence_threshold;
      double w_maxiter;
      double feature_drop_threshold;
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      //constructor with default values
      Opt_Param(){
        lambda_init_max = .25;
        lambda_init_min = .005;
        lambda_init_stepsize = .5;
        lambda_init_maxiter = 100;
        lambda_fill_mindist = .2;
        lambda_fill_maxiter = 10;
        lambda_bo_maxiter = 20;
        h_huber_loss = .001;
        hinge_max = 1;
        w_convergence_threshold = 1e-7;
        w_maxiter = 100;
        feature_drop_threshold = .001; //this is applied to v, so should be set to the square root of what seems like a reasonable value 
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      void serialize(std::ifstream& infile)
      {
        //read in from stream
        infile.read((char*) (&lambda_init_max),sizeof(lambda_init_max));
        infile.read((char*) (&lambda_init_min),sizeof(lambda_init_min));
        infile.read((char*) (&lambda_init_stepsize),sizeof(lambda_init_stepsize));
        infile.read((char*) (&lambda_init_maxiter),sizeof(lambda_init_maxiter));
        infile.read((char*) (&lambda_fill_mindist),sizeof(lambda_fill_mindist));
        infile.read((char*) (&lambda_fill_maxiter),sizeof(lambda_fill_maxiter));
        infile.read((char*) (&lambda_bo_maxiter),sizeof(lambda_bo_maxiter));
        infile.read((char*) (&h_huber_loss),sizeof(h_huber_loss));
        infile.read((char*) (&hinge_max),sizeof(hinge_max));
        infile.read((char*) (&w_convergence_threshold),sizeof(w_convergence_threshold));
        infile.read((char*) (&w_maxiter),sizeof(w_maxiter));
        infile.read((char*) (&feature_drop_threshold),sizeof(feature_drop_threshold));
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      void serialize(std::ofstream& outfile) const
      {
        outfile.write((char*) (&lambda_init_max),sizeof(lambda_init_max));
        outfile.write((char*) (&lambda_init_min),sizeof(lambda_init_min));
        outfile.write((char*) (&lambda_init_stepsize),sizeof(lambda_init_stepsize));
        outfile.write((char*) (&lambda_init_maxiter),sizeof(lambda_init_maxiter));
        outfile.write((char*) (&lambda_fill_mindist),sizeof(lambda_fill_mindist));
        outfile.write((char*) (&lambda_fill_maxiter),sizeof(lambda_fill_maxiter));
        outfile.write((char*) (&lambda_bo_maxiter),sizeof(lambda_bo_maxiter));
        outfile.write((char*) (&h_huber_loss),sizeof(h_huber_loss));
        outfile.write((char*) (&hinge_max),sizeof(hinge_max));
        outfile.write((char*) (&w_convergence_threshold),sizeof(w_convergence_threshold));
        outfile.write((char*) (&w_maxiter),sizeof(w_maxiter));
        outfile.write((char*) (&feature_drop_threshold),sizeof(feature_drop_threshold));
      }
    };
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
  }
}
#endif