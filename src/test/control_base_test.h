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
#include <deker/io/control.h>
#include <gtest/gtest.h>
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
TEST (ControlBaseTest, ControlBaseTextWriteRead){
  deker::io::Control_Base ref_test;
  //constructors fill Opt_Param with default values
  ref_test.data_file="test_data_loc.csv";
  ref_test.lambda_iter_per_job = 5;
  //fake data for in_res
  ref_test.in_res.push_back(0);
  ref_test.in_res.push_back(1);
  ref_test.in_res.push_back(2);
  //fake data for ex_res
  ref_test.ex_res.push_back(0);
  ref_test.ex_res.push_back(1);
  ref_test.ex_res.push_back(2);
  //fake data for in_pred
  ref_test.in_pred.push_back(0);
  ref_test.in_pred.push_back(1);
  ref_test.in_pred.push_back(2);
  //fake data for ex_pred
  //ref_test.ex_pred.push_back(0);
  //ref_test.ex_pred.push_back(1);
  //ref_test.ex_pred.push_back(2);
  //location to write test file to
  std::string write_test_loc = "write_test_loc.control";
  ref_test.write_control_file(write_test_loc); //writes to text
  //read from location, building new control file 
  deker::io::Control_Base text_wr;
  text_wr.read_from_txt_file(write_test_loc);
  //compare values in ref_test to text_wr
  ASSERT_EQ(ref_test.data_file,text_wr.data_file)<<"data_file doesn't match";
  ASSERT_EQ(ref_test.lambda_iter_per_job,text_wr.lambda_iter_per_job)<<"lambda_iter_per_job doesn't match";
  ASSERT_EQ(ref_test.in_res,text_wr.in_res)<<"in_res doesn't match";
  ASSERT_EQ(ref_test.ex_res,text_wr.ex_res)<<"ex_res doesn't match";
  ASSERT_EQ(ref_test.in_pred,text_wr.in_pred)<<"in_pred doesn't match";
  ASSERT_EQ(ref_test.ex_pred,text_wr.ex_pred)<<"ex_pred doesn't match";
  //opt_param members
  ASSERT_EQ(ref_test.opt_param.h_huber_loss,text_wr.opt_param.h_huber_loss)<<"h_huber_loss doesn't match";
  ASSERT_EQ(ref_test.opt_param.hinge_max,text_wr.opt_param.hinge_max)<<"hinge_max doesn't match";
  ASSERT_EQ(ref_test.opt_param.w_convergence_threshold,text_wr.opt_param.w_convergence_threshold)<<"w_convergence_threshold doesn't match";
  ASSERT_EQ(ref_test.opt_param.w_maxiter,text_wr.opt_param.w_maxiter)<<"w_maxiter doesn't match";
  ASSERT_EQ(ref_test.opt_param.feature_drop_threshold,text_wr.opt_param.feature_drop_threshold)<<"feature_drop_threshold doesn't match";
  //remove control file
  remove(write_test_loc.c_str());
}
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
TEST (ControlBaseTest, ControlBaseBinWriteRead){
  deker::io::Control_Base ref_test;
  //constructors fill Opt_Param/HOpt_Param with default values
  ref_test.data_file="test_data_loc.csv";
  ref_test.lambda_iter_per_job = 5;
  //fake data for in_res
  ref_test.in_res.push_back(0);
  ref_test.in_res.push_back(1);
  ref_test.in_res.push_back(2);
  //fake data for ex_res
  ref_test.ex_res.push_back(0);
  ref_test.ex_res.push_back(1);
  ref_test.ex_res.push_back(2);
  //fake data for in_pred
  ref_test.in_pred.push_back(0);
  ref_test.in_pred.push_back(1);
  ref_test.in_pred.push_back(2);
  //fake data for ex_pred
  //ref_test.ex_pred.push_back(0);
  //ref_test.ex_pred.push_back(1);
  //ref_test.ex_pred.push_back(2);
  //location to write test file to
  std::string write_test_loc = "write_test_loc.control";
  //write to file
  deker::io::write_to_binary<deker::io::Control_Base>(write_test_loc,ref_test);
  //std::ofstream outfile(write_test_loc, std::ios::out | std::ios::binary | std::ios::trunc);
  //ref_test.serialize(outfile);
  //outfile.close();
  //read from location, building new control file 
  deker::io::Control_Base bin_wr;
  deker::io::read_from_binary<deker::io::Control_Base>(write_test_loc,bin_wr);
  //std::ifstream infile(write_test_loc, std::ios::in | std::ios::binary);
  //bin_wr.serialize(infile);
  //infile.close();
  //compare values in ref_test to text_wr
  ASSERT_EQ(ref_test.data_file,bin_wr.data_file)<<"data_file doesn't match";
  ASSERT_EQ(ref_test.lambda_iter_per_job,bin_wr.lambda_iter_per_job)<<"lambda_iter_per_job doesn't match";
  ASSERT_EQ(ref_test.in_res,bin_wr.in_res)<<"in_res doesn't match";
  ASSERT_EQ(ref_test.ex_res,bin_wr.ex_res)<<"ex_res doesn't match";
  ASSERT_EQ(ref_test.in_pred,bin_wr.in_pred)<<"in_pred doesn't match";
  ASSERT_EQ(ref_test.ex_pred,bin_wr.ex_pred)<<"ex_pred doesn't match";
  //opt_param members
  ASSERT_EQ(ref_test.opt_param.h_huber_loss,bin_wr.opt_param.h_huber_loss)<<"h_huber_loss doesn't match";
  ASSERT_EQ(ref_test.opt_param.hinge_max,bin_wr.opt_param.hinge_max)<<"hinge_max doesn't match";
  ASSERT_EQ(ref_test.opt_param.w_convergence_threshold,bin_wr.opt_param.w_convergence_threshold)<<"w_convergence_threshold doesn't match";
  ASSERT_EQ(ref_test.opt_param.w_maxiter,bin_wr.opt_param.w_maxiter)<<"w_maxiter doesn't match";
  ASSERT_EQ(ref_test.opt_param.feature_drop_threshold,bin_wr.opt_param.feature_drop_threshold)<<"feature_drop_threshold doesn't match";
  //remove control file
  remove(write_test_loc.c_str());
}
