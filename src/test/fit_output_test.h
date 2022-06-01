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
#include <deker/io/output.h>
#include <gtest/gtest.h>
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
TEST (OutputTest, OutputWriteRead){
  deker::fit::Output_deker_IO ref_test;
  //Fill with fake data
  ref_test.lambda_regularization = 1.0525;
  ref_test.sigma_kernel_width = 500.1230;
  ref_test.BIC = 239.1123;
  ref_test.null_BIC = 654321.012;
  ref_test.return_status = 0;
  ref_test.data_filename = "dummy_data_location.bin";
  ref_test.which_response = 100;
  ref_test.response_name = "dummy_response_name";
  //
  ref_test.feature_index.push_back(0);
  ref_test.feature_index.push_back(1);
  ref_test.feature_index.push_back(2);
  ref_test.feature_index.push_back(4);
  ref_test.key_to_original_data.push_back(0);
  ref_test.key_to_original_data.push_back(1);
  ref_test.key_to_original_data.push_back(2);
  ref_test.key_to_original_data.push_back(4);
  ref_test.key_to_original_data.push_back(8);
  ref_test.key_to_original_data.push_back(9);
  ref_test.key_to_original_data.push_back(10);
  ref_test.feature_names.push_back("A0");
  ref_test.feature_names.push_back("B1");
  ref_test.feature_names.push_back("C2");
  ref_test.feature_names.push_back("E4");
  ref_test.v.resize(4);
  for(unsigned i = 0; i<ref_test.v.size(); i++){
    ref_test.v(i) = -i+.005*i;
  }
  //location to write test file to
  std::string write_test_loc = "write_test_loc.dout";
  deker::io::write_to_binary<deker::fit::Output_deker_IO>(write_test_loc,ref_test);
  //read from location, building new control file 
  deker::fit::Output_deker_IO bin_wr; 
  deker::io::read_from_binary<deker::fit::Output_deker_IO>(write_test_loc,bin_wr);
  //compare values in ref_test to bin_wr
  ASSERT_EQ(ref_test.lambda_regularization,bin_wr.lambda_regularization)<<"lambda doesn't match";
  ASSERT_EQ(ref_test.sigma_kernel_width,bin_wr.sigma_kernel_width)<<"sigma doesn't match";
  ASSERT_EQ(ref_test.BIC,bin_wr.BIC)<<"BIC doesn't match";
  ASSERT_EQ(ref_test.null_BIC,bin_wr.null_BIC)<<"null_BIC doesn't match";
  ASSERT_EQ(ref_test.return_status,bin_wr.return_status)<<"return_status doesn't match";
  ASSERT_EQ(ref_test.data_filename,bin_wr.data_filename)<<"data_filename doesn't match";
  ASSERT_EQ(ref_test.which_response,bin_wr.which_response)<<"which_response doesn't match";
  ASSERT_EQ(ref_test.response_name,bin_wr.response_name)<<"response_name doesn't match";
  //
  ASSERT_EQ(ref_test.v.size(),bin_wr.v.size())<<"v of unequal sizes";
  for(unsigned i = 0;i<ref_test.v.size();i++)
    ASSERT_EQ(ref_test.v(i),bin_wr.v(i))<<"v differs at index "<<i;
  //
  ASSERT_EQ(ref_test.feature_index.size(),bin_wr.feature_index.size())<<"feature_index of unequal sizes";
  for(unsigned i = 0;i<ref_test.feature_index.size();i++)
    ASSERT_EQ(ref_test.feature_index.at(i),bin_wr.feature_index.at(i))<<"feature_index differs at index "<<i;
  //
  ASSERT_EQ(ref_test.key_to_original_data.size(),bin_wr.key_to_original_data.size())<<"key_to_original_data of unequal sizes";
  for(unsigned i = 0;i<ref_test.key_to_original_data.size();i++)
    ASSERT_EQ(ref_test.key_to_original_data.at(i),bin_wr.key_to_original_data.at(i))<<"key_to_original_data differs at index "<<i;
  //
  ASSERT_EQ(ref_test.feature_names.size(),bin_wr.feature_names.size())<<"feature_names of unequal sizes";
  for(unsigned i = 0;i<ref_test.feature_names.size();i++)
    ASSERT_EQ(ref_test.feature_names.at(i),bin_wr.feature_names.at(i))<<"feature_names differs at index "<<i;
  //remove control file
  remove(write_test_loc.c_str());
}
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////