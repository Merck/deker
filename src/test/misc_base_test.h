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
#include <deker/io/misc.h>
#include <gtest/gtest.h>
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
TEST (MiscBaseTest, Eigen_Matrix_IOTextWriteRead){
  deker::io::Eigen_Matrix_IO ref_test;
  //Add dummy data
  ref_test.matrix.resize(10,10);
  for(unsigned i = 0;i<ref_test.matrix.rows();i++)
    for(unsigned j = 0; j<ref_test.matrix.cols();j++)
      ref_test.matrix(i,j) = i+j;
  ref_test.rows = ref_test.matrix.rows();
  ref_test.cols = ref_test.matrix.cols();
  ref_test.feature_names.push_back("A0"); 
  ref_test.feature_names.push_back("B1"); 
  ref_test.feature_names.push_back("C2"); 
  ref_test.feature_names.push_back("D3"); 
  ref_test.feature_names.push_back("E4"); 
  ref_test.feature_names.push_back("F5"); 
  ref_test.feature_names.push_back("G6"); 
  ref_test.feature_names.push_back("H7"); 
  ref_test.feature_names.push_back("I9"); 
  ref_test.feature_names.push_back("J9"); 
  //location to write test file to
  std::string write_test_loc = "write_test_loc.csv";
  ref_test.write_to_text(write_test_loc); //writes to text
  //read from location as text, building new matrix data
  deker::io::Eigen_Matrix_IO text_wr(write_test_loc,false);
  //compare values in ref_test to text_wr
  ASSERT_TRUE(ref_test.matrix.isApprox(text_wr.matrix))<<"matrix doesn't match";
  ASSERT_EQ(ref_test.rows,text_wr.rows)<<"rows doesn't match";
  ASSERT_EQ(ref_test.cols,text_wr.cols)<<"cols doesn't match";
  ASSERT_EQ(ref_test.feature_names.size(),text_wr.feature_names.size())<<"feature_names of different sizes";
  for(unsigned i = 0; i<ref_test.feature_names.size();i++)
    ASSERT_EQ(ref_test.feature_names.at(i),text_wr.feature_names.at(i))<<"feature_names differs at index "<<i;
  //read from location as text, sizes only 
  deker::io::Eigen_Matrix_IO sizes_only_data(write_test_loc,false,true);
  ASSERT_EQ(sizes_only_data.rows,ref_test.rows)<<"rows doesn't match for sizes_only_data";
  ASSERT_EQ(sizes_only_data.cols,ref_test.cols)<<"rows doesn't match for sizes_only_data";
  //remove test file
  remove(write_test_loc.c_str());
}
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
TEST (MiscBaseTest, Eigen_Matrix_IOBinaryWriteRead){
  deker::io::Eigen_Matrix_IO ref_test;
  //Add dummy data
  ref_test.matrix.resize(10,10);
  for(unsigned i = 0;i<ref_test.matrix.rows();i++)
    for(unsigned j = 0; j<ref_test.matrix.cols();j++)
      ref_test.matrix(i,j) = i+j;
  ref_test.rows = ref_test.matrix.rows();
  ref_test.cols = ref_test.matrix.cols();
  ref_test.feature_names.push_back("A0"); 
  ref_test.feature_names.push_back("B1"); 
  ref_test.feature_names.push_back("C2"); 
  ref_test.feature_names.push_back("D3"); 
  ref_test.feature_names.push_back("E4"); 
  ref_test.feature_names.push_back("F5"); 
  ref_test.feature_names.push_back("G6"); 
  ref_test.feature_names.push_back("H7"); 
  ref_test.feature_names.push_back("I9"); 
  ref_test.feature_names.push_back("J9"); 
  //location to write test file to
  std::string write_test_loc = "write_test_loc.csv";
  ref_test.write_to_binary(write_test_loc); //writes to text
  //read from location as text, building new matrix data
  deker::io::Eigen_Matrix_IO bin_wr(write_test_loc);
  //compare values in ref_test to bin_wr
  ASSERT_TRUE(ref_test.matrix.isApprox(bin_wr.matrix))<<"matrix doesn't match";
  ASSERT_EQ(ref_test.rows,bin_wr.rows)<<"rows doesn't match";
  ASSERT_EQ(ref_test.cols,bin_wr.cols)<<"cols doesn't match";
  ASSERT_EQ(ref_test.feature_names.size(),bin_wr.feature_names.size())<<"feature_names of different sizes";
  for(unsigned i = 0; i<ref_test.feature_names.size();i++)
    ASSERT_EQ(ref_test.feature_names.at(i),bin_wr.feature_names.at(i))<<"feature_names differs at index "<<i;
  deker::io::Eigen_Matrix_IO sizes_only_data(write_test_loc,true,true);
  ASSERT_EQ(sizes_only_data.rows,ref_test.rows)<<"rows doesn't match for sizes_only_data";
  ASSERT_EQ(sizes_only_data.cols,ref_test.cols)<<"rows doesn't match for sizes_only_data";
  //remove test file
  remove(write_test_loc.c_str());
}
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////