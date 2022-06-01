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
#include <deker/fit/split_data.h>
#include <deker/io/misc.h>
#include <gtest/gtest.h>
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
class SplitDataTest : public ::testing::Test{
protected:
    SplitDataTest(){
        input_data_0 = Eigen::MatrixXd::Random(20,5);
        input_data_1 = Eigen::MatrixXd::Random(20,5);
        input_data_2 = Eigen::MatrixXd::Random(20,5);
        input_data_3 = Eigen::MatrixXd::Random(20,5);
        input_data_4 = Eigen::MatrixXd::Random(20,5);
        which_response = 0;
        predictors_use_index.push_back(1);
        predictors_use_index.push_back(2);
        predictors_use_index.push_back(3);
        predictors_use_index.push_back(4);
    }
    Eigen::MatrixXd input_data_0;
    Eigen::MatrixXd input_data_1;
    Eigen::MatrixXd input_data_2;
    Eigen::MatrixXd input_data_3;
    Eigen::MatrixXd input_data_4;
    std::vector<unsigned> predictors_use_index;
    unsigned which_response;
};
TEST_F(SplitDataTest, SortPredictorsTest){
    bool is_sorted;
    
    deker::fit::Split_Data proc_dat_0(input_data_0,which_response,predictors_use_index);
    is_sorted = true;
    for(unsigned i=1; i<proc_dat_0.response.size();i++){
        if(proc_dat_0.response(i-1)>proc_dat_0.response(i)){
            is_sorted=false;
            break;
        }
    }
    ASSERT_TRUE(is_sorted)<<"Split data is not sorted by response (case 0)";
    //
    deker::fit::Split_Data proc_dat_1(input_data_1,which_response,predictors_use_index);
    is_sorted=true;
    for(unsigned i=1; i<proc_dat_1.response.size();i++){
        if(proc_dat_1.response(i-1)>proc_dat_1.response(i)){
            is_sorted=false;
            break;
        }
    }
    ASSERT_TRUE(is_sorted)<<"Split data is not sorted by response (case 1)";
    //
    deker::fit::Split_Data proc_dat_2(input_data_2,which_response,predictors_use_index);
    is_sorted=true;
    for(unsigned i=1; i<proc_dat_2.response.size();i++){
        if(proc_dat_2.response(i-1)>proc_dat_2.response(i)){
            is_sorted=false;
            break;
        }
    }
    ASSERT_TRUE(is_sorted)<<"Split data is not sorted by response (case 2)";
    ///
    deker::fit::Split_Data proc_dat_3(input_data_3,which_response,predictors_use_index);
    is_sorted=true;
    for(unsigned i=1; i<proc_dat_3.response.size();i++){
        if(proc_dat_3.response(i-1)>proc_dat_3.response(i)){
            is_sorted=false;
            break;
        }
    }
    ASSERT_TRUE(is_sorted)<<"Split data is not sorted by response (case 3)";
    ///
    deker::fit::Split_Data proc_dat_4(input_data_4,which_response,predictors_use_index);
    is_sorted=true;
    for(unsigned i=1; i<proc_dat_4.response.size();i++){
        if(proc_dat_4.response(i-1)>proc_dat_4.response(i)){
            is_sorted=false;
            break;
        }
    }
    ASSERT_TRUE(is_sorted)<<"Split data is not sorted by response (case 4)";
}
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
TEST_F(SplitDataTest, SubsetPredictorsTest){
    std::vector<unsigned> subset_predictors_342{3,4,2};
    std::vector<unsigned> subset_predictors_20{2,0};
    //
    deker::fit::Split_Data proc_dat_0(input_data_0,which_response,predictors_use_index);
    Eigen::MatrixXd sorted_data_0 = input_data_0;
    proc_dat_0.subset_predictors(subset_predictors_342);
    //
    ASSERT_EQ(proc_dat_0.predictors.cols(),3)<<"Split_Data predictors has the wrong number of columns post-subsetting";
    ASSERT_EQ(proc_dat_0.key_to_original_data.size(),3)<<"Split_Data key_to_original_data has the wrong size post-subsetting";
    ASSERT_EQ(proc_dat_0.predictors_index.size(),3)<<"Split_Data predictors_index has the wrong size post-subsetting";
    //
    ASSERT_EQ(proc_dat_0.key_to_original_data.at(0),3)<<"Split_Data key_to_original_data.at(0) has the wrong value post-subsetting";
    ASSERT_NEAR((proc_dat_0.predictors.col(0).array()-sorted_data_0.col(3).array()).abs().sum(),0,1e-10)<<"Split_Data predictors(0) post-subsetting does not match data pre-subsetting";
    ASSERT_EQ(proc_dat_0.predictors_index.at(0),0)<<"Split_Data predictors_index.at(0) has the wrong value post-subsetting";
    //
    ASSERT_EQ(proc_dat_0.key_to_original_data.at(1),4)<<"Split_Data key_to_original_data.at(1) has the wrong value post-subsetting";
    ASSERT_NEAR((proc_dat_0.predictors.col(1).array()-sorted_data_0.col(4).array()).abs().sum(),0,1e-10)<<"Split_Data predictors(1) post-subsetting does not match data pre-subsetting";
    ASSERT_EQ(proc_dat_0.predictors_index.at(1),1)<<"Split_Data predictors_index.at(0) has the wrong value post-subsetting";
    //
    ASSERT_EQ(proc_dat_0.key_to_original_data.at(2),2)<<"Split_Data key_to_original_data.at(2) has the wrong value post-subsetting";
    ASSERT_NEAR((proc_dat_0.predictors.col(2).array()-sorted_data_0.col(2).array()).abs().sum(),0,1e-10)<<"Split_Data predictors(2) post-subsetting does not match data pre-subsetting";
    ASSERT_EQ(proc_dat_0.predictors_index.at(2),2)<<"Split_Data predictors_index.at(0) has the wrong value post-subsetting";
    //
    Eigen::MatrixXd subset_once_data_0 = input_data_0;
    proc_dat_0.subset_predictors(subset_predictors_20);
    //
    ASSERT_EQ(proc_dat_0.predictors.cols(),2)<<"Split_Data predictors has the wrong number of columns post-subsetting (2nd iter)";
    ASSERT_EQ(proc_dat_0.key_to_original_data.size(),2)<<"Split_Data key_to_original_data has the wrong size post-subsetting (2nd iter)";
    ASSERT_EQ(proc_dat_0.predictors_index.size(),2)<<"Split_Data predictors_index has the wrong size post-subsetting (2nd iter)";
    //
    ASSERT_EQ(proc_dat_0.key_to_original_data.at(0),2)<<"Split_Data key_to_original_data.at(0) has the wrong value post-subsetting (2nd iter)";
    ASSERT_NEAR((proc_dat_0.predictors.col(0).array()-sorted_data_0.col(2).array()).abs().sum(),0,1e-10)<<"Split_Data predictors(0) post-subsetting does not match data pre-subsetting (2nd iter)";
    ASSERT_EQ(proc_dat_0.predictors_index.at(0),0)<<"Split_Data predictors_index.at(0) has the wrong value post-subsetting (2nd iter)";
    //
    ASSERT_EQ(proc_dat_0.key_to_original_data.at(1),3)<<"Split_Data key_to_original_data.at(1) has the wrong value post-subsetting (2nd iter)";
    ASSERT_NEAR((proc_dat_0.predictors.col(1).array()-sorted_data_0.col(3).array()).abs().sum(),0,1e-10)<<"Split_Data predictors(1) post-subsetting does not match data pre-subsetting (2nd iter)";
    ASSERT_EQ(proc_dat_0.predictors_index.at(1),1)<<"Split_Data predictors_index.at(0) has the wrong value post-subsetting (2nd iter)";
}
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
TEST_F(SplitDataTest, HashCompareTest){
    Eigen::MatrixXd test_input_data_0 = input_data_0;
    deker::fit::Split_Data proc_data_0(input_data_0,which_response,predictors_use_index);
    deker::fit::Split_Data test_proc_data_0(test_input_data_0,which_response,predictors_use_index);
    ASSERT_TRUE(proc_data_0.parent_hash==test_proc_data_0.parent_hash)<<"parent_hash not equal for data_0";
}
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
TEST_F(SplitDataTest, WriteReadTest){
    std::vector<unsigned> subset_predictors_342{3,4,2};
    std::string write_test_loc = "write_test_loc_split_data.bin";
    //
    Eigen::MatrixXd input_data_0_original = input_data_0;
    deker::fit::Split_Data proc_data_0(input_data_0,which_response,predictors_use_index);
    proc_data_0.subset_predictors(subset_predictors_342);
    deker::io::write_to_binary<deker::fit::Split_Data>(write_test_loc,proc_data_0);
    deker::fit::Split_Data proc_data_0_test(input_data_0_original,which_response);
    deker::io::read_from_binary<deker::fit::Split_Data>(write_test_loc,proc_data_0_test);
    //
    remove(write_test_loc.c_str());
    //
    ASSERT_TRUE(proc_data_0.parent_hash==proc_data_0_test.parent_hash)<<"parent_hash not equal for data_0 original vs. write-read";
    ASSERT_EQ(proc_data_0.which_response,proc_data_0_test.which_response)<<"which_response not equal for data_0 original vs write-read";
    ASSERT_EQ(proc_data_0.original_response_sort_index.size(),proc_data_0_test.original_response_sort_index.size())<<"original_response_sort_index.size() not equal for data_0 original vs write-read";
    ASSERT_EQ(proc_data_0.predictors_index.size(),proc_data_0_test.predictors_index.size())<<"predictors_index.size() not equal for data_0 original vs write-read";
    ASSERT_EQ(proc_data_0.key_to_original_data.size(),proc_data_0_test.key_to_original_data.size())<<"key_to_original_data.size() not equal for data_0 original vs write-read";
    //
    bool all_equal = true;
    for(unsigned i = 0; i<proc_data_0.original_response_sort_index.size();i++){
        if(proc_data_0.original_response_sort_index.at(i)!=proc_data_0_test.original_response_sort_index.at(i)){
            all_equal = false;
            break;
        }
    }
    ASSERT_TRUE(all_equal)<<"original_response_sort_index not equal for data_0 original vs write-read";
    //
    all_equal = true;
    for(unsigned i = 0; i<proc_data_0.predictors_index.size();i++){
        if(proc_data_0.predictors_index.at(i)!=proc_data_0_test.predictors_index.at(i)){
            all_equal = false;
            break;
        }
    }
    ASSERT_TRUE(all_equal)<<"predictors_index not equal for data_0 original vs write-read";
    //
    all_equal = true;
    for(unsigned i = 0; i<proc_data_0.key_to_original_data.size();i++){
        if(proc_data_0.key_to_original_data.at(i)!=proc_data_0_test.key_to_original_data.at(i)){
            all_equal = false;
            break;
        }
    }
    ASSERT_TRUE(all_equal)<<"key_to_original_data not equal for data_0 original vs write-read";
    //
    ASSERT_EQ(proc_data_0.response.size(),proc_data_0_test.response.size())<<"response.size() not equal for data_0 original vs. write-read";
    ASSERT_EQ(proc_data_0.predictors.rows(),proc_data_0_test.predictors.rows())<<"predictors.rows() not equal for data_0 original vs. write-read";
    ASSERT_EQ(proc_data_0.predictors.cols(),proc_data_0_test.predictors.cols())<<"predictors.cols() not equal for data_0 original vs. write-read";
    //
    all_equal = true;
    for(unsigned i = 0; i<proc_data_0.response.size();i++){
        if(proc_data_0.response(i)!=proc_data_0_test.response(i)){
            all_equal = false;
            break;
        }
    }
    ASSERT_TRUE(all_equal)<<"response not equal for data_0 original vs write-read";
    //
    all_equal = true;
    for(unsigned i = 0; i<proc_data_0.predictors.rows();i++){
        for(unsigned j = 0; j<proc_data_0.predictors.cols();j++){
            if(proc_data_0.predictors(i,j)!=proc_data_0_test.predictors(i,j)){
                all_equal = false;
                break;
            }   
        }
    }
    ASSERT_TRUE(all_equal)<<"predictors not equal for data_0 original vs write-read";
}