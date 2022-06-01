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
#include <deker/fit/response_labels.h>
#include <deker/io/misc.h>
#include <gtest/gtest.h>

class ResponseLabelsTest : public ::testing::Test{
protected:
    ResponseLabelsTest(){
        known_response.resize(20);
        known_response<<-5,-5,-4,-1,-1,-1,0,0,0,0,1,1,2,2,3,3,4,5,5,9;
    }
    Eigen::VectorXd known_response;
};
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
TEST_F(ResponseLabelsTest, CheckResponseLabelsCalcs){
    deker::fit::Response_Labels response_labels_test(known_response,true);
    Eigen::MatrixXd y_lab_true(20,9);
    y_lab_true << 0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,
                  1,0,0,0,0,0,0,0,0,
                  1,1,0,0,0,0,0,0,0,
                  1,1,0,0,0,0,0,0,0,
                  1,1,0,0,0,0,0,0,0,
                  1,1,1,0,0,0,0,0,0,
                  1,1,1,0,0,0,0,0,0,
                  1,1,1,0,0,0,0,0,0,
                  1,1,1,0,0,0,0,0,0,
                  1,1,1,1,0,0,0,0,0,
                  1,1,1,1,0,0,0,0,0,
                  1,1,1,1,1,0,0,0,0,
                  1,1,1,1,1,0,0,0,0,
                  1,1,1,1,1,1,0,0,0,
                  1,1,1,1,1,1,0,0,0,
                  1,1,1,1,1,1,1,0,0,
                  1,1,1,1,1,1,1,1,0,
                  1,1,1,1,1,1,1,1,0,
                  1,1,1,1,1,1,1,1,1;
    double y_lab_check = (y_lab_true-response_labels_test.get_y_labels(false)).array().abs().sum();
    ASSERT_NEAR(y_lab_check,0,1e-10)<<"y_labels does not match hand-calculated matrix";
    // Eigen::MatrixXd y_pseudo_true(9,20);
    // y_pseudo_true << 0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    //                  0,0,-1,.333,.333,.333,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    //                  0,0,0,-.333,-.333,-.333,.25,.25,.25,.25,0,0,0,0,0,0,0,0,0,0,
    //                  0,0,0,0,0,0,-.25,-.25,-.25,-.25,.5,.5,0,0,0,0,0,0,0,0,
    //                  0,0,0,0,0,0,0,0,0,0,-.5,-.5,.5,.5,0,0,0,0,0,0,
    //                  0,0,0,0,0,0,0,0,0,0,0,0,-.5,-.5,.5,.5,0,0,0,0,
    //                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,-.5,-.5,1,0,0,0,
    //                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,.5,.5,0,
    //                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-.5,-.5,1;
    Eigen::MatrixXd y_pseudo_true(20,20);
    y_pseudo_true << 1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                     -1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                     0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                     0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                     0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                     0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                     0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
                     0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,
                     0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,
                     0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,
                     0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,
                     0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,
                     0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,
                     0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,
                     0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,
                     0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,
                     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,
                     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,
                     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0,
                     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1;
    double y_pseudo_check = (y_pseudo_true-response_labels_test.get_y_labels_pseudoinv()).array().abs().sum();
    ASSERT_NEAR(y_pseudo_check,0,1e-2)<<"y_pseudoinv does not match hand-calculated matrix";
}
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
TEST_F(ResponseLabelsTest, WriteReadTest){
    deker::fit::Response_Labels response_labels_test(known_response,true);
    std::string write_test_loc = "write_test_loc_response_labels.bin";
    deker::io::write_to_binary<deker::fit::Response_Labels>(write_test_loc,response_labels_test);
    deker::fit::Response_Labels response_labels_wr(known_response);
    deker::io::read_from_binary<deker::fit::Response_Labels>(write_test_loc,response_labels_wr);
    remove(write_test_loc.c_str());
    //
    ASSERT_EQ(response_labels_test.get_threshold_inds(true).size(),response_labels_wr.get_threshold_inds(true).size())<<"threshold_inds.size() does not match between original and write-read data";
    ASSERT_EQ(response_labels_test.get_threshold_inds(false).size(),response_labels_wr.get_threshold_inds(false).size())<<"threshold_inds.size() does not match between original and write-read data";
    ASSERT_EQ(response_labels_test.get_threshold_mult(true).size(),response_labels_wr.get_threshold_mult(true).size())<<"threshold_mult.size() does not match between original and write-read data";
    ASSERT_EQ(response_labels_test.get_threshold_mult(false).size(),response_labels_wr.get_threshold_mult(false).size())<<"threshold_mult.size() does not match between original and write-read data";
    ASSERT_EQ(response_labels_test.get_y_labels(true).rows(),response_labels_wr.get_y_labels(true).rows())<<"y_labels.rows() does not match between original and write-read data";
    ASSERT_EQ(response_labels_test.get_y_labels(false).rows(),response_labels_wr.get_y_labels(false).rows())<<"y_labels.rows() does not match between original and write-read data";
    ASSERT_EQ(response_labels_test.get_y_labels(true).cols(),response_labels_wr.get_y_labels(true).cols())<<"y_labels.cols() does not match between original and write-read data";
    ASSERT_EQ(response_labels_test.get_y_labels_pseudoinv().rows(),response_labels_wr.get_y_labels_pseudoinv().rows())<<"y_labels_pseudoinv.rows() does not match between original and write-read data";
    ASSERT_EQ(response_labels_test.get_y_labels_pseudoinv().cols(),response_labels_wr.get_y_labels_pseudoinv().cols())<<"y_labels_pseudoinv.cols() does not match between original and write-read data";
    //
    bool all_equal = true;
    for(unsigned i = 0; i<response_labels_test.get_threshold_inds(true).size();i++){
        if(response_labels_test.get_threshold_inds(true).at(i)!=response_labels_wr.get_threshold_inds(true).at(i)){
            all_equal = false;
            break;
        }
    }
    ASSERT_TRUE(all_equal)<<"threshold_inds does not match between original and write-read data";
    //
    all_equal = true;
    for(unsigned i = 0; i<response_labels_test.get_threshold_inds(false).size();i++){
        if(response_labels_test.get_threshold_inds(false).at(i)!=response_labels_wr.get_threshold_inds(false).at(i)){
            all_equal = false;
            break;
        }
    }
    ASSERT_TRUE(all_equal)<<"threshold_inds does not match between original and write-read data";
    //
    all_equal = true;
    for(unsigned i = 0; i<response_labels_test.get_threshold_mult(true).size();i++){
        if(response_labels_test.get_threshold_mult(true).at(i)!=response_labels_wr.get_threshold_mult(true).at(i)){
            all_equal = false;
            break;
        }
    }
    ASSERT_TRUE(all_equal)<<"threshold_mult does not match between original and write-read data";
    //
    all_equal = true;
    for(unsigned i = 0; i<response_labels_test.get_threshold_mult(false).size();i++){
        if(response_labels_test.get_threshold_mult(false).at(i)!=response_labels_wr.get_threshold_mult(false).at(i)){
            all_equal = false;
            break;
        }
    }
    ASSERT_TRUE(all_equal)<<"threshold_mult does not match between original and write-read data";
    //
    all_equal = true;
    for(unsigned i = 0; i<response_labels_test.get_y_labels(true).rows();i++){
        for(unsigned j = 0; j<response_labels_test.get_y_labels(true).cols();j++){
            if(response_labels_test.get_y_labels(true)(i,j)!=response_labels_wr.get_y_labels(true)(i,j)){
                all_equal = false;
                break;
            }
        }
    }
    ASSERT_TRUE(all_equal)<<"y_labels does not match between original and write-read data";
    //
    all_equal = true;
    for(unsigned i = 0; i<response_labels_test.get_y_labels(false).rows();i++){
        for(unsigned j = 0; j<response_labels_test.get_y_labels(false).cols();j++){
            if(response_labels_test.get_y_labels(false)(i,j)!=response_labels_wr.get_y_labels(false)(i,j)){
                all_equal = false;
                break;
            }
        }
    }
    ASSERT_TRUE(all_equal)<<"y_labels does not match between original and write-read data";
    //
    all_equal = true;
    for(unsigned i = 0; i<response_labels_test.get_y_labels_pseudoinv().rows();i++){
        for(unsigned j = 0; j<response_labels_test.get_y_labels_pseudoinv().cols();j++){
            if(response_labels_test.get_y_labels_pseudoinv()(i,j)!=response_labels_wr.get_y_labels_pseudoinv()(i,j)){
                all_equal = false;
                break;
            }
        }
    }
    ASSERT_TRUE(all_equal)<<"y_labels_pseudoinv does not match between original and write-read data";
    // //
    // all_equal = true;
    // for(unsigned i = 0; i<response_labels_test.y_labels_pseudoident.rows();i++){
    //     for(unsigned j = 0; j<response_labels_test.y_labels_pseudoident.cols();j++){
    //         if(response_labels_test.y_labels_pseudoident(i,j)!=response_labels_wr.y_labels_pseudoident(i,j)){
    //             all_equal = false;
    //             break;
    //         }
    //     }
    // }
    // ASSERT_TRUE(all_equal)<<"y_labels_pseudoident does not match between original and write-read data";
}