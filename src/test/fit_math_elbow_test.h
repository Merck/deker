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
#include <deker/fit/math.h>
#include <gtest/gtest.h>
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
TEST (FitMathElbowTest, FitMathElbowSimple){
    Eigen::VectorXd vals_1(7);
    vals_1 << 0,.1,.15,.2,1,2,5;
    double val_index_1_test = deker::fit::math::find_elbow(vals_1);
    EXPECT_EQ(val_index_1_test,vals_1(4))<<"find_elbow test 1 (simple) does not match expected value";
    ////
    Eigen::VectorXd vals_1shuff(7);
    vals_1shuff << 5,1,2,.2,.1,.15,0;
    double val_index_1shuff_test = deker::fit::math::find_elbow(vals_1shuff);
    EXPECT_EQ(val_index_1shuff_test,vals_1shuff(1))<<"find_elbow test 1 shuffled (simple) does not match expected value";
    ////
    Eigen::VectorXd vals_2(10);
    vals_2 << 0.047425873,0.119202922,0.268941421,0.500000000,0.731058579,
              0.880797078,0.952574127,0.982013790,0.993307149,0.997527377;
    double val_index_2_test = deker::fit::math::find_elbow(vals_2);
    EXPECT_EQ(val_index_2_test,vals_2(1))<<"find_elbow test 2 (sigmoid) does not match expected value";
    ////
    Eigen::VectorXd vals_2shuff(10);
    vals_2shuff << 0.99330715,0.04742587,0.88079708,0.95257413,0.73105858,
                   0.26894142,0.98201379,0.50000000,0.99752738,0.11920292;
    double val_index_2shuff_test = deker::fit::math::find_elbow(vals_2shuff);
    EXPECT_EQ(val_index_2shuff_test,vals_2shuff(9))<<"find_elbow test 2 shuffled (sigmoid) does not match expected value";
    ////
    Eigen::VectorXd vals_3(14);
    vals_3 << 0,0,0,0,0,0,0,0,0,0,4,9,16,25;
    double val_index_3_test = deker::fit::math::find_elbow(vals_3);
    EXPECT_EQ(val_index_3_test,vals_3(9))<<"find_elbow test 3 (0 inflated) does not match expected value";
    ////
    Eigen::VectorXd vals_3shuff(14);
    vals_3shuff << 0,0,9,0,16,0,0,4,0,25,0,0,0,0;
    double val_index_3shuff_test = deker::fit::math::find_elbow(vals_3shuff);
    EXPECT_EQ(val_index_3shuff_test,vals_3shuff(0))<<"find_elbow test 3 shuffled (0 inflated) does not match expected value";
}