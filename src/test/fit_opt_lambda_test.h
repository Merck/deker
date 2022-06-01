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
#include <deker/opt_lambda.h>
#include <gtest/gtest.h>
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
class OptLambdaTest : public ::testing::Test{
protected:
    //dummy class for testing opt_lambda
    class Dummy_Opt{
    public:
        double cutoff_lambda;
        double parabola_opt_lambda;
        Dummy_Opt():cutoff_lambda(.2),parabola_opt_lambda(.15){}
        void Lambda_Obs(const double& lambda_regularization,
                        deker::fit::Output_deker& sol) const
        {
            sol.lambda_regularization = lambda_regularization;
            sol.feature_index.push_back(0);
            if(lambda_regularization>.2){
                sol.return_status = 2;
                sol.BIC = 0;
            }else{
                sol.return_status = 1;
                sol.BIC = pow(lambda_regularization-parabola_opt_lambda,2.0)+.003*sin(lambda_regularization*90.0);
            }
            sol.sigma_kernel_width = 1;
        }
    };
    //
    Dummy_Opt dummy_opt;
    deker::fit::Opt_Param control;
};
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
TEST_F(OptLambdaTest, WriteReadTest){
    std::string write_loc = "write_test_loc_opt_lambda.bin";
    deker::fit::deker_Lambda_Optimizer<Dummy_Opt,&Dummy_Opt::Lambda_Obs> lambda_opt(control);
    bool still_opting = true;
    while(still_opting){
        deker::io::write_to_binary<deker::fit::deker_Lambda_Optimizer<Dummy_Opt,&Dummy_Opt::Lambda_Obs>>(write_loc,lambda_opt);
        deker::fit::deker_Lambda_Optimizer<Dummy_Opt,&Dummy_Opt::Lambda_Obs> lambda_opt_wr;
        deker::io::read_from_binary<deker::fit::deker_Lambda_Optimizer<Dummy_Opt,&Dummy_Opt::Lambda_Obs>>(write_loc,lambda_opt_wr);
        //
        bool rangefind_complete, fill_complete, bo_complete;
        std::tie(rangefind_complete, fill_complete, bo_complete) = lambda_opt.get_status();
        bool rangefind_complete_wr, fill_complete_wr, bo_complete_wr;
        std::tie(rangefind_complete_wr, fill_complete_wr, bo_complete_wr) = lambda_opt_wr.get_status();
        ASSERT_EQ(rangefind_complete,rangefind_complete_wr)<<"rangefind_complete not equal between original and write-read";
        ASSERT_EQ(fill_complete,fill_complete_wr)<<"fill_complete not equal between original and write-read";
        ASSERT_EQ(bo_complete,bo_complete_wr)<<"bo_complete not equal between original and write-read";
        ASSERT_EQ(lambda_opt.get_lambda_max(),lambda_opt_wr.get_lambda_max())<<"lambda_max not equal between original and write-read";
        ASSERT_EQ(lambda_opt.get_sampled_lambda().size(),lambda_opt_wr.get_sampled_lambda().size())<<"sampled_lambda.size() not equal between original and write-read";
        ASSERT_EQ(lambda_opt.get_sampled_BIC().size(),lambda_opt_wr.get_sampled_BIC().size())<<"sampled_lambda.size() not equal between original and write-read";
        //
        bool all_equal = true;
        for(unsigned i = 0; i<lambda_opt.get_sampled_lambda().size();i++){
            if(lambda_opt.get_sampled_lambda().at(i)!=lambda_opt_wr.get_sampled_lambda().at(i)){
                all_equal = false;
                break;
            }
        }
        ASSERT_TRUE(all_equal)<<"sampled_lambda not equal between original and write-read";
        //
        all_equal = true;
        for(unsigned i = 0; i<lambda_opt.get_sampled_BIC().size();i++){
            if(lambda_opt.get_sampled_BIC().at(i)!=lambda_opt_wr.get_sampled_BIC().at(i)){
                all_equal = false;
                break;
            }
        }
        ASSERT_TRUE(all_equal)<<"sampled_BIC not equal between original and write-read";
        //
        still_opting = lambda_opt.opt_iteration(dummy_opt);
        lambda_opt_wr.opt_iteration(dummy_opt);
        //
        if(!fill_complete){
            //rangefind and fill are deterministic, so we can compare solution output to confirm both did the same optimization step
            //bayesian optimization is not deterministic, so skip comparison
            ASSERT_EQ(lambda_opt.get_sampled_lambda().size(),lambda_opt.get_sampled_BIC().size())<<"sampled_lambda.size() not equal to sampled_BIC.size() after opt_iteration";
            ASSERT_EQ(lambda_opt.get_sampled_lambda().size(),lambda_opt_wr.get_sampled_lambda().size())<<"sampled_lambda.size() not equal between original and write-read after opt_iteration";
            ASSERT_EQ(lambda_opt.get_sampled_BIC().size(),lambda_opt_wr.get_sampled_BIC().size())<<"sampled_lambda.size() not equal between original and write-read after opt_iteration";
            if(lambda_opt.get_sampled_lambda().size()>0){
                ASSERT_EQ(lambda_opt.get_sampled_lambda().back(),lambda_opt_wr.get_sampled_lambda().back())<<"sampled_lambda not equal between original and write-read after opt_iteration (pre BO)";
                ASSERT_EQ(lambda_opt.get_sampled_BIC().back(),lambda_opt_wr.get_sampled_BIC().back())<<"sampled_BIC not equal between original and write-read after opt_iteration (pre BO)";
            }
            ASSERT_EQ(lambda_opt.get_lambda_max(),lambda_opt_wr.get_lambda_max())<<"lambda_max not equal between original and write-read after opt_iteration";
        }
    }
    EXPECT_NEAR(lambda_opt.get_lambda_max(),.2,.01)<<"final lambda_max not near expected value";
    EXPECT_NEAR(lambda_opt.get_best_sol().lambda_regularization,.124,.01)<<"minimizing lambda not near expected value";;
    EXPECT_NEAR(lambda_opt.get_best_sol().BIC,-.0023,.0001)<<"minimum BIC not near expected value";
    remove(write_loc.c_str());
}
