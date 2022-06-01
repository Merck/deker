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
#include <deker/io/digest.h>
#include <gtest/gtest.h>
#include <experimental/filesystem>
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
TEST(DigestTest, CheckResponseLabelsCalcs){
    std::string digest_write_path = "digest_test/digest_test_loc";
    std::experimental::filesystem::path p(digest_write_path);
    std::string digest_dir = p.parent_path().string();
    std::experimental::filesystem::create_directory(digest_dir);
    //
    std::vector<unsigned> use_response{0,1,3,5,6,4,8,11,19};
    deker::io::Control_Base control;
    deker::io::Control_Digest digest(digest_dir,control);
    for(unsigned i = 0; i<use_response.size();i++){
        digest.add_record(use_response.at(i));
    }
    deker::io::write_to_binary<deker::io::Control_Digest>(digest_write_path,digest);
    deker::io::Control_Digest digest_wr;
    deker::io::read_from_binary<deker::io::Control_Digest>(digest_write_path,digest_wr);
    //
    ASSERT_EQ(digest.get_digest_dir(),digest_wr.get_digest_dir())<<"get_digest_dir() is not equal after write-read";
    ASSERT_EQ(digest.get_control().data_file,digest_wr.get_control().data_file)<<"get_control().data_file is not equal after write-read";
    ASSERT_EQ(digest.get_control().lambda_iter_per_job,digest_wr.get_control().lambda_iter_per_job)<<"get_control().lambda_iter_per_job is not equal after write-read";
    ASSERT_EQ(digest.get_control().opt_param.lambda_init_max,digest_wr.get_control().opt_param.lambda_init_max)<<"get_control().opt_param.lambda_init_max is not equal after write-read";
    ASSERT_EQ(digest.get_control().opt_param.lambda_init_min,digest_wr.get_control().opt_param.lambda_init_min)<<"get_control().opt_param.lambda_init_max is not equal after write-read";
    ASSERT_EQ(digest.get_control().opt_param.lambda_init_stepsize,digest_wr.get_control().opt_param.lambda_init_stepsize)<<"get_control().opt_param.lambda_init_stepsize is not equal after write-read";
    ASSERT_EQ(digest.get_control().opt_param.lambda_init_maxiter,digest_wr.get_control().opt_param.lambda_init_maxiter)<<"get_control().opt_param.lambda_init_maxiter is not equal after write-read";
    ASSERT_EQ(digest.get_control().opt_param.lambda_fill_mindist,digest_wr.get_control().opt_param.lambda_fill_mindist)<<"get_control().opt_param.lambda_fill_mindist is not equal after write-read";
    ASSERT_EQ(digest.get_control().opt_param.lambda_fill_maxiter,digest_wr.get_control().opt_param.lambda_fill_maxiter)<<"get_control().opt_param.lambda_fill_maxiter is not equal after write-read";
    ASSERT_EQ(digest.get_control().opt_param.lambda_bo_maxiter,digest_wr.get_control().opt_param.lambda_bo_maxiter)<<"get_control().opt_param.lambda_bo_maxiter is not equal after write-read";
    ASSERT_EQ(digest.get_control().opt_param.h_huber_loss,digest_wr.get_control().opt_param.h_huber_loss)<<"get_control().opt_param.h_huber_loss is not equal after write-read";
    ASSERT_EQ(digest.get_control().opt_param.hinge_max,digest_wr.get_control().opt_param.hinge_max)<<"get_control().opt_param.hinge_max is not equal after write-read";
    ASSERT_EQ(digest.get_control().opt_param.w_convergence_threshold,digest_wr.get_control().opt_param.w_convergence_threshold)<<"get_control().opt_param.w_convergence_threshold is not equal after write-read";
    ASSERT_EQ(digest.get_control().opt_param.w_maxiter,digest_wr.get_control().opt_param.w_maxiter)<<"get_control().opt_param.w_maxiter is not equal after write-read";
    ASSERT_EQ(digest.get_control().opt_param.feature_drop_threshold,digest_wr.get_control().opt_param.feature_drop_threshold)<<"get_control().opt_param.feature_drop_threshold is not equal after write-read";
    //
    ASSERT_EQ(digest.get_record_size(),digest_wr.get_record_size())<<"get_record_size() is not equal after write-read";
    //
    for(size_t i = 0; i<digest.get_record_size();i++){
        ASSERT_EQ(digest_wr.get_digest_record(i).which_response,digest.get_digest_record(i).which_response)<<"which response not equal for digest & digest_wr after write-read";
        //
        std::string sol_file_loc = std::string(digest.get_digest_record(i).solve_file);
        std::string sol_file_loc_test = std::string(digest_wr.get_digest_record(i).solve_file);
        ASSERT_EQ(sol_file_loc,sol_file_loc_test)<<"solve_file not equal after write-read";
        //
        std::string lopt_file_loc = digest.get_digest_dir()+std::string(digest.get_digest_record(i).opt_file);
        std::string lopt_file_loc_test = digest.get_digest_dir()+std::string(digest_wr.get_digest_record(i).opt_file);
        ASSERT_EQ(lopt_file_loc,lopt_file_loc_test)<<"opt_file not equal after write-read";
        //
        std::string run_file_loc = digest.get_digest_dir()+std::string(digest.get_digest_record(i).run_file);
        std::string run_file_loc_test = digest.get_digest_dir()+std::string(digest_wr.get_digest_record(i).run_file);
        ASSERT_EQ(run_file_loc,run_file_loc_test)<<"run_file not equal after write-read";
        //
        std::string out_file_loc = digest.get_digest_dir()+std::string(digest.get_digest_record(i).output_file);
        std::string out_file_loc_test = digest.get_digest_dir()+std::string(digest_wr.get_digest_record(i).output_file);
        ASSERT_EQ(out_file_loc,out_file_loc_test)<<"output_file not equal after write-read";
    }
    deker::io::Control_Digest single_digest;
    for(unsigned i = 0; i<use_response.size();i++){
        single_digest.read_single_record(digest_write_path,i);
        //
        ASSERT_EQ(digest.get_digest_dir(),single_digest.get_digest_dir())<<"get_digest_dir() is not equal after write-read";
        ASSERT_EQ(digest.get_control().data_file,single_digest.get_control().data_file)<<"get_control().data_file is not equal after write-read";
        ASSERT_EQ(digest.get_control().lambda_iter_per_job,single_digest.get_control().lambda_iter_per_job)<<"get_control().lambda_iter_per_job is not equal after write-read";
        ASSERT_EQ(digest.get_control().opt_param.lambda_init_max,single_digest.get_control().opt_param.lambda_init_max)<<"get_control().opt_param.lambda_init_max is not equal after write-read";
        ASSERT_EQ(digest.get_control().opt_param.lambda_init_min,single_digest.get_control().opt_param.lambda_init_min)<<"get_control().opt_param.lambda_init_max is not equal after write-read";
        ASSERT_EQ(digest.get_control().opt_param.lambda_init_stepsize,single_digest.get_control().opt_param.lambda_init_stepsize)<<"get_control().opt_param.lambda_init_stepsize is not equal after write-read";
        ASSERT_EQ(digest.get_control().opt_param.lambda_init_maxiter,single_digest.get_control().opt_param.lambda_init_maxiter)<<"get_control().opt_param.lambda_init_maxiter is not equal after write-read";
        ASSERT_EQ(digest.get_control().opt_param.lambda_fill_mindist,single_digest.get_control().opt_param.lambda_fill_mindist)<<"get_control().opt_param.lambda_fill_mindist is not equal after write-read";
        ASSERT_EQ(digest.get_control().opt_param.lambda_fill_maxiter,single_digest.get_control().opt_param.lambda_fill_maxiter)<<"get_control().opt_param.lambda_fill_maxiter is not equal after write-read";
        ASSERT_EQ(digest.get_control().opt_param.lambda_bo_maxiter,single_digest.get_control().opt_param.lambda_bo_maxiter)<<"get_control().opt_param.lambda_bo_maxiter is not equal after write-read";
        ASSERT_EQ(digest.get_control().opt_param.h_huber_loss,single_digest.get_control().opt_param.h_huber_loss)<<"get_control().opt_param.h_huber_loss is not equal after write-read";
        ASSERT_EQ(digest.get_control().opt_param.hinge_max,single_digest.get_control().opt_param.hinge_max)<<"get_control().opt_param.hinge_max is not equal after write-read";
        ASSERT_EQ(digest.get_control().opt_param.w_convergence_threshold,single_digest.get_control().opt_param.w_convergence_threshold)<<"get_control().opt_param.w_convergence_threshold is not equal after write-read";
        ASSERT_EQ(digest.get_control().opt_param.w_maxiter,single_digest.get_control().opt_param.w_maxiter)<<"get_control().opt_param.w_maxiter is not equal after write-read";
        ASSERT_EQ(digest.get_control().opt_param.feature_drop_threshold,single_digest.get_control().opt_param.feature_drop_threshold)<<"get_control().opt_param.feature_drop_threshold is not equal after write-read";
        ASSERT_EQ(digest.get_digest_record(i).which_response,single_digest.get_digest_record(0).which_response)<<"which response not equal for digest & single_digest after write-read";
        //
        std::string sol_file_loc = digest.get_digest_dir()+std::string(digest.get_digest_record(i).solve_file);
        std::string sol_file_loc_test = digest.get_digest_dir()+std::string(single_digest.get_digest_record(0).solve_file);
        ASSERT_EQ(sol_file_loc,sol_file_loc_test)<<"solve_file not equal for digest & single_digest after write-read";
        //
        std::string lopt_file_loc = digest.get_digest_dir()+std::string(digest.get_digest_record(i).opt_file);
        std::string lopt_file_loc_test = digest.get_digest_dir()+std::string(single_digest.get_digest_record(0).opt_file);
        ASSERT_EQ(lopt_file_loc,lopt_file_loc_test)<<"opt_file not equal for digest & single_digest after write-read";
        //
        std::string run_file_loc = digest.get_digest_dir()+std::string(digest.get_digest_record(i).run_file);
        std::string run_file_loc_test = digest.get_digest_dir()+std::string(single_digest.get_digest_record(0).run_file);
        ASSERT_EQ(run_file_loc,run_file_loc_test)<<"run_file not equal for digest & single_digest after write-read";
        //
        std::string out_file_loc = digest.get_digest_dir()+std::string(digest.get_digest_record(i).output_file);
        std::string out_file_loc_test = digest.get_digest_dir()+std::string(single_digest.get_digest_record(0).output_file);
        ASSERT_EQ(out_file_loc,out_file_loc_test)<<"out_file not equal for digest & single_digest after write-read";
    }
    std::experimental::filesystem::remove_all(digest_dir);
}
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////