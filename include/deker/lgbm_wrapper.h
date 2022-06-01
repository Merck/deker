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
#ifndef DEKER_FIT_LGBM
#define DEKER_FIT_LGBM
///////////////////////////////////////////
#include <LightGBM/c_api.h>
#include <Eigen/Core>
#include <vector>
#include <fstream>
///////////////////////////////////////////
namespace deker{
    struct LGBM_out{
        int booster_finished;
        int iter;
        double fit_metric;
        Eigen::VectorXd feature_imp_split;
        Eigen::VectorXd feature_imp_gain;
        /////////////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////
        LGBM_out(){
            booster_finished = 0;
            iter = 0;
        }
        /////////////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////
        void serialize(std::ifstream& infile)
        {
            infile.read((char*) (&booster_finished),sizeof(booster_finished));
            infile.read((char*) (&iter),sizeof(iter));
            infile.read((char*) (&fit_metric),sizeof(fit_metric));
            //
            Eigen::MatrixXd::Index vec_size;
            infile.read((char*) (&vec_size),sizeof(Eigen::MatrixXd::Index));
            feature_imp_split.resize(vec_size);
            infile.read((char*) feature_imp_split.data(),vec_size*sizeof(Eigen::MatrixXd::Scalar));
            infile.read((char*) (&vec_size),sizeof(Eigen::MatrixXd::Index));
            feature_imp_gain.resize(vec_size);
            infile.read((char*) feature_imp_gain.data(),vec_size*sizeof(Eigen::MatrixXd::Scalar));
        }
        /////////////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////
        void serialize(std::ofstream& outfile) const
        {
            outfile.write((char*) (&booster_finished),sizeof(booster_finished));
            outfile.write((char*) (&iter),sizeof(iter));
            outfile.write((char*) (&fit_metric),sizeof(fit_metric));
            //write out feature_imp_split
            Eigen::MatrixXd::Index vec_size=feature_imp_split.size();
            outfile.write((char*) (&vec_size), sizeof(Eigen::MatrixXd::Index));
            outfile.write((char*) feature_imp_split.data(), vec_size*sizeof(Eigen::MatrixXd::Scalar));
            //write out feature_gain_split
            vec_size=feature_imp_gain.size();
            outfile.write((char*) (&vec_size), sizeof(Eigen::MatrixXd::Index));
            outfile.write((char*) feature_imp_gain.data(), vec_size*sizeof(Eigen::MatrixXd::Scalar));
        }
    };
    
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    class LGBM_wrapper{
    protected:
        const Eigen::Ref<const Eigen::MatrixXd> predictors;
        const Eigen::Ref<const Eigen::VectorXd> response;
        char *parameters;
        //
        DatasetHandle lgbm_data;
        BoosterHandle lgbm_booster;
        //
        int booster_check;
        // int booster_finished;
        // int iter;
        // double fit_metric;
        // Eigen::VectorXd feature_imp_split;
        // Eigen::VectorXd feature_imp_gain;
        LGBM_out output;
    public:
        LGBM_wrapper(const Eigen::Ref<const Eigen::MatrixXd> predictors,
                     const Eigen::Ref<const Eigen::VectorXd> response,
                     char *parameters):
        predictors(predictors),response(response),parameters(parameters){
            int dataset_create = LGBM_DatasetCreateFromMat(predictors.data(),C_API_DTYPE_FLOAT64,predictors.rows(),predictors.cols(),0,parameters,nullptr,&lgbm_data);
            char *label_field = "label";
            Eigen::VectorXf response_f = response.cast<float>();
            dataset_create = LGBM_DatasetSetField(lgbm_data, label_field, response_f.data(), response_f.size(), C_API_DTYPE_FLOAT32);
            if(dataset_create==-1){
                throw std::runtime_error("LightGBM failed to create the dataset");
            }
            booster_check = -1;
            output.booster_finished = 0;
            output.iter = 0;
        }
        /////////////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////
        void build_booster(){
            int booster_create = LGBM_BoosterCreate(lgbm_data, parameters, &lgbm_booster);
            if(booster_create==-1){
                throw std::runtime_error("LightGBM failed to create the boosting model");
            }
        }
        /////////////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////
        void run(unsigned maxiter)
        {
            for(output.iter = 0; output.iter<maxiter;output.iter++){
                booster_check = LGBM_BoosterUpdateOneIter(lgbm_booster, &output.booster_finished);
                if(output.booster_finished==1||booster_check==-1){
                    break;
                }
            }
            int n_eval;
            LGBM_BoosterGetEvalCounts(lgbm_booster, &n_eval);
            size_t buffer_len = 10;
            size_t buffer_len_ind[n_eval];
            //
            char** out_names;
            out_names = (char**)malloc(sizeof(char*)*n_eval+1);
            for(unsigned i = 0;i < n_eval ; i++)
                out_names[i] = (char*)malloc(sizeof(char)*buffer_len);
            //
            LGBM_BoosterGetEvalNames(lgbm_booster, n_eval, &n_eval, buffer_len, buffer_len_ind, out_names);
            double eval_out[n_eval];
            LGBM_BoosterGetEval(lgbm_booster, 0, &n_eval, eval_out);
            for(unsigned i=0;i<n_eval; i++) 
                free(out_names[i]);
            free(out_names);
            //
            output.fit_metric = eval_out[0];
            /////
            int n_features;
            LGBM_BoosterGetNumFeature(lgbm_booster, &n_features);
            double importance[n_features];
            /////
            int feature_imp = LGBM_BoosterFeatureImportance(lgbm_booster, 0, C_API_FEATURE_IMPORTANCE_SPLIT, importance);
            if(feature_imp==-1){
                throw std::runtime_error("LightGBM failed to calculate importance");
            }
            output.feature_imp_split = Eigen::Map<Eigen::VectorXd>(importance,n_features);
            /////
            feature_imp = LGBM_BoosterFeatureImportance(lgbm_booster, 0, C_API_FEATURE_IMPORTANCE_GAIN, importance);
            if(feature_imp==-1){
                throw std::runtime_error("LightGBM failed to calculate importance");
            }
            output.feature_imp_gain = Eigen::Map<Eigen::VectorXd>(importance,n_features);
        }
        /////////////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////
        const LGBM_out& get_output() const {return output;}
        // const Eigen::VectorXd& get_feature_imp_split() const {return feature_imp_split;}
        // const Eigen::VectorXd& get_feature_imp_gain() const {return feature_imp_gain;}
        // const double& get_fit_metric() const {return fit_metric;}
        // const int& get_iter() const {return iter;}
        // //
        // void set_feature_imp_split(const Eigen::Ref<const Eigen::VectorXd> in_imp_split) {feature_imp_split = in_imp_split;}
        // void set_feature_imp_gain(const Eigen::Ref<const Eigen::VectorXd> in_imp_gain) {feature_imp_gain = in_imp_gain;}
        // void set_fit_metric(const double& in_fit) {fit_metric = in_fit;}
        /////////////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////
        void serialize(std::ifstream& infile)
        {
            //add parent hash check
            infile.read((char*) (&booster_check),sizeof(booster_check));
            output.serialize(infile);
            // infile.read((char*) (&booster_finished),sizeof(booster_finished));
            // infile.read((char*) (&iter),sizeof(iter));
            // infile.read((char*) (&fit_metric),sizeof(fit_metric));
            // //
            // Eigen::MatrixXd::Index vec_size;
            // infile.read((char*) (&vec_size),sizeof(Eigen::MatrixXd::Index));
            // feature_imp_split.resize(vec_size);
            // infile.read((char*) feature_imp_split.data(),vec_size*sizeof(Eigen::MatrixXd::Scalar));
            // infile.read((char*) (&vec_size),sizeof(Eigen::MatrixXd::Index));
            // feature_imp_gain.resize(vec_size);
            // infile.read((char*) feature_imp_gain.data(),vec_size*sizeof(Eigen::MatrixXd::Scalar));
        }
        /////////////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////
        void serialize(std::ofstream& outfile) const
        {
            //add parent hash check
            //figure out how to load model
            outfile.write((char*) (&booster_check),sizeof(booster_check));
            output.serialize(outfile);
            // outfile.write((char*) (&booster_finished),sizeof(booster_finished));
            // outfile.write((char*) (&iter),sizeof(iter));
            // outfile.write((char*) (&fit_metric),sizeof(fit_metric));
            // //write out feature_imp_split
            // Eigen::MatrixXd::Index vec_size=feature_imp_split.size();
            // outfile.write((char*) (&vec_size), sizeof(Eigen::MatrixXd::Index));
            // outfile.write((char*) feature_imp_split.data(), vec_size*sizeof(Eigen::MatrixXd::Scalar));
            // //write out feature_gain_split
            // vec_size=feature_imp_gain.size();
            // outfile.write((char*) (&vec_size), sizeof(Eigen::MatrixXd::Index));
            // outfile.write((char*) feature_imp_gain.data(), vec_size*sizeof(Eigen::MatrixXd::Scalar));
        }
    };
}
#endif