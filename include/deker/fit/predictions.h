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
#ifndef DEKER_FIT_PREDICTIONS
#define DEKER_FIT_PREDICTIONS
///////////////////////////////////////////
#include <deker/fit/math.h>
#include <deker/fit/platt_scaling.h>
#include <deker/io/misc.h>
///////////////////////////////////////////
namespace deker{
    namespace fit{
        class Predictions_deker{
        protected:
            const Split_Data& split_data;
            const Response_Labels& response_labels;
            //
            double sigma_kernel_width;
            Eigen::VectorXd v;
            std::vector<unsigned> feature_index;
            //
            double A;
            double B;
            Eigen::VectorXd raw_predicted;
            Eigen::VectorXd platt_probs; 
            //
            double df;
            double logli;
        public:
            /////////////////////////////////////////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////////////////
            Predictions_deker(const Split_Data& split_data,
                              const Response_Labels& response_labels):
            split_data(split_data),
            response_labels(response_labels){}
            /////////////////////////////////////////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////////////////
            Predictions_deker(const Split_Data& split_data,
                              const Response_Labels& response_labels,
                              const double& sigma_kernel_width_in,
                              const Eigen::Ref<const Eigen::VectorXd> v_in,
                              const std::vector<unsigned>& feature_index_in):
            split_data(split_data),
            response_labels(response_labels){
                build(sigma_kernel_width_in,v_in,feature_index_in);
            }
            /////////////////////////////////////////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////////////////
            void build(const double& sigma_kernel_width_in,
                       const Eigen::Ref<const Eigen::VectorXd> v_in,
                       const std::vector<unsigned>& feature_index_in)
            {
                sigma_kernel_width = sigma_kernel_width_in;
                v = v_in;
                feature_index = feature_index_in;
                ///
                Eigen::MatrixXd tmp_z_cache = std::get<0>(math::calc_z_cache(split_data,response_labels,sigma_kernel_width,v,feature_index,false,false));
                //calculate properties = flat, samples in blocks with thresholds in blocks
                //raw_predicted = (tmp_z_cache.array().colwise()*(v.array().square()/v.array().square().sum())).colwise().sum();
                ////////////////////////////////
                raw_predicted.resize(split_data.response.size()+tmp_z_cache.cols());
                raw_predicted << Eigen::VectorXd::Constant(split_data.response.size(),1.0), ((tmp_z_cache.array().colwise()*(v.array().square()/v.array().square().sum())).colwise().sum()).transpose();
                //
                Eigen::Map<Eigen::MatrixXd> pred_mat(raw_predicted.data(),
                                                     response_labels.get_y_labels(true).rows(),
                                                     response_labels.get_y_labels(true).cols()+1);
                
                ////////////////////////////////
                Eigen::ArrayXd raw_predicted_good((response_labels.get_threshold_inds(false).size()-1)*response_labels.get_y_labels(true).rows());
                for(unsigned i=1; i<response_labels.get_threshold_inds(false).size(); i++){
                    raw_predicted_good.segment((i-1)*response_labels.get_y_labels(true).rows(),response_labels.get_y_labels(true).rows()) = pred_mat.col(response_labels.get_threshold_inds(false).at(i));
                }
                //
                Platt_Params platt_params;
                Eigen::Map<const Eigen::ArrayXd> lab_array(response_labels.get_y_labels(false).data(),
                                                           response_labels.get_y_labels(false).size());
                std::tie(A,B) = platt_scaling(raw_predicted_good,lab_array,platt_params);
                ///
                platt_probs = rescale_deci_by_platt(raw_predicted_good,A,B);
                logli = (platt_probs.array()*lab_array+(1-platt_probs.array())*(1-lab_array)).log().sum();
                logli /= (double) response_labels.get_y_labels(false).cols();
                ////////////////////////////////
                for(unsigned i=1; i<pred_mat.cols();i++){
                    pred_mat.col(i) = rescale_deci_by_platt(pred_mat.col(i),A,B);
                }
                Eigen::MatrixXd pseudohat = pred_mat*response_labels.get_y_labels_pseudoinv();
                df = pseudohat.trace();
                if(df < 0){df = std::numeric_limits<double>::max();}
                ////////////////////////////////
                ////do Platt scaling to get probabilities
                // Platt_Params platt_params;
                // Eigen::Map<const Eigen::ArrayXd> lab_array(response_labels.get_y_labels(true).data(),
                //                                            response_labels.get_y_labels(true).size());
                // std::tie(A,B) = platt_scaling(raw_predicted,lab_array,platt_params);
                // ///
                // // platt_probs = rescale_deci_by_platt(raw_predicted,A,B);
                // platt_probs.resize(split_data.response.size()+raw_predicted.size());
                // platt_probs << Eigen::VectorXd::Constant(split_data.response.size(),1.0), rescale_deci_by_platt(raw_predicted,A,B);
                // Eigen::Map<const Eigen::MatrixXd> prob_mat(platt_probs.data(),
                //                                            response_labels.get_y_labels(true).rows(),
                //                                            response_labels.get_y_labels(true).cols()+1);
                ////calculate log likelihood
                ////only use prob_mat columns from true thresholds, not imputed ones
                // logli = 0;
                // for(unsigned i=1; i<response_labels.get_threshold_inds(false).size(); i++){
                //     logli+= (prob_mat.col(response_labels.get_threshold_inds(false).at(i)).array()*response_labels.get_y_labels(false).col(i-1).array()+
                //         (1-prob_mat.col(response_labels.get_threshold_inds(false).at(i)).array())*(1-response_labels.get_y_labels(false).col(i-1).array())).log().sum();
                // }
                // logli /= (double) response_labels.get_y_labels(false).cols();
                ////calculate degrees of freedom
                // Eigen::MatrixXd pseudohat = prob_mat*response_labels.get_y_labels_pseudoinv();
                // df = pseudohat.trace();
                // if(df < 0){df = std::numeric_limits<double>::max();}
            }
            /////////////////////////////////////////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////////////////
            const Eigen::Ref<const Eigen::VectorXd> get_raw_predicted(){return raw_predicted;}
            const Eigen::Ref<const Eigen::VectorXd> get_platt_probs(){return platt_probs;}
            std::tuple<double,double> get_fit(){return std::make_tuple(df,logli);}
            /////////////////////////////////////////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////////////////
            void csv_header(std::ofstream& outfile,
                            const std::vector<std::string>& extra_identifiers)
            {
                for(unsigned i = 0; i<extra_identifiers.size();i++)
                    outfile<<extra_identifiers.at(i)<<",";
                outfile<<"response"<<","<<"threshold"<<","<<"true"<<","<<"raw_margin"<<","<<"platt_prob"<<"\n";
            }
            /////////////////////////////////////////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////////////////
            void csv(std::ofstream& outfile,
                     const std::vector<std::string>& extra_identifiers)
            {
                //prediction data
                //first threshold in response_labels.threshold_inds ignored
                unsigned nrow = response_labels.get_y_labels(true).rows();
                for(unsigned j=1; j<response_labels.get_threshold_inds(false).size(); j++){
                    for(unsigned i=0; i<nrow; i++){
                        for(unsigned k = 0; k<extra_identifiers.size();k++)
                            outfile<<extra_identifiers.at(k)<<",";
                        outfile<<response_labels.get_response()(i)<<",";//<<
                        outfile<<response_labels.get_response()(response_labels.get_threshold_inds(false).at(j))<<",";//<<
                        outfile<<response_labels.get_y_labels(false)(i,j-1)<<",";//<<
                        outfile<<raw_predicted((response_labels.get_threshold_inds(false).at(j)-1)*nrow+i)<<",";//<<
                        outfile<<platt_probs((response_labels.get_threshold_inds(false).at(j)-1)*nrow+i)<<"\n";
                    }
                }
            }
            /////////////////////////////////////////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////////////////
            void write_to_csv(const std::string& file_location,
                              const std::vector<std::string>& extra_identifiers_header,
                              const std::vector<std::string>& extra_identifiers)
            {
                bool write_header = !io::file_exists_test(file_location);
                std::ofstream outfile(file_location,std::ios_base::app);
                if(write_header){
                    csv_header(outfile,extra_identifiers_header);
                }
                csv(outfile,extra_identifiers);
                outfile.close();
            }
            /////////////////////////////////////////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////////////////
            void write_to_csv(const std::string& file_location)
            {
                std::vector<std::string> extra_identifiers_header;
                std::vector<std::string> extra_identifiers;
                write_to_csv(file_location,extra_identifiers_header,extra_identifiers);
            }
        };
        /////
        class Predictions_deker_holdout{
        protected:
            const Split_Data& split_data;
            const Response_Labels& response_labels;
            //
            double sigma_kernel_width;
            Eigen::VectorXd v;
            std::vector<unsigned> feature_index;
            //
            Eigen::VectorXd raw_predicted;
            Eigen::VectorXd platt_probs; 
            //
            const Eigen::Ref<const Eigen::MatrixXd> test_predictors;
            const Eigen::Ref<const Eigen::VectorXd> test_response;
            
        public:
            /////////////////////////////////////////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////////////////
            Predictions_deker_holdout(const Split_Data& split_data,
                                      const Response_Labels& response_labels,
                                      const Eigen::Ref<const Eigen::MatrixXd> test_predictors,
                                      const Eigen::Ref<const Eigen::VectorXd> test_response):
            split_data(split_data),
            response_labels(response_labels),
            test_predictors(test_predictors),
            test_response(test_response){}
            /////////////////////////////////////////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////////////////
            Predictions_deker_holdout(const Split_Data& split_data,
                              const Response_Labels& response_labels,
                              const double& sigma_kernel_width_in,
                              const Eigen::Ref<const Eigen::VectorXd> v_in,
                              const std::vector<unsigned>& feature_index_in,
                              const Eigen::Ref<const Eigen::MatrixXd> test_predictors,
                              const Eigen::Ref<const Eigen::VectorXd> test_response):
            split_data(split_data),
            response_labels(response_labels),
            test_predictors(test_predictors),
            test_response(test_response){
                build(sigma_kernel_width_in,v_in,feature_index_in);
            }
            /////////////////////////////////////////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////////////////
            void build(const double& sigma_kernel_width_in,
                       const Eigen::Ref<const Eigen::VectorXd> v_in,
                       const std::vector<unsigned>& feature_index_in)
            {
                sigma_kernel_width = sigma_kernel_width_in;
                v = v_in;
                feature_index = feature_index_in;
                //
                Eigen::MatrixXd tmp_z_cache = math::calc_z_cache_holdout(split_data,response_labels,sigma_kernel_width,v,feature_index,test_predictors);
                //calculate properties = flat, samples in blocks with thresholds in blocks
                raw_predicted = (tmp_z_cache.array().colwise()*(v.array().square()/v.array().square().sum())).colwise().sum();
            }
            /////////////////////////////////////////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////////////////
            const Eigen::Ref<const Eigen::VectorXd> get_raw_predicted(){return raw_predicted;}
            /////////////////////////////////////////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////////////////
            void csv_header(std::ofstream& outfile,
                            const std::vector<std::string>& extra_identifiers)
            {
                for(unsigned i = 0; i<extra_identifiers.size();i++)
                    outfile<<extra_identifiers.at(i)<<",";
                outfile<<"response"<<","<<"threshold"<<","<<"true"<<","<<"raw_margin"<<"\n";
            }
            /////////////////////////////////////////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////////////////
            void csv(std::ofstream& outfile,
                     const std::vector<std::string>& extra_identifiers)
            {
                //prediction data
                //first threshold in response_labels.threshold_inds ignored
                //unsigned nrow = response_labels.get_y_labels(true).rows();
                unsigned nrow = test_response.size();
                for(unsigned j=1; j<response_labels.get_threshold_inds(false).size(); j++){
                    for(unsigned i=0; i<test_response.size(); i++){
                        for(unsigned k = 0; k<extra_identifiers.size();k++){
                            outfile<<extra_identifiers.at(k)<<",";
                        }
                        outfile<<test_response(i)<<","<<
                            response_labels.get_response()(response_labels.get_threshold_inds(false).at(j))<<","<<
                                (test_response(i)>=response_labels.get_response()(response_labels.get_threshold_inds(false).at(j)))<<","<<
                                    raw_predicted((j-1)*nrow+i)<<"\n";
                    }
                }
            }
            /////////////////////////////////////////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////////////////
            void write_to_csv(const std::string& file_location,
                              const std::vector<std::string>& extra_identifiers_header,
                              const std::vector<std::string>& extra_identifiers)
            {
                bool write_header = !io::file_exists_test(file_location);
                std::ofstream outfile(file_location,std::ios_base::app);
                if(write_header){
                    csv_header(outfile,extra_identifiers_header);
                }
                csv(outfile,extra_identifiers);
                outfile.close();
            }
            /////////////////////////////////////////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////////////////
            void write_to_csv(const std::string& file_location)
            {
                std::vector<std::string> extra_identifiers_header;
                std::vector<std::string> extra_identifiers;
                write_to_csv(file_location,extra_identifiers_header,extra_identifiers);
            }
        };
    }
}
#endif