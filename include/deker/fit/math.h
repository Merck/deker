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
#ifndef DEKER_FIT_MATH
#define DEKER_FIT_MATH
///////////////////////////////////////////
#include <deker/fit/split_data.h>
#include <deker/fit/response_labels.h>
#include <deker/fit/opt_param.h>
//
#include <limits>
///////////////////////////////////////////
namespace deker{
  namespace fit{
    namespace math{
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      Eigen::MatrixXd calc_kernel_matrix(const Split_Data& input_data,
                                         const double& sigma_kernel_width,
                                         const Eigen::Ref<const Eigen::VectorXd> v,
                                         const std::vector<unsigned>& feature_index)
      {
        double kernel_multiplier = 1/(2*pow(exp(sigma_kernel_width),2));
        //constuct kernel matrix
        Eigen::MatrixXd kernel_matrix = Eigen::MatrixXd::Zero(input_data.response.size(),input_data.response.size());
        //entries equal kernel(i,j)
        for(unsigned i = 0; i<kernel_matrix.rows();i++){
          kernel_matrix(i,i) = 0;
          for(unsigned j = i+1; j<kernel_matrix.cols();j++){
            for(unsigned k = 0; k<feature_index.size(); k++){
              //distance calculation for kernel function
              kernel_matrix(i,j) += pow((input_data.predictors(i,feature_index.at(k))-input_data.predictors(j,feature_index.at(k)))*v(k)*v(k),2);
            }
            kernel_matrix(i,j) = -kernel_matrix(i,j)*kernel_multiplier;
            kernel_matrix(j,i) = kernel_matrix(i,j);
          } 
        }
        return kernel_matrix;
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      Eigen::MatrixXd calc_weighting_piece(const unsigned n,
                                           const Response_Labels& response_labels,
                                           const bool& use_lab_fix,
                                           const Eigen::Ref<const Eigen::MatrixXd> kernel_matrix,
                                           const bool& label_type)
      {
        Eigen::MatrixXd wt_mat;
        if(label_type){
          //take only kernel_matrix where y_labels == 1
          wt_mat = kernel_matrix.bottomRows(response_labels.get_y_labels(use_lab_fix).rows()-response_labels.get_threshold_inds(use_lab_fix).at(n+1));
        }else{
          //take only kernel_matrix where y_labels == 0
          wt_mat = kernel_matrix.topRows(response_labels.get_threshold_inds(use_lab_fix).at(n+1));
        }
        wt_mat.rowwise() -= wt_mat.colwise().maxCoeff();
        wt_mat = wt_mat.array().exp();
        for(unsigned i = 0; i<wt_mat.cols();i++){
          double col_sum = wt_mat.col(i).sum();
          if(col_sum>0){
            wt_mat.col(i).array() /= col_sum;
          }
        }
        if(label_type){
          wt_mat *= -1;
        }
        return wt_mat;
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      Eigen::MatrixXd calc_weighting_matrix(const unsigned n, 
                                            const Response_Labels& response_labels,
                                            const bool& use_lab_fix,
                                            const Eigen::Ref<const Eigen::MatrixXd> kernel_matrix)
      {
        Eigen::MatrixXd split_false = calc_weighting_piece(n,response_labels,use_lab_fix,kernel_matrix,false);
        Eigen::MatrixXd split_true = calc_weighting_piece(n,response_labels,use_lab_fix,kernel_matrix,true);
        Eigen::MatrixXd to_return(split_false.rows()+split_true.rows(),split_false.cols());
        to_return << split_false,split_true;
        return to_return;
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      std::tuple<Eigen::MatrixXd,Eigen::VectorXd> calc_z_cache(const Split_Data& input_data,
                                                               const Response_Labels& response_labels,
                                                               const double& sigma_kernel_width,
                                                               const Eigen::Ref<const Eigen::VectorXd> v,
                                                               const std::vector<unsigned>& feature_index,
                                                               const bool& exclude_diag = true,
                                                               const bool& opt_w_prep = false)
      {
        Eigen::MatrixXd kernel_matrix = calc_kernel_matrix(input_data,
                                                           sigma_kernel_width,
                                                           v,
                                                           feature_index);
        //
        if(exclude_diag){
          //set diagonals to (sortof) infinite distance to remove them from calculation
          kernel_matrix.diagonal() = Eigen::VectorXd::Constant(kernel_matrix.rows(),std::numeric_limits<double>::lowest());
        }
        //
        std::vector<Eigen::MatrixXd> threshold_diff_mat;
        for(unsigned n = 0; n<response_labels.get_threshold_inds(!opt_w_prep).size()-1; n++){
          threshold_diff_mat.push_back(calc_weighting_matrix(n,response_labels,!opt_w_prep,kernel_matrix));
        }
        //
        Eigen::MatrixXd tmp_z_cache = Eigen::MatrixXd::Zero(feature_index.size(),response_labels.get_y_labels(!opt_w_prep).size());
        Eigen::VectorXd tmp_z_cache_col_multiple;
        tmp_z_cache_col_multiple = Eigen::VectorXd::Constant(tmp_z_cache.cols(),1);
        //
        for(unsigned feature_i = 0; feature_i<feature_index.size();feature_i++){
          unsigned z_cache_col_index = 0; //the next column to start filling
          Eigen::MatrixXd feature_diff_mat = input_data.predictors.col(feature_index.at(feature_i)).replicate(1,input_data.predictors.rows());
          feature_diff_mat.rowwise() -= input_data.predictors.col(feature_index.at(feature_i)).transpose();
          feature_diff_mat = feature_diff_mat.array().abs();
          for(unsigned n = 0; n<response_labels.get_threshold_inds(!opt_w_prep).size()-1; n++){
            Eigen::MatrixXd all_diff_mat = feature_diff_mat.array()*threshold_diff_mat.at(n).array();
            tmp_z_cache.block(feature_i,z_cache_col_index,1,threshold_diff_mat.at(n).cols()) = all_diff_mat.colwise().sum();
            z_cache_col_index += threshold_diff_mat.at(n).cols();
          }
        }
        if(opt_w_prep){
          //multiply z_cache columns by their label as -1/1 (so positive is correct classification, and negative is missclassification)
          //once this is done, ordering of columns does not need to be preserved
          Eigen::Map<const Eigen::ArrayXd> lab_array(response_labels.get_y_labels(!opt_w_prep).data(),
                                                     response_labels.get_y_labels(!opt_w_prep).size());
          tmp_z_cache.array().rowwise() *= (lab_array.transpose()*2-1);
          //walk backwards through columns and consolidate columns with identical multiples
          unsigned last_column = tmp_z_cache.cols()-1;
          for(int i = tmp_z_cache.cols()-1; i>=0; i--){
            int which_sample = i%input_data.response.size();
            int which_threshold = 1+i/((int)input_data.response.size());
            //check if the column sample is the value the column y_label was thresholded on
            //if so, remove it - for w opt we're doing LOOCV
            //also drop if the column is all 0's
            if(((which_sample>=response_labels.get_threshold_inds(!opt_w_prep).at(which_threshold))&&
               (which_sample<(response_labels.get_threshold_inds(!opt_w_prep).at(which_threshold)+response_labels.get_threshold_mult(!opt_w_prep).at(which_threshold)))&&
               (response_labels.get_threshold_mult(!opt_w_prep).at(which_threshold)==1))||(
               tmp_z_cache.col(i).minCoeff()==0&&tmp_z_cache.col(i).maxCoeff()==0)){
              //drops column by replacing it with the last good column
              //will be removed with conservativeResize once all columns have been sorted to the end
              tmp_z_cache.col(i) = tmp_z_cache.col(last_column);
              tmp_z_cache_col_multiple(i) = tmp_z_cache_col_multiple(last_column);
              last_column--;
            }else if(i>=input_data.response.size()){
              //check of column is unique to next column with the same sample i
              //remove and add mult if not
              double maximum_difference_to_next_threshold = (tmp_z_cache.col(i)-tmp_z_cache.col(i-input_data.response.size())).maxCoeff();
              if(maximum_difference_to_next_threshold<1e-3){
                tmp_z_cache_col_multiple(i-input_data.response.size())+=tmp_z_cache_col_multiple(i);
                tmp_z_cache.col(i) = tmp_z_cache.col(last_column);
                tmp_z_cache_col_multiple(i) = tmp_z_cache_col_multiple(last_column);
                last_column--;
              }
            }
          }
          tmp_z_cache.conservativeResize(tmp_z_cache.rows(),last_column+1);
          tmp_z_cache_col_multiple.conservativeResize(last_column+1);
        }
        return std::make_tuple(tmp_z_cache,tmp_z_cache_col_multiple);
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      Eigen::MatrixXd calc_z_cache_holdout(const Split_Data& input_data,
                                           const Response_Labels& response_labels,
                                           const double& sigma_kernel_width,
                                           const Eigen::Ref<const Eigen::VectorXd> v,
                                           const std::vector<unsigned>& feature_index,
                                           const Eigen::Ref<const Eigen::MatrixXd> test_predictors)
      {
        double kernel_multiplier = 1/(2*pow(exp(sigma_kernel_width),2));
        //constuct kernel matrix
        Eigen::MatrixXd kernel_matrix = Eigen::MatrixXd::Zero(input_data.response.size(),test_predictors.rows());
        //entries equal kernel(i,j)
        for(unsigned i = 0; i<kernel_matrix.rows();i++){
          for(unsigned j = 0; j<test_predictors.rows();j++){
            for(unsigned k = 0; k<feature_index.size(); k++){
              //distance calculation for kernel function
              kernel_matrix(i,j) += pow((input_data.predictors(i,feature_index.at(k))-test_predictors(j,feature_index.at(k)))*v(k)*v(k),2);
            }
            kernel_matrix(i,j) = -kernel_matrix(i,j)*kernel_multiplier;
          } 
        }
        std::vector<Eigen::MatrixXd> threshold_diff_mat;
        for(unsigned n = 0; n<response_labels.get_threshold_inds(false).size()-1; n++){
          threshold_diff_mat.push_back(calc_weighting_matrix(n,response_labels,false,kernel_matrix));
        }
        Eigen::MatrixXd tmp_z_cache = Eigen::MatrixXd::Zero(feature_index.size(),threshold_diff_mat.size()*test_predictors.rows());
        for(unsigned feature_i = 0; feature_i<feature_index.size();feature_i++){
          unsigned z_cache_col_index = 0; //the next column to start filling
          Eigen::MatrixXd feature_diff_mat = input_data.predictors.col(feature_index.at(feature_i)).replicate(1,test_predictors.rows());
          feature_diff_mat.rowwise() -= test_predictors.col(feature_index.at(feature_i)).transpose();
          feature_diff_mat = feature_diff_mat.array().abs();
          for(unsigned n = 0; n<threshold_diff_mat.size(); n++){
            Eigen::MatrixXd all_diff_mat = feature_diff_mat.array()*threshold_diff_mat.at(n).array();
            tmp_z_cache.block(feature_i,z_cache_col_index,1,threshold_diff_mat.at(n).cols()) = all_diff_mat.colwise().sum();
            z_cache_col_index += threshold_diff_mat.at(n).cols();
          } 
        }
        return tmp_z_cache;
      } 
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      //z_cache_ptr.calculate_v_gradient(v,control.h_huber_los,control.hinge_max,lambda_regularization,fx,gradient);
      void calculate_v_gradient(const Opt_Param& control,
                                const Eigen::Ref<const Eigen::VectorXd> v,
                                const Eigen::Ref<const Eigen::MatrixXd> z_cache,
                                const Eigen::Ref<const Eigen::VectorXd> z_cache_col_multiple,
                                const double& lambda_regularization,
                                double& fx, 
                                Eigen::VectorXd& gradient)
      {
        fx = 0;
        gradient = Eigen::VectorXd::Zero(v.size());
        //
        Eigen::VectorXd v_sq = v.array().square();
        //first, #columns * (#rows multiplications)
        //then, #columns * (#rows sum)
        Eigen::VectorXd margin = (z_cache.array().colwise()*v_sq.array()).colwise().sum();
        //#columns operations
        for(unsigned i = 0; i<margin.size();i++){
          if(margin(i)<=control.hinge_max+control.h_huber_loss){
            if(margin(i)<control.hinge_max-control.h_huber_loss){
              fx+= z_cache_col_multiple(i) * (control.hinge_max - margin(i));
              gradient.array() += z_cache_col_multiple(i) * -1 * z_cache.col(i).array();
            }else{
              fx+= z_cache_col_multiple(i) * pow((control.hinge_max+control.h_huber_loss-margin(i)),2)/(4*control.h_huber_loss);
              gradient.array() += z_cache_col_multiple(i) * (margin(i)-control.h_huber_loss-control.hinge_max)/(2*control.h_huber_loss) * z_cache.col(i).array();
            }
          }
        }
        gradient.array() /= z_cache_col_multiple.sum();
        gradient.array() += lambda_regularization;
        gradient.array() *= 2*v.array();
        fx /= z_cache_col_multiple.sum();
        fx += lambda_regularization*v_sq.array().sum();
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      //thanks to Oscar Veliz's youtube + github for this method and examples of code
      Eigen::VectorXd steffensen_fixed_point_acceleration(const Eigen::Ref<const Eigen::MatrixXd> v_i) 
      {
        //v_i needs four columns
        // column (0) is the current point
        // column (1)-(3) are the fixed-point iterations of each previous column
        Eigen::MatrixXd d_v_i(v_i.rows(),2);
        d_v_i.col(0) = v_i.col(1) - v_i.col(0);
        d_v_i.col(1) = v_i.col(2) - v_i.col(1);
        Eigen::MatrixXd d_sq_v_i(v_i.rows(),2);
        d_sq_v_i.col(0) = v_i.col(2) - v_i.col(1) - d_v_i.col(0);
        d_sq_v_i.col(1) = v_i.col(3) - v_i.col(2) - d_v_i.col(1);
        return v_i.col(0) - d_v_i*(d_sq_v_i.colPivHouseholderQr().solve(d_v_i.col(0)));
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      inline double diff_test(const Eigen::Ref<const Eigen::VectorXd> c_base, const Eigen::Ref<const Eigen::VectorXd> c_new)
      {
        return (c_base - c_new).squaredNorm()/std::max(c_base.squaredNorm(),1.0);
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      //update v with v_plus_one and drop features below threshold
      void update_v_with_drop(const double& drop_threshold, 
                              std::vector<unsigned>& feature_index, //updates
                              const Eigen::Ref<const Eigen::VectorXd> v, 
                              Eigen::VectorXd& v_plus_one) //updates
      {
        //identify features with v and v_plus_one below threshold
        Eigen::VectorXd check_vector = v.cwiseMax(v_plus_one);
        std::vector<unsigned> drop_features_index;
        for(unsigned i = 0; i<check_vector.size();i++){
          if(check_vector(i)<drop_threshold){
            drop_features_index.push_back(i);
          }
        }
        if(drop_features_index.size()>0){
          //start dropping
          if(feature_index.size()<=drop_features_index.size()){
            feature_index.clear();
          }
          else{
            for(int i = drop_features_index.size()-1; i>=0; i--){
              if(drop_features_index.at(i) < feature_index.size()-1){
                feature_index.at(drop_features_index.at(i)) = feature_index.back();
                v_plus_one(drop_features_index.at(i)) = v_plus_one(feature_index.size()-1);
              }
              feature_index.pop_back();
            }
          }
          v_plus_one.conservativeResize(feature_index.size());
        }
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      double calc_null_BIC(const Response_Labels& response_labels)
      {
        int keep_y_cols = response_labels.get_y_labels(false).cols();
        int keep_y_size = response_labels.get_y_labels(false).rows()*keep_y_cols;
        double y_lab_sum = response_labels.get_y_labels(false).rightCols(keep_y_cols).sum();
        double y_lab_mean = y_lab_sum/(keep_y_size);
        double logli = std::log(y_lab_mean)*y_lab_sum+std::log(1.0-y_lab_mean)*(keep_y_size-y_lab_sum);
        logli /= (double) keep_y_cols;
        
        return log((double)response_labels.get_y_labels(false).rows()) * 1.0/keep_y_cols - 2*logli;
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      double find_elbow(const Eigen::Ref<const Eigen::VectorXd> vals) 
      {
        std::vector<unsigned> val_sort = sort_indexes(vals);
        if(vals.size()<=2){
          throw std::invalid_argument("for deker::math::find_elbow, vals.size() must be > 2");
        }
        Eigen::MatrixXd val_mat(vals.size(),2);
        //val_mat.col(1) = vals;
        for(unsigned i = 0;i<val_mat.rows();i++){
          val_mat(i,0) = (double) i;
          val_mat(i,1) = vals(val_sort.at(i));
        }
        Eigen::RowVectorXd line_vec_n = val_mat.bottomRows(1)-val_mat.topRows(1);
        line_vec_n.array() /= line_vec_n.norm();
        val_mat = (val_mat.rowwise() - val_mat.row(0)).eval();
        Eigen::VectorXd dot_product = (val_mat.array().rowwise() * line_vec_n.array()).rowwise().sum();
        val_mat -= (dot_product * line_vec_n);
        Eigen::VectorXd dist_to_line = val_mat.rowwise().squaredNorm();
        unsigned max_val_index = 0;
        double max_dist = 0;
        double tolerance = 1e-8;
        for(int i = dist_to_line.size()-2; i>=0; i--){
          if(val_mat(i,0)>=0&&val_mat(i,1)<=0&&
             dist_to_line(i)-max_dist>tolerance)
          {
            max_val_index = (unsigned) i;
            max_dist = dist_to_line(i);
          }
        }
        return vals(val_sort.at(max_val_index));
      }
    }
  }
}
#endif