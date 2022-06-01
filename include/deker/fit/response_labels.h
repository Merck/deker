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
#ifndef DEKER_FIT_RESPONSE_LABELS
#define DEKER_FIT_RESPONSE_LABELS
///////////////////////////////////////////
#include <deker/io/crypto.h>
//
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <Eigen/LU>
#include <vector>
#include <fstream>
#include <iostream>
///////////////////////////////////////////
namespace deker{
  namespace fit{
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    class Response_Labels{
    protected:
      std::string parent_hash;
      const Eigen::Ref<const Eigen::VectorXd> response;
      std::vector<unsigned> threshold_inds;
      std::vector<unsigned> threshold_mult;
      Eigen::MatrixXd y_labels; //samples x samples
      Eigen::MatrixXd y_labels_pseudoinv;
      bool use_fix_labs;
      bool class_fix;
      Eigen::MatrixXd p_y_labels;
      std::vector<unsigned> p_threshold_inds;
      std::vector<unsigned> p_threshold_mult;
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
    public:
      Response_Labels(const Eigen::Ref<const Eigen::VectorXd> response):response(response){}
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      Response_Labels(const Eigen::Ref<const Eigen::VectorXd> response, bool dummy_build):response(response)
      {
        build();
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      //!! Response_Labels assumes & requires that response be sorted least->greatest
      //void build(const Eigen::Ref<const Eigen::VectorXd> response)
      void build()
      {
        parent_hash = crypto::SHA256_Eigen(response);
        //find number of classifiers to build (number of unique values in response)
        //ignore the lowest point - will not have two classes to split on
        //threshold_vals.push_back(response(0));
        threshold_inds.push_back(0);
        threshold_mult.push_back(1);
        class_fix = false;
        use_fix_labs = false;
        double last_threshold_val = response(0);
        for(unsigned i = 1; i<response.size();i++){
          if(response(i)>last_threshold_val){
            last_threshold_val = response(i);
            threshold_inds.push_back(i);
            threshold_mult.push_back(1);
          }else{
            threshold_mult.back()++;
            class_fix = true;
          }
        }
        //thresholds start with lowest, which will not have a valid split (all >=)
        //so skip lowest split in y_labels
        y_labels = Eigen::MatrixXd::Zero(response.size(),threshold_inds.size());
        for(unsigned i = 0;i<threshold_inds.size();i++){
          for(unsigned j = threshold_inds.at(i);j<response.size();j++){
            y_labels(j,i) = 1.0;
          }
        }
        Eigen::MatrixXd* y_lab_ptr = &y_labels;
        if(class_fix){
          p_y_labels = Eigen::MatrixXd::Zero(response.size(),response.size());
          for(unsigned i = 0;i<p_y_labels.cols();i++){
            p_threshold_inds.push_back(i);
            p_threshold_mult.push_back(1);
            for(unsigned j = i;j<p_y_labels.rows();j++){
              p_y_labels(j,i) = 1.0;
            }
          }
          y_lab_ptr = &p_y_labels;
        }
        //guarunteed that y_lab_ref is invertible
        y_labels_pseudoinv = (*y_lab_ptr).colPivHouseholderQr().solve(Eigen::MatrixXd::Identity((*y_lab_ptr).cols(),(*y_lab_ptr).cols()));
        /////////////////////////////
        //drop first column
        y_labels.leftCols(threshold_inds.size()-1) = y_labels.rightCols(threshold_inds.size()-1).eval();
        y_labels.conservativeResize(response.size(),threshold_inds.size()-1);
        if(class_fix){
          p_y_labels.leftCols(response.size()-1) = p_y_labels.rightCols(response.size()-1).eval();
          p_y_labels.conservativeResize(response.size(),response.size()-1);
        }
        /////////////////////////////
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      const Eigen::Ref<const Eigen::VectorXd> get_response() const {return response;}
      const std::vector<unsigned>& get_threshold_inds(const bool& use_fix) const
      {
        if(use_fix&&class_fix){
          return p_threshold_inds;
        }else{
          return threshold_inds;
        }
      }
      const std::vector<unsigned>& get_threshold_mult(const bool& use_fix) const
      {
        if(use_fix&&class_fix){
          return p_threshold_mult;
        }else{
          return threshold_mult;
        }
      }
      const Eigen::Ref<const Eigen::MatrixXd> get_y_labels(const bool& use_fix) const
      {
        if(use_fix&&class_fix){
          return p_y_labels;
        }else{
          return y_labels;
        }
      }
      const Eigen::Ref<const Eigen::MatrixXd> get_y_labels_pseudoinv() const {return y_labels_pseudoinv;}
      const bool& get_class_fix() const {return class_fix;}
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      void serialize(std::ifstream& infile)
      {
        //
        threshold_inds.clear();
        threshold_mult.clear();
        p_threshold_inds.clear();
        p_threshold_mult.clear();
        //read in parent_hash
        size_t tmp_size;
        infile.read((char*) (&tmp_size), sizeof(tmp_size));
        parent_hash.resize(tmp_size);
        infile.read(&parent_hash[0],tmp_size);
        //check that response matches parent hash - using same data as the Response_Labels that was written from
        std::string input_data_hash(crypto::SHA256_Eigen(response));
        if(input_data_hash!=parent_hash){
          throw std::invalid_argument("input_data hash does not match parent_hash");
        }
        //proceed with read-in as usual
        
        //read in bools
        infile.read((char*) (&class_fix), sizeof(class_fix));
        infile.read((char*) (&use_fix_labs), sizeof(use_fix_labs));
        //read in threshold_inds & threshold_mults
        infile.read((char*) (&tmp_size),sizeof(tmp_size));
        unsigned tmp_index;
        for(size_t i = 0; i<tmp_size; i++){
          infile.read((char*) (&tmp_index), sizeof(tmp_index));
          threshold_inds.push_back(tmp_index);
          infile.read((char*) (&tmp_index), sizeof(tmp_index));
          threshold_mult.push_back(tmp_index);
        }
        //read in p_threshold_inds & p_threshold_mults
        infile.read((char*) (&tmp_size),sizeof(tmp_size));
        for(size_t i = 0; i<tmp_size; i++){
          infile.read((char*) (&tmp_index), sizeof(tmp_index));
          p_threshold_inds.push_back(tmp_index);
          infile.read((char*) (&tmp_index), sizeof(tmp_index));
          p_threshold_mult.push_back(tmp_index);
        }
        //
        Eigen::MatrixXd::Index rows=0, cols=0;
        //read in y_labels
        infile.read((char*) (&rows),sizeof(Eigen::MatrixXd::Index));
        infile.read((char*) (&cols),sizeof(Eigen::MatrixXd::Index));
        y_labels.resize(rows, cols);
        infile.read( (char *) y_labels.data() , rows*cols*sizeof(Eigen::MatrixXd::Scalar));
        //read in p_y_labels
        infile.read((char*) (&rows),sizeof(Eigen::MatrixXd::Index));
        infile.read((char*) (&cols),sizeof(Eigen::MatrixXd::Index));
        p_y_labels.resize(rows, cols);
        infile.read( (char *) p_y_labels.data() , rows*cols*sizeof(Eigen::MatrixXd::Scalar));
        //read in y_labels_pseudoinv
        infile.read((char*) (&rows),sizeof(Eigen::MatrixXd::Index));
        infile.read((char*) (&cols),sizeof(Eigen::MatrixXd::Index));
        y_labels_pseudoinv.resize(rows, cols);
        infile.read( (char *) y_labels_pseudoinv.data() , rows*cols*sizeof(Eigen::MatrixXd::Scalar));
        // //read in y_labels_pseudoident
        // infile.read((char*) (&rows),sizeof(Eigen::MatrixXd::Index));
        // infile.read((char*) (&cols),sizeof(Eigen::MatrixXd::Index));
        // y_labels_pseudoident.resize(rows, cols);
        // infile.read( (char *) y_labels_pseudoident.data() , rows*cols*sizeof(Eigen::MatrixXd::Scalar));
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      void serialize(std::ofstream& outfile) const
      {
        size_t tmp_size = parent_hash.size();
        outfile.write((char*) (&tmp_size), sizeof(tmp_size));
        outfile.write(parent_hash.c_str(),tmp_size);
        //write bools
        outfile.write((char*) (&class_fix), sizeof(class_fix));
        outfile.write((char*) (&use_fix_labs), sizeof(use_fix_labs));
        //write threshold_inds & threshold_mults
        size_t index_size = threshold_inds.size();
        outfile.write((char*) (&index_size),sizeof(index_size));
        for(size_t i = 0; i<index_size; i++){
          outfile.write((char*) (&threshold_inds.at(i)),sizeof(unsigned));
          outfile.write((char*) (&threshold_mult.at(i)),sizeof(unsigned));
        }
        //write p_threshold_inds & p_threshold_mults
        index_size = p_threshold_inds.size();
        outfile.write((char*) (&index_size),sizeof(index_size));
        for(size_t i = 0; i<index_size; i++){
          outfile.write((char*) (&p_threshold_inds.at(i)),sizeof(unsigned));
          outfile.write((char*) (&p_threshold_mult.at(i)),sizeof(unsigned));
        }
        //write out y_labels
        Eigen::MatrixXd::Index rows=y_labels.rows(), cols=y_labels.cols();
        outfile.write((char*) (&rows), sizeof(Eigen::MatrixXd::Index));
        outfile.write((char*) (&cols), sizeof(Eigen::MatrixXd::Index));
        outfile.write((char*) y_labels.data(), rows*cols*sizeof(Eigen::MatrixXd::Scalar));
        //write out p_y_labels
        rows=p_y_labels.rows();
        cols=p_y_labels.cols();
        outfile.write((char*) (&rows), sizeof(Eigen::MatrixXd::Index));
        outfile.write((char*) (&cols), sizeof(Eigen::MatrixXd::Index));
        outfile.write((char*) p_y_labels.data(), rows*cols*sizeof(Eigen::MatrixXd::Scalar));
        //write out y_labels_pseudoinv
        rows = y_labels_pseudoinv.rows();
        cols = y_labels_pseudoinv.cols();
        outfile.write((char*) (&rows), sizeof(Eigen::MatrixXd::Index));
        outfile.write((char*) (&cols), sizeof(Eigen::MatrixXd::Index));
        outfile.write((char*) y_labels_pseudoinv.data(), rows*cols*sizeof(Eigen::MatrixXd::Scalar));
        // //write out y_labels_pseudoident
        // rows = y_labels_pseudoident.rows();
        // cols = y_labels_pseudoident.cols();
        // outfile.write((char*) (&rows), sizeof(Eigen::MatrixXd::Index));
        // outfile.write((char*) (&cols), sizeof(Eigen::MatrixXd::Index));
        // outfile.write((char*) y_labels_pseudoident.data(), rows*cols*sizeof(Eigen::MatrixXd::Scalar));
      }
    };
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
  }
}
#endif