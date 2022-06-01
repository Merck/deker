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
#ifndef DEKER_FIT_SPLIT_DATA
#define DEKER_FIT_SPLIT_DATA
///////////////////////////////////////////
#include <deker/io/crypto.h>
#include <Eigen/Core>
#include <vector>
#include <fstream>
///////////////////////////////////////////
namespace deker{
  namespace fit{
  /////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////
  //stolen from stackoverflow, credit Lukasz Wiklendt
  std::vector<unsigned> sort_indexes(const std::vector<double>& v) {
    // initialize original index locations
    std::vector<unsigned> idx(v.size());
    for(unsigned i = 0; i<idx.size();i++) idx.at(i) = i;
    // sort indexes based on comparing values in v
    sort(idx.begin(), idx.end(), [&v](unsigned i1, unsigned i2) {return v.at(i1) < v.at(i2);});
    return idx;
  }
  /////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////
  //stolen from stackoverflow, credit Lukasz Wiklendt
  std::vector<unsigned> sort_indexes(const Eigen::Ref<const Eigen::VectorXd> v) {
    // initialize original index locations
    std::vector<unsigned> idx(v.size());
    for(unsigned i = 0; i<idx.size();i++) idx.at(i) = i;
    // sort indexes based on comparing values in v
    sort(idx.begin(), idx.end(), [&v](unsigned i1, unsigned i2) {return v(i1) < v(i2);});
    return idx;
  }
  /////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////
  struct Split_Data{
    std::string parent_hash;
    Eigen::MatrixXd& predictors;
    Eigen::VectorXd response;
    unsigned which_response;
    std::vector<unsigned> original_response_sort_index;
    std::vector<unsigned> predictors_index;
    std::vector<unsigned> key_to_original_data;
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    Split_Data(Eigen::MatrixXd& input_data, 
               unsigned which_response):
      predictors(input_data),
      which_response(which_response),
      response(predictors.col(which_response)){}
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    Split_Data(Eigen::MatrixXd& input_data,
               unsigned which_response,
               const std::vector<unsigned>& predictors_use_index):
      predictors(input_data),
      which_response(which_response),
      response(predictors.col(which_response))
    {
      build(predictors_use_index);
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    void build(const std::vector<unsigned>& predictors_use_index)
    {
      predictors_index = predictors_use_index;
      parent_hash = crypto::SHA256_Eigen(predictors);
      //
      //populate original data key - current order of 'predictors'
      key_to_original_data.clear();
      for(unsigned i = 0; i<predictors.cols();i++){
        key_to_original_data.push_back(i);
      }
      //get index to resort
      original_response_sort_index = sort_indexes(predictors.col(which_response));
      //resort
      sort_predictors(original_response_sort_index);
      //extract response 
      extract_response(which_response);
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    void sort_predictors(const std::vector<unsigned>& sort_index)
    {
      std::vector<unsigned> position_index;
      for(unsigned i = 0; i<sort_index.size();i++){
        position_index.push_back(i);
      }
      std::vector<unsigned> element_index;
      for(unsigned i = 0; i<sort_index.size();i++){
        element_index.push_back(i);
      } 
      for(unsigned move_to = 0; move_to<sort_index.size(); move_to++){
        unsigned move_from = position_index.at(sort_index.at(move_to));
        if(move_to!=move_from){
          std::swap(position_index.at(element_index.at(move_to)),position_index.at(element_index.at(move_from)));
          std::swap(element_index.at(move_to),element_index.at(move_from));
          predictors.row(move_to).swap(predictors.row(move_from));
        }
      } 
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    void extract_response(const unsigned& which_response)
    {
      response = predictors.col(which_response);
      for(unsigned i = 0; i<predictors_index.size();i++){
        if(which_response==predictors_index.at(i)){
          predictors_index.at(i) = predictors_index.back();
          predictors_index.pop_back();
          break;
        }
      }
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    void subset_predictors(const std::vector<unsigned>& keep_predictors)
    {
      /////////////////
      std::vector<unsigned> position_index;
      for(unsigned i = 0; i<predictors.cols();i++){
        position_index.push_back(i);
      }
      std::vector<unsigned> element_index;
      for(unsigned i = 0; i<predictors.cols();i++){
        element_index.push_back(i);
      }
      /////////////////
      for(unsigned move_to = 0;  move_to<keep_predictors.size(); move_to++){
        unsigned move_from = position_index.at(keep_predictors.at(move_to));
        if(move_to!=move_from){
          std::swap(position_index.at(element_index.at(move_to)),position_index.at(element_index.at(move_from)));
          std::swap(element_index.at(move_to),element_index.at(move_from));
          predictors.col(move_to).swap(predictors.col(move_from));
        }
      }
      predictors.conservativeResize(predictors.rows(),keep_predictors.size());
      std::vector<unsigned> new_data_key;
      predictors_index.clear();
      for(unsigned i = 0; i<keep_predictors.size(); i++){
        new_data_key.push_back(key_to_original_data.at(keep_predictors.at(i)));
        predictors_index.push_back(i);
      }
      key_to_original_data = new_data_key;
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    void serialize(std::ifstream& infile)
    {
      original_response_sort_index.clear();
      predictors_index.clear();
      key_to_original_data.clear();
      //
      size_t temp_name_size;
      size_t tmp_size;
      unsigned tmp_index;
      //read in parent_hash
      infile.read((char*) (&temp_name_size), sizeof(temp_name_size));
      parent_hash.resize(temp_name_size);
      infile.read(&parent_hash[0],temp_name_size);
      //check that input_data matches parent hash - using same data as the Split_Data that was written from
      std::string input_data_hash(crypto::SHA256_Eigen(predictors));
      if(input_data_hash!=parent_hash){
        throw std::invalid_argument("input_data hash does not match parent_hash");
      }
      //read in which_response
      unsigned tmp_which_response;
      infile.read((char*) (&tmp_which_response),sizeof(tmp_which_response));
      //check that which_response matches input
      if(tmp_which_response!=which_response){
        throw std::invalid_argument("which_response does not match read in value");
      }
      //read in original_response_sort_index
      infile.read((char*) (&tmp_size),sizeof(tmp_size));
      for(unsigned i = 0; i<tmp_size;i++){
        infile.read((char*) (&tmp_index), sizeof(tmp_index));
        original_response_sort_index.push_back(tmp_index);
      }
      //read in predictors_index
      infile.read((char*) (&tmp_size),sizeof(tmp_size));
      for(unsigned i = 0; i<tmp_size;i++){
        infile.read((char*) (&tmp_index), sizeof(tmp_index));
        predictors_index.push_back(tmp_index);
      }
      //read in key_to_original_data
      infile.read((char*) (&tmp_size),sizeof(tmp_size));
      for(unsigned i = 0; i<tmp_size;i++){
        infile.read((char*) (&tmp_index), sizeof(tmp_index));
        key_to_original_data.push_back(tmp_index);
      }
      ///////////////////////////////////////////////////
      //rebuild by formatting predictors
      //make response/predictors from input_data
      sort_predictors(original_response_sort_index);
      extract_response(which_response);
      //
      std::vector<unsigned> new_key = key_to_original_data;
      key_to_original_data.clear();
      for(unsigned i = 0; i<predictors.cols();i++){
        key_to_original_data.push_back(i);
      }
      subset_predictors(new_key);
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    void serialize(std::ofstream& outfile) const
    {
      //write out parent_hash
      size_t tmp_size = parent_hash.size();
      outfile.write((char*) (&tmp_size), sizeof(tmp_size));
      outfile.write(parent_hash.c_str(),tmp_size);
      //write out which_response
      outfile.write((char*) (&which_response),sizeof(which_response));
      //read in original_response_sort_index
      tmp_size = original_response_sort_index.size();
      outfile.write((char*) (&tmp_size),sizeof(tmp_size));
      for(unsigned i = 0; i<tmp_size;i++){
        outfile.write((char*) (&original_response_sort_index.at(i)),sizeof(unsigned));
      }
      //write out predictors_index
      tmp_size = predictors_index.size();
      outfile.write((char*) (&tmp_size),sizeof(tmp_size));
      for(unsigned i = 0; i<tmp_size;i++){
        outfile.write((char*) (&predictors_index.at(i)),sizeof(unsigned));
      }
      //write out key_to_original_data
      tmp_size = key_to_original_data.size();
      outfile.write((char*) (&tmp_size),sizeof(tmp_size));
      for(unsigned i = 0; i<tmp_size;i++){
        outfile.write((char*) (&key_to_original_data.at(i)),sizeof(unsigned));
      }
    }
  };
  /////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////
  }
}
#endif