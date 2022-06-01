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
#ifndef DEKER_IO_CONTROL
#define DEKER_IO_CONTROL
#include <deker/fit/opt_param.h>
#include <deker/io/misc.h>
////////////////////////////////
namespace deker{
  namespace io{
  /////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////
    struct Control_Base{
      fit::Opt_Param opt_param;
      std::vector<unsigned> in_res, ex_res, in_pred, ex_pred;
      std::vector<unsigned> exclude_sample;
      std::string data_file;
      int lambda_iter_per_job;
      double max_job_time;
      ///////////////////////////////////////////////////////
      Control_Base()
      {
        lambda_iter_per_job = -1;
        max_job_time = -1;
      }
      ///////////////////////////////////////////////////////
      void read_from_txt_file(const std::string& input_file)
      {
        if(!file_exists_test(input_file)||input_file.empty())
          throw std::invalid_argument("input file not found\n");
        in_res.clear();
        ex_res.clear();
        in_pred.clear();
        ex_pred.clear();
        std::ifstream source;
        source.open(input_file);
        std::string line;
        while(std::getline(source,line)){
          std::size_t arg_pos = line.find(" ");
          std::string identifier = line.substr(0,arg_pos);
          std::string arg = line.substr(arg_pos+1);
          if(identifier.compare("lambda_init_max") == 0){
            opt_param.lambda_init_max = read_arg_to_val<double>(arg);
          }
          else if(identifier.compare("lambda_init_min") == 0){
            opt_param.lambda_init_min = read_arg_to_val<double>(arg);
          }
          else if(identifier.compare("lambda_init_stepsize") == 0){
            opt_param.lambda_init_stepsize = read_arg_to_val<double>(arg);
          }
          else if(identifier.compare("lambda_init_maxiter") == 0){
            opt_param.lambda_init_maxiter = read_arg_to_val<int>(arg);
          }
          else if(identifier.compare("lambda_fill_mindist") == 0){
            opt_param.lambda_fill_mindist = read_arg_to_val<double>(arg);
          }
          else if(identifier.compare("lambda_fill_maxiter") == 0){
            opt_param.lambda_fill_maxiter = read_arg_to_val<int>(arg);
          }
          else if(identifier.compare("lambda_bo_maxiter") == 0){
            opt_param.lambda_bo_maxiter = read_arg_to_val<int>(arg);
          }
          else if(identifier.compare("h_huber_loss") == 0){
            opt_param.h_huber_loss = read_arg_to_val<double>(arg);
          }
          else if(identifier.compare("w_convergence_threshold") == 0){
            opt_param.w_convergence_threshold = read_arg_to_val<double>(arg);
          }
          else if(identifier.compare("w_maxiter") == 0){
            opt_param.w_maxiter = read_arg_to_val<double>(arg);
          }
          else if(identifier.compare("feature_drop_threshold") == 0){
            opt_param.feature_drop_threshold = read_arg_to_val<double>(arg);
          }
          else if(identifier.compare("data_file") == 0){
            data_file = arg;
          }
          else if(identifier.compare("include_response") == 0){
            in_res = read_arg_to_vector<unsigned>(arg);
          }
          else if(identifier.compare("exclude_response") == 0){
            ex_res = read_arg_to_vector<unsigned>(arg);
          }
          else if(identifier.compare("include_predictor") == 0){
            in_pred = read_arg_to_vector<unsigned>(arg);
          }
          else if(identifier.compare("exclude_predictor") == 0){
            ex_pred = read_arg_to_vector<unsigned>(arg);
          }
          else if(identifier.compare("exclude_sample") == 0){
            exclude_sample = read_arg_to_vector<unsigned>(arg);
          }
          else if(identifier.compare("lambda_iter_per_job") == 0){
            lambda_iter_per_job = read_arg_to_val<int>(arg);
          }
          else if(identifier.compare("max_job_time") == 0){
            max_job_time = read_arg_to_val<double>(arg);
          }
        }
      }
      ///////////////////////////////////////////////////////
      void serialize(std::ifstream& infile)
      {
        //read data_file
        size_t data_file_size;
        infile.read((char*) (&data_file_size), sizeof(data_file_size));
        data_file.resize(data_file_size);
        infile.read(&data_file[0], data_file_size);
        //response include/exclude
        serialize_vector<unsigned>(infile,in_res);
        serialize_vector<unsigned>(infile,ex_res);
        //predictor include/exclude
        serialize_vector<unsigned>(infile,in_pred);
        serialize_vector<unsigned>(infile,ex_pred);
        //exclude sample
        serialize_vector<unsigned>(infile,exclude_sample);
        //lambda_iter_per_job
        infile.read((char*) (&lambda_iter_per_job), sizeof(lambda_iter_per_job));
        //max_job_time
        infile.read((char*) (&max_job_time), sizeof(max_job_time));
        //read h_opt & opt param structs
        opt_param.serialize(infile);
      }
      ///////////////////////////////////////////////////////
      void serialize(std::ofstream& outfile) const
      {
        //write data_file 
        size_t data_file_size = data_file.size();
        outfile.write((char*) (&data_file_size),sizeof(data_file_size));
        outfile.write(data_file.c_str(), data_file_size);
        //response include/exclude
        serialize_vector<unsigned>(outfile,in_res);
        serialize_vector<unsigned>(outfile,ex_res);
        //predictor include/exclude
        serialize_vector<unsigned>(outfile,in_pred);
        serialize_vector<unsigned>(outfile,ex_pred);
        //exclude sample
        serialize_vector<unsigned>(outfile,exclude_sample);
        //lambda_iter_per_job
        outfile.write((char*) (&lambda_iter_per_job),sizeof(lambda_iter_per_job));
        //max_job_time
        outfile.write((char*) (&max_job_time), sizeof(max_job_time));
        //write h_opt & opt param structs
        opt_param.serialize(outfile);
      }
      ///////////////////////////////////////////////////////
      void text(std::ofstream& outfile) const
      {
        outfile<<"h_huber_loss"<<" "<<opt_param.h_huber_loss<<"\n";
        outfile<<"w_convergence_threshold"<<" "<<opt_param.w_convergence_threshold<<"\n";
        outfile<<"w_maxiter"<<" "<<opt_param.w_maxiter<<"\n";
        outfile<<"feature_drop_threshold"<<" "<<opt_param.feature_drop_threshold<<"\n";
        outfile<<"data_file"<<" "<<data_file<<"\n";
        outfile<<"lambda_iter_per_job"<<" "<<lambda_iter_per_job<<"\n";
        outfile<<"max_job_time"<<" "<<max_job_time<<"\n";
        outfile<<"include_response"<<" "; write_vector_to_csv(outfile,in_res); outfile<<"\n";
        outfile<<"exclude_response"<<" "; write_vector_to_csv(outfile,ex_res); outfile<<"\n";
        outfile<<"include_predictor"<<" "; write_vector_to_csv(outfile,in_pred); outfile<<"\n";
        outfile<<"exclude_predictor"<<" "; write_vector_to_csv(outfile,ex_pred); outfile<<"\n";
        outfile<<"exclude_sample"<<" "; write_vector_to_csv(outfile,exclude_sample); outfile<<"\n";
      }
      ///////////////////////////////////////////////////////
      void write_control_file(const std::string& output_file) const
      {
        std::ofstream outfile(output_file,std::ios_base::app);
        text(outfile);
        outfile.close();
      }
    };
    /////////////////////////////////////////////////////
    //Credit stackoverflow user DShook
    template<typename T>
    void remove_duplicates(std::vector<T>& vec)
    {
      std::sort(vec.begin(), vec.end());
      vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
    }
    ////////////////////////////////////////////////////
    std::vector<unsigned> build_use_vectors(const std::vector<unsigned>& exclude_elements, const std::vector<unsigned>& include_elements,const unsigned& total_length=0)
    {
      if(include_elements.size()==0&&total_length==0)
        throw std::invalid_argument("if include_elements is empty, total length must be given and nonzero\n");
      std::vector<unsigned> to_use(include_elements);
      std::vector<unsigned> to_remove(exclude_elements);
      //if include is empty, fill with a vector of 0:(n_features-1)
      if(to_use.empty()){
        for(unsigned i = 0; i<total_length; i++)
          to_use.push_back(i);
      }
      else{
        remove_duplicates<unsigned>(to_use);
      }
      if(!to_remove.empty()){
        remove_duplicates<unsigned>(to_remove);
        //iterate through to_use, removing any value in exclude
        unsigned i=0;
        while(i<to_use.size()){
          bool increase_i = true;
          for(unsigned j = 0; j<to_remove.size(); j++){
            if(to_use.at(i)==to_remove.at(j)){
              //remove value from to_remove once it has been found
              if(j<to_remove.size()-1)
                to_remove.at(j)=to_remove.back();
              to_remove.pop_back();
              //remove value from to_use once it has been found
              if(i<to_use.size()-1)
                to_use.at(i)=to_use.back();
              to_use.pop_back();
              increase_i = false;
              break;
            }
          }
          if(increase_i) i++;
        }
      }
      return(to_use);
    }
  ////////////////////////////////////////////////////
  }
}
#endif