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
#include <deker/fit.h>
#include <deker/io/control.h>
#include <deker/io/misc.h>
#include <deker/io/opt_param.h>
#include <deker/io/output.h>
//
#include <unistd.h>
////////////////////////////////////////////////////////////////////////////////////////////////
//main program control
int main(int argc,  char **argv){
  //setup for options
  //options for binary conversion
  bool binary_convert = false;
  //options for compiling edge list from summary file
  bool make_edge_list = false;
  //options for doing fit
  bool do_fit = false;
  bool do_fit_manual = false;
  //options for finding sigma
  bool start_control = false;
  //options for writing output
  std::string output_file;
  /////////////////////////////////////
  /////////////////////////////////////
  //parse options
  int opt;
  while((opt = getopt(argc, argv, "befsw:")) != -1){  
    switch(opt){
    case 'b': //switch to use deker_job to convert a csv file to binary
      binary_convert=true;
      break;
    case 'e': //switch for creating edgelist output
      make_edge_list = true;
      break;
    case 'f': //switch for doing fit 
      do_fit = true;
      break;
    case 's': //switch for starting opt control
      start_control = true;
      break;
    case 'w': //where to write output
      output_file = optarg;
      break;
    }  
  }
  /////////////////////////////////////
  /////////////////////////////////////
  if(binary_convert){
    /////////////////////////////////////
    if(argc-optind<2)
      throw std::invalid_argument("At least 2 non-option arguments are required: data file location & location to write binary file\n");
    deker::io::Eigen_Matrix_IO to_convert(argv[optind],false);
    to_convert.write_to_binary(argv[optind+1]);
    
    /////////////////////////////////////
  }
  else if(start_control){
    /////////////////////////////////////
    if(argc-optind<1)
      throw std::invalid_argument("At least 1 non-option argument is required: control file location\n");
    if(!deker::io::file_exists_test(argv[optind]))
      throw std::invalid_argument("Control file not found at provided location\n");
    deker::io::Control_Base control(argv[optind]);
    deker::io::Eigen_Matrix_IO sizes_only_data(control.data_file,true,true);
    std::vector<unsigned> use_response = deker::io::build_use_vectors(control.ex_res,control.in_res,(unsigned) sizes_only_data.cols);
    std::cout<<use_response.size()<<"\n";
    /////////////////////////////////////
  }
  else if(make_edge_list){
    /////////////////////////////////////
    //deker_edge_list(argc, argv, optind, dump_all_models);
    if(argc-optind<2)
      throw std::invalid_argument("At least 2 non-option arguments are required: output file location & location to edge list file\n");
    //
    deker::fit::Output_deker_IO to_write(argv[optind]);
    std::vector<std::string> extra_id_header;
    extra_id_header.push_back("edgelist_filename");
    extra_id_header.push_back("source_filename");
    std::vector<std::string> extra_id;
    extra_id.push_back(argv[optind+1]);
    extra_id.push_back(argv[optind]);
    //
    to_write.write_to_csv(argv[optind+1],extra_id_header,extra_id);
    /////////////////////////////////////
  }
  else if(do_fit){
    /////////////////////////////////////
    //values needed to run fit
    if(argc-optind<2)
      throw std::invalid_argument("2 non-option arguments are required: which response to use and control file");
    deker::io::Control_Base control(argv[optind+1]);
    deker::fit::Opt_Param control_params = (deker::fit::Opt_Param) control.opt_param;
    deker::io::Eigen_Matrix_IO input_data(control.data_file);
    std::vector<unsigned> use_response = deker::io::build_use_vectors(control.ex_res,control.in_res,(unsigned) input_data.cols);
    unsigned input_response(deker::io::read_arg_to_val <unsigned> (argv[optind]));
    unsigned which_response = use_response.at(input_response);
    std::vector<unsigned> use_predictors = deker::io::build_use_vectors(control.ex_pred,control.in_pred,(unsigned) input_data.cols);
    //if use_predictors is empty, make a vector of indicies for data columns - 0:(input_data.matrix.cols()-1)
    if(use_predictors.size()==0){
      for(unsigned i = 0; i<input_data.cols; i++)
        use_predictors.push_back(i);
    }
    //do the fitting
    deker::fit::Solve_deker sol_deker(input_data.matrix,
                                      use_predictors,
                                      which_response,
                                      control_params);
    deker::fit::Output_deker opt_deker_sol;
    sol_deker.opt_lambda(opt_deker_sol);
    //write output to file 
    deker::fit::Output_deker_IO to_write(opt_deker_sol,
                                         input_data.feature_names,
                                         control.data_file,
                                         which_response);
    to_write.write_binary(output_file);
    /////////////////////////////////////
  }else{
    /////////////////////////////////////
    std::cout<<"usage/help TODO\n";
    /////////////////////////////////////
  }
  /////////////////////////////////////
  /////////////////////////////////////
  return 0;
}