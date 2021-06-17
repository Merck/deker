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
#ifndef DEKER_IO_MISC
#define DEKER_IO_MISC
#include <Eigen/Core>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
////////////////////////////////
namespace deker {
  namespace io {
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    //From stackoverflow, thanks PherricOxide
    inline bool file_exists_test (const std::string& name) {
      if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
      } else {
        return false;
      }   
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    template <typename T>
    void serialize_vector(std::ifstream& infile,std::vector<T>& to_read)
    {
      to_read.clear();
      size_t size;
      infile.read((char*) (&size),sizeof(size));
      for(size_t i=0;i<size;i++){
        T temp_val;
        infile.read((char*) (&temp_val), sizeof(T));
        to_read.push_back(temp_val);
      }
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    template <typename T>
    void serialize_vector(std::ofstream& outfile,const std::vector<T>& to_write)
    {
      size_t size;
      size = to_write.size();
      outfile.write((char*) (&size),sizeof(size));
      for(size_t i=0;i<size;i++)
        outfile.write((char*) (&to_write.at(i)), sizeof(T));
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    struct Eigen_Matrix_IO{
      Eigen::MatrixXd::Index rows, cols;
      Eigen::MatrixXd matrix;
      std::vector<std::string> feature_names;
      ///////////////////////////////
      void serialize(std::ifstream& infile, const bool& sizes_only=false)
      {
        //Eigen::MatrixXd::Index rows=0, cols=0;
        infile.read((char*) (&rows),sizeof(Eigen::MatrixXd::Index));
        infile.read((char*) (&cols),sizeof(Eigen::MatrixXd::Index));
        if(!sizes_only){
          matrix.resize(rows, cols);
          infile.read( (char *) matrix.data() , rows*cols*sizeof(Eigen::MatrixXd::Scalar) );
          //read in row names data
          feature_names.clear();
          size_t size;
          std::string str;
          while(infile.read((char*) (&size), sizeof(size))){
            str.resize(size);
            infile.read(&str[0], size);
            feature_names.push_back(str);
          }
        }
      }
      ///////////////////////////////
      void serialize(std::ofstream& outfile)
      {
        //Eigen::MatrixXd::Index rows=matrix.rows(), cols=matrix.cols();
        outfile.write((char*) (&rows), sizeof(Eigen::MatrixXd::Index));
        outfile.write((char*) (&cols), sizeof(Eigen::MatrixXd::Index));
        outfile.write((char*) matrix.data(), rows*cols*sizeof(Eigen::MatrixXd::Scalar));
        //write in row names data
        for(unsigned i = 0; i<feature_names.size();i++){
          size_t size=feature_names.at(i).size();
          outfile.write((char*) (&size),sizeof(size));
          outfile.write(feature_names.at(i).c_str(), size);
        }
      }
      ///////////////////////////////
      void text(std::ifstream& infile)
      {
        std::string line;
        //take first row as column names
        if(std::getline(infile,line)){
          //remove the end of line so we can parse by commas without catching an endline in the last element
          //line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
          if (line[line.length()-1] == '\n'||line[line.length()-1] == '\r') {
            line.erase(line.length()-1);
          }
          std::stringstream stream(line);
          feature_names.clear();
          std::string d;
          while (std::getline(stream,d,',')){
            feature_names.push_back(d);
          }
        }
        //prep the rest for read in to matrix
        std::vector<std::vector<double>> elements;
        while(std::getline(infile,line)){
          std::stringstream stream(line);
          elements.push_back(std::vector<double>());
          double d;
          while (stream >> d){
            elements.back().push_back(d);
            stream.ignore(1, ',');
          }
        }
        rows = elements.size();
        //TODO: Error checking that all lines are equally sized
        cols = elements.back().size();
        matrix.resize(rows,cols);
        for(unsigned i=0; i<matrix.rows();i++)
          for(unsigned j=0;  j<matrix.cols();j++)
            matrix(i,j) = elements.at(i).at(j);
      }
      ///////////////////////////////
      void text(std::ofstream& outfile)
      {
        for(unsigned i=0;i<feature_names.size();i++){
          outfile<<feature_names.at(i);
          if(i<feature_names.size()-1)
            outfile<<",";
          else
            outfile<<"\n";
        }
        
        for(unsigned i=0;i<rows;i++)
          for(unsigned j=0;j<cols;j++){
            outfile<<matrix(i,j);
            if(j<cols-1)
              outfile<<",";
            else
              outfile<<"\n";
          }
      }
      ///////////////////////////////
      Eigen_Matrix_IO(){}
      ///////////////////////////////
      Eigen_Matrix_IO(const std::string& filename, const bool& is_binary=true, const bool& sizes_only=false)
      {
        if(!file_exists_test(filename)||filename.empty())
          throw std::invalid_argument("Data file not found at provided location\n");
        if(is_binary){
          std::ifstream infile(filename, std::ios::in | std::ios::binary);
          serialize(infile,sizes_only);
          infile.close();
        }else{
          std::ifstream infile(filename, std::ios::in);
          text(infile);
          infile.close();
        }
      }
      ///////////////////////////////
      void write_to_binary(const std::string& filename){
        std::ofstream outfile(filename, std::ios::out | std::ios::binary | std::ios::trunc);
        serialize(outfile);
        outfile.close();
      }
      ///////////////////////////////
      void write_to_text(const std::string& filename){
        std::ofstream outfile(filename, std::ios::out | std::ios::binary | std::ios::trunc);
        text(outfile);
        outfile.close();
      }
    };
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    template <class Type>
    std::vector<Type> read_arg_to_vector(const std::string& arg)
    {
      std::vector<Type> to_return;
      std::stringstream stream(arg);
      Type u;
      while (stream >> u){
        to_return.push_back(u);
        stream.ignore(1, ',');
      }
      return to_return;
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    template <typename Type>
    void write_vector_to_csv(std::ofstream& outfile,const std::vector<Type>& to_write)
    {
      for(unsigned i = 0; i<to_write.size();i++){
        outfile<<to_write.at(i);
        if(i<to_write.size()-1)
          outfile<<",";
      }
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    template <class Type>
    Type read_arg_to_val(const std::string& arg)
    {
      Type to_return;
      std::stringstream stream(arg);
      if(!(stream >> to_return)) {
        std::cerr << "Invalid argument: "<<arg<<'\n';
      }
      return to_return;
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
  }
}
////////////////////////////////
#endif