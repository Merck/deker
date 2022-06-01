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
#ifndef DEKER_IO_DIGEST
#define DEKER_IO_DIGEST
////////////////////////////////
#include <deker/io/control.h>
////////////////////////////////
namespace deker {
    namespace io {
        struct Digest_Record{
            unsigned which_response;
            char solve_file[11];
            char lgbm_file[11];
            char opt_file[11];
            char run_file[11];
            char output_file[11];
            ///////////////////////////////
            Digest_Record(){
                solve_file[0]='\0';
                lgbm_file[0]='\0';
                opt_file[0]='\0';
                run_file[0]='\0';
                output_file[0]='\0';
            }
            ///////////////////////////////
            void serialize(std::ifstream& infile)
            {
                infile.read((char*) (&which_response), sizeof(which_response));
                infile.read((char*) (&solve_file), sizeof(solve_file));
                infile.read((char*) (&lgbm_file), sizeof(lgbm_file));
                infile.read((char*) (&opt_file), sizeof(opt_file));
                infile.read((char*) (&run_file), sizeof(run_file));
                infile.read((char*) (&output_file), sizeof(output_file));
            }
            ///////////////////////////////
            void serialize(std::ofstream& outfile) const
            {
                outfile.write((char*) (&which_response), sizeof(which_response));
                outfile.write((char*) (&solve_file), sizeof(solve_file));
                outfile.write((char*) (&lgbm_file), sizeof(lgbm_file));
                outfile.write((char*) (&opt_file), sizeof(opt_file));
                outfile.write((char*) (&run_file), sizeof(run_file));
                outfile.write((char*) (&output_file), sizeof(output_file));
            }
            ///////////////////////////////
            size_t serial_size() const
            {
                return sizeof(which_response)+sizeof(solve_file)+sizeof(lgbm_file)+sizeof(opt_file)+sizeof(run_file)+sizeof(output_file);
            }
        };
        class Control_Digest{
        protected:
            size_t full_record_size;
            Control_Base control;
            std::string digest_dir;
            std::vector<Digest_Record> record_vec;
        public:
            ///////////////////////////////
            Control_Digest(){
                full_record_size = 0;
            }
            Control_Digest(const std::string& digest_dir,
                           const Control_Base& control):
            digest_dir(digest_dir),control(control){
                full_record_size = 0;
            }
            ///////////////////////////////
            void add_record(unsigned which_response)
            {
                Digest_Record new_record;
                int which_file = record_vec.size();
                if(which_file>999999){
                    throw std::invalid_argument("Too many files for digest; maximum is 999,999\n");
                }
                char basename[7];
                for(int i=6;i>0;i--){
                    int new_digit = ((int)which_file/(int)pow(10,i-1))%10;
                    char new_char = '0'+new_digit;
                    int new_pos = 6-i;
                    basename[new_pos]=new_char;
                }
                strcat(new_record.solve_file,basename);
                strcat(new_record.solve_file,".sol");
                strcat(new_record.lgbm_file,basename);
                strcat(new_record.lgbm_file,".lgb");
                strcat(new_record.opt_file,basename);
                strcat(new_record.opt_file,".opt");
                strcat(new_record.run_file,basename);
                strcat(new_record.run_file,".run");
                strcat(new_record.output_file,basename);
                strcat(new_record.output_file,".out");
                new_record.which_response = which_response;
                record_vec.push_back(new_record);
                full_record_size++;
            }
            ///////////////////////////////
            void serialize(std::ifstream& infile)
            {
                infile.read((char*) (&full_record_size),sizeof(full_record_size));
                size_t size;
                infile.read((char*) (&size),sizeof(size));
                //
                control.serialize(infile);
                serialize_string(infile,digest_dir);
                //
                record_vec.clear();
                for(size_t i=0;i<size;i++){
                    record_vec.push_back(Digest_Record());
                    record_vec.back().serialize(infile);
                }
            }
            ///////////////////////////////
            void serialize(std::ofstream& outfile) const
            {
                //check to see that the files we're pointing to aren't written yet
                for(unsigned i = 0;i<record_vec.size();i++){
                    if(file_exists_test(digest_dir+std::string(record_vec.at(i).solve_file))||
                       file_exists_test(digest_dir+std::string(record_vec.at(i).lgbm_file))||
                       file_exists_test(digest_dir+std::string(record_vec.at(i).opt_file))||
                       file_exists_test(digest_dir+std::string(record_vec.at(i).run_file))){
                        throw std::invalid_argument("Digest-associated files found in digest directory - remove files or write to a clean directory\n");
                    }
                }
                //
                outfile.write((char*) (&full_record_size),sizeof(full_record_size));
                size_t size = record_vec.size();
                outfile.write((char*) (&size),sizeof(size));
                //write data
                control.serialize(outfile);
                serialize_string(outfile,digest_dir);
                //
                for(size_t i=0;i<size;i++)
                    record_vec.at(i).serialize(outfile);
            }
            ///////////////////////////////
            void read_single_record(const std::string& input_file, unsigned which_response)
            {
                if(!file_exists_test(input_file)||input_file.empty()){
                    throw std::invalid_argument("input file "+input_file+" not found");
                }
                std::ifstream infile(input_file, std::ios::in | std::ios::binary);
                //
                infile.read((char*) (&full_record_size),sizeof(full_record_size));
                size_t size;
                infile.read((char*) (&size),sizeof(size));
                //
                control.serialize(infile);
                serialize_string(infile,digest_dir);
                //
                if(which_response>=size){
                    throw std::invalid_argument("which_response is greater than the size of the digest\n");
                }
                //
                record_vec.clear();
                record_vec.push_back(Digest_Record());
                infile.seekg(which_response*record_vec.back().serial_size(),infile.cur);
                record_vec.back().serialize(infile);
                ///
                infile.close();
            }
            ///////////////////////////////
            void read_digest_size(const std::string& input_file)
            {
                if(!file_exists_test(input_file)||input_file.empty()){
                    throw std::invalid_argument("input file "+input_file+" not found");
                }
                std::ifstream infile(input_file, std::ios::in | std::ios::binary);
                //
                infile.read((char*) (&full_record_size),sizeof(full_record_size));
                ///
                infile.close();
            }
            ///////////////////////////////
            const Control_Base& get_control() const{return control;}
            const std::string& get_digest_dir() const{return digest_dir;}
            const Digest_Record& get_digest_record(unsigned i) const{return record_vec.at(i);}
            size_t get_full_record_size() const{return full_record_size;}
            size_t get_record_size() const{return record_vec.size();}
        };
    }
}
#endif