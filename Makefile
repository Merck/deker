all: test build

test: src/deker_test.cpp src/test
	g++ -w -O3 -march=native -mno-fma -fopenmp -o deker_test src/deker_test.cpp /SFS/user/ctc/hayesse/deker_include_libs/LightGBM/lib_lightgbm.a /SFS/user/ctc/hayesse/deker_include_libs/cryptopp/libcryptopp.a /SFS/user/ctc/hayesse/googletest/build/lib/libgtest.a /SFS/user/ctc/hayesse/googletest/build/lib/libgmock.a -I/SFS/user/ctc/hayesse/googletest/googletest/include -I/SFS/user/ctc/hayesse/googletest/googlemock/include -pthread -I/SFS/product/eigen/3.2.2/include -I/SFS/product/nlopt/2.4.2/rh68_x86_64/include -I/SFS/user/ctc/hayesse/deker_include_libs/LightGBM/include -I/SFS/user/ctc/hayesse/deker_command_line/include -I/SFS/user/ctc/hayesse/deker_include_libs -I/SFS/product/boost/1.59.0/include -lboost_system -lboost_filesystem -L/SFS/product/nlopt/2.4.2/rh68_x86_64/lib -lnlopt -lm -lstdc++fs  
	./deker_test

debug: src/deker.cpp
	g++ -w -Og -g -ggdb -march=native -fopenmp -o deker src/deker.cpp /SFS/user/ctc/hayesse/deker_include_libs/LightGBM/lib_lightgbm.a /SFS/user/ctc/hayesse/deker_include_libs/cryptopp/libcryptopp.a -I/SFS/product/eigen/3.2.2/include -I/SFS/product/nlopt/2.4.2/rh68_x86_64/include -I/SFS/user/ctc/hayesse/deker_include_libs/LightGBM/include -I/SFS/user/ctc/hayesse/deker_command_line/include -I/SFS/user/ctc/hayesse/deker_include_libs -I/SFS/product/boost/1.59.0/include -lboost_system -L/SFS/product/nlopt/2.4.2/rh68_x86_64/lib -lnlopt -lm -lstdc++fs 

build: src/deker.cpp
	g++ -w -O3 -march=native -mno-fma -fopenmp -o deker src/deker.cpp /SFS/user/ctc/hayesse/deker_include_libs/LightGBM/lib_lightgbm.a /SFS/user/ctc/hayesse/deker_include_libs/cryptopp/libcryptopp.a -I/SFS/product/eigen/3.2.2/include -I/SFS/product/nlopt/2.4.2/rh68_x86_64/include -I/SFS/user/ctc/hayesse/deker_include_libs/LightGBM/include -I/SFS/user/ctc/hayesse/deker_command_line/include -I/SFS/user/ctc/hayesse/deker_include_libs -I/SFS/product/boost/1.59.0/include -lboost_system -L/SFS/product/nlopt/2.4.2/rh68_x86_64/lib -lnlopt -lm -lstdc++fs 
