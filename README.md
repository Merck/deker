# DEKER feature selection
Copyright (c) 2021 Merck Sharp & Dohme Corp. a subsidiary of Merck & Co., Inc., Kenilworth, NJ, USA.

by Sean M.S. Hayes

## Introduction
This library is made to perform feature selection based on a method originally proposed in by Sun et al. [1]. This library specifically relates to the methodology described in [2], named DEKER for **de**composed **ke**rnel **r**egression, which includes methods for identifying optimal hyperparameter values. This library was also designed for use in the context of network inference, also described in [2], by iteratively reapplying the DEKER method for feature selection accross all features of a dataset. 

The library may be used either as a header-only C++ library using the files in (include/deker) or compiled into a simple executable providing a command-line interface to the library with some simple I/O operations for reading in data and extracting output.

## Executable usage (src/deker.cpp):
The executable is built with the expectation that many iterations of feature selection will be done on a single data set, as in the context of network inference, most likely in a high performance computing (HPC) environment. The expected workflow is: 
1. Convert data to a binary format that can be quickly read accross all iterations of feature selection.
2. Check the number of responses to perform feature selection on to determine how many batch jobs to provision.
3. Execute each iteration of feature selection on a given response as a separate batch job.
4. Write output (contained in binary files) to a single 'edge-list'.
#### Convert input data (.csv) to binary file:
``` 
deker -b [input file location] [location to write output file]
```
#### Check number of responses:
```
deker -s [control file location]
```
* returns the number of responses specified by include_response/exclude_response to console
#### Do feature selection:
```
deker -w [location to write output file] -f [index of feature to use as response] [control file location]
```
* Program writes results to a binary file, use edgelist function to write binary file results to .csv
* Program expects control file to be plain text with each line in the format of [parameter] [value]
  * Control file must contain a line "data_file [data file location]". The data file should be a binary file produced by converting a .csv using the -b option, above. 
  * Optional parameters to include in control file: 
    * include_response: which features in the dataset allowed to be responses; expected value integer or integers seperated by comma (no spaces), defaults to all features
    * exclude_response: which features in the dataset not allowed to be responses; expected value integer or integers seperated by comma (no spaces), defaults to no features
    * include_predictor: which features in the dataset allowed to be predictors; expected value integer or integers seperated by comma (no spaces), defaults to all features
    * exclude_predictor: which features in the dataset not allowed to be predictors; expected value integer or integers seperated by comma (no spaces), defaults to no features
    * h_huber_loss: size of the quadratically smoothed elbow in the loss function; expected value double, defaults to .001
    * w_convergence_threshold: threshold of difference in feature weights (sqrt(w)) determining convergence; expected value double, defaults to .0000001
    * w_maxiter: maximum number of fixed-point iterations in solving for feature weights; expected value integer, defaults to 100
    * feature_drop_threshold: threshold of feature weight (sqrt(w)) to drop feature from calculation; expected value double, defaults to .001
    * lambda_init_max: maximum allowed lambda value; expected value double, defaults to .25
    * lambda_init_min: minimum allowed lambda value; expected value double, defaults to .005
    * lambda_init_stepsize: stepsize (fraction of maximum to decrease) when searching for 'good' lambda range (highest lambda without eliminating all features in [lambda_init_min,lambda_init_max]); expected value double, defaults to .5
    * lambda_init_maxiter: maximum iterations when searching for 'good' lambda range; expected value integer, defaults to 100
    * lambda_bo_maxiter: maximum iterations of Bayesian optimization for lambda within 'good' lambda range; expected value integer, defaults to 20
    
#### Write results to edge list:
```
deker -e [deker output file location] [location to write/append to edge list file]
```
## Header-only library usage (include/deker):

The library is primarily built around the 'deker::fit::Solve_deker' class defined in include/deker/fit.h, with additional functions providing support or I/O functions related to this class. Each instance of the class is constructed using the source data on which the feature selection model will be fit.
```
Solve_deker(const Eigen::Ref<const Eigen::MatrixXd> data_matrix,
                  const std::vector<unsigned>& predictors_use_index,
                  const unsigned& which_response,
                  const Opt_Param& control)
```
For the inputs, 'data_matrix' is expected to have samples in rows and features in columns, 'predictors_use_index' should indicate which columns are to be used as predictors, 'which_response' should indicate which column is to be used as response, and 'control' is a struct 'Opt_Param' defined in 'deker::fit::opt_param.h'. The class contains several methods related to solving to the heavily nested feature selection & hyperparameter problems.

#### Solve_deker method inputs
```
const Eigen::Ref<const Eigen::VectorXd> v
```
* An Eigen vector of double values, indicating weights on feature columns. More precisely, v are the square roots of the weights w, such that v^2 = w. When specified as a 'const' value in a given optimization step, weights are fixed and will not be updated.  
```
const std::vector<unsigned>& feature_index
```
* A vector of unsigned integers, which is of the same length as 'v' with each element indicating which column in Solve_deker's internal store of predictor data ('input_data.predictors') each element of 'v' corresponds to. When specified as a 'const' value in a given optimization step, weights are fixed and will not be updated in the optimization step. 
```
const double& lambda_regularization
const double& sigma_kernel_width
```
* Hyperparameter values; when specified as a 'const' value in a given optimization step, hyperparameter is fixed and will not be updated.
```
Output_deker& sol
```
* Struct containing fields for output of Solve_deker. Values provided in the struct may be used and will be overwritten and updated with the results of the optimization step. Full definition of the struct is in include/deker/fit/output.h. 

#### Solve_deker methods
```
std::tuple<double,double,double,double> calc_model_fit(const double& sigma_kernel_width,
                                                       const Eigen::Ref<const Eigen::VectorXd> v,
                                                       const std::vector<unsigned>& feature_index)
```
* 'calc_model_fit' determines the fit of the model given a value of 'sigma_kernel_width' (characteristic lengthscale of the kernel function), and feature weights ('v' and 'feature_index'). 'calc_model_fit' returns a tuple of four doubles with different quantities regarding model fit: (in order) estimated degrees of freedom, log-likelihood of fit, parameters A and B of the Platt scaling fit to condition classifier output to probabilities.
```
double opt_sigma(const Eigen::Ref<const Eigen::VectorXd> initial_v,
                       const std::vector<unsigned>& feature_index)
```
* 'opt_sigma' determines the optimal 'sigma_kernel_width' (characteristic lengthscale of the kernel function) value for a given set of feature weights ('v' and 'feature_index'). 'opt_sigma' returns a single double value, the optimal 'sigma_kernel_width' value. 
```
Eigen::VectorXd opt_v(const double& lambda_regularization,
                      const double& sigma_kernel_width,
                      const Eigen::Ref<const Eigen::VectorXd> v, 
                      const std::vector<unsigned>& feature_index)
```
* 'opt_v' determines the optimal 'v' (feature weight) values for given hyperparameter values ('lambda_regularization' and 'sigma_kernel_width'). 'opt_v' returns an Eigen vector, the new 'v' weights. 
```
void steff_iter_v(const double& lambda_regularization,
                  Output_deker& sol)
```
* 'steff_iter_v' performs one iteration of Steffensen's method for fixed point iteration [3] of the combined hyperparameter ('sigma_kernel_width') and feature weight ('v' & 'feature_index') optimization problem using 'opt_sigma' and 'opt_v' for a given value of the 'lambda_regularization' hyperparameter. 'steff_iter_v' returns values by updating the Output_deker struct 'sol.'
```
void inner_full_opt(const double& lambda_regularization,
                    Output_deker& sol)
```
* 'inner_full_opt' performs repeated fixed-point iteration until convergence of the combined hyperparameter ('sigma_kernel_width') and feature weight ('v' & 'feature_index') optimization problem using 'steff_iter_v' for a given value of the 'lambda_regularization' hyperparameter. 'inner_full_opt' returns values by updating the Output_deker struct 'sol.'
```
void opt_lambda(Output_deker& sol)
```
* 'opt_lambda' performs a complete optimization of both hyperparameters ('lambda_regularization' and 'sigma_kernel_width') and feature weights ('v' & 'feature_index') using Bayesian optimization [4] of the optimized model fit by 'inner_full_opt' for given values of 'lambda_regularization'. 'opt_lambda' returns values by updating the Output_deker struct 'sol.'

## References
1. Sun, Y.; Yao, J.; Goodison, S.: Feature selection for nonlinear regression and its application to cancer research. In: Proceedings of the 2015 SIAM International Conference on Data Mining, pp. 73-81. SIAM (2015)
2. Hayes, S.M.; Sachs, J.R.; Cho, C.R.: From complex data to biological insight: 'DEKER' feature selection and network inference. In Review. 
3. Johnson, L.W.; Scholz, D.R.: On Steffensen's method. SIAM Journal on Numerical Analysis. 5 (2): 296–302. (1968) doi:10.1137/0705026
4. Cully, A., Chatzilygeroudis, K., Allocati, F., Mouret, J.B.: Limbo: A flexible high-performance library for gaussian processes modeling and data-ecient optimization. Journal of Open Source Software 3(26) (2018)
