### Artifacts for K-BEKEM$q_s$ ###

This repository contains the MATLAB implementation for the K-BEKEM$q_s$ framework.

1. Getting Started
The main execution point of the project is 'main.m'.
* To run the K-BEKEM$q_s$ algorithm, execute 'main.m' in the MATLAB environment.

2. Configuration and Environment Setup
To modify LWE parameters or the simulation environment, please open and edit 'dosetup.m'.
* LWE Parameters: You can configure $q, n, k, \sigma$ for generating LWE samples.
* Algorithm Parameters: Settings for $q_s, \alpha$, and other specific variables can be adjusted here.

3. Core Source Files
'main.m' requires the following dependency files:
* dosetup.m: Environment and parameter configuration.
* sample_from_LWE.m: LWE sample generation.
* allow_pm.m: plus-minus modulos handling. Ex) allow_pm(v=7, q=9) returns -2 
* getH.m: Finding H matrix 
* minv.m: modular inverse operation for scalar.
* find_primitive_equation.m: Finds primitive polynomials/equations.
* create_knapsack_map.m, trace_knapsack_items.m: Knapsack-based solving modules.
* Gaussian Elimination Modules:
    * `modular_GausElimination.m`: Returns results in the form $[I, D_1, b_1]$.
    * `GaussElimModulo.m`: Returns results in the form $[D_2, I, b_2]$.

4. Reproducing Figures (Visualization)
To reproduce the figures presented in the paper, use 'visualization_of_results.m'. This script can be run independently or as part of the main pipeline.

### How to use for the Visualization ###
Set the 'experimentType' variable within the script to 1, 2, or 3 to run specific simulations:

*  experimentType = 1: Runs visualize_ESS_and_RDD(). Compares the Enumeration Error Search (EES) algorithm with RDD. 
    This requires: compareEES_RDD.m, rdd_error_search.m, enumeration_error_search.m.
*  experimentType = 2: Runs compare_TopK_with_varying_topK_and_qs().
*  experimentType = 3: Runs compare_betamax_with_q().
    This requires calculate_betamax_with_varying_qs.m.
    In addition, experimentType 3 require symbolic execution for large q. so, we add extra functions which is based on symbolic execution: 
     *** getH_symbolic.m, new_modInv.m, modGaussianElimination_symbolic.m, find_primitive_equation_symbolic.m  

### Dependencies for Visualization:   ###
The following files must be in the same directory:
* inv_modulo.m, compareEES_RDD.m, rdd_error_search.m, enumeration_error_search.m, calculate_betamax_with_varying_qs.m.

5. Pre-processed Data  
The repository includes several ".mat" files containing pre-processed data to facilitate faster execution of the simulations and visualizations.
