# SEGRED-NCO
Official Code for the Paper "Sample Efficient Geometry Reconstruction from Euclidean Distances using Non-Convex Optimization"

# Execution 
For executing our algorithm run the following command:
Protein Data 
```
problem = struct;problem.type = '1AX8';N0=500;N0_inner=500;N0_firstorder=5000;tolerance = 1e-12;r_list=[2,3];oversampling_list=[0.5,1.1];instancesize=1;alg_name={'MatrixIRLS','ScaledSGD','AL_BurerMonteiro','ReiEDG'};snippet_EDG_phasetransitions;
```
Gaussian Setup 
```
problem = struct;problem.type = 'GaussianData';problem.n = 50;problem.r = 2;N0=500;N0_inner=500;N0_firstorder=5000;tolerance = 1e-12;r_list=[2,3];oversampling_list=[0.5,1.1];instancesize=1;alg_name={'MatrixIRLS','ScaledSGD','AL_BurerMonteiro','ReiEDG'};snippet_EDG_phasetransitions;
```
Ill-conditioned Setup
```
problem = struct;problem.type ='GaussianDataIllcond';problem.n = 500;problem.modeX0 = 'condition_control_1/x2';problem.cond_nr =sqrt(1e5);N0=500;N0_inner=5000;N0_firstorder=2000;tolerance =1e-10;r_list=[2,3,4,5];oversampling_list=[1,1.5,2,2.5,3,3.5,4];instancesize=8;alg_name={'MatrixIRLS','ScaledSGD','AL_BurerMonteiro','ReiEDG'};snippet_EDG_phasetransitions
```
USCities Data:
```
problem = struct;problem.type = 'USCities';N0=500;N0_inner=500;N0_firstorder=5000;tolerance = 1e-12;r_list=[2,3];oversampling_list=[0.5,1.1];instancesize=1;alg_name={'MatrixIRLS','ScaledSGD','AL_BurerMonteiro','ReiEDG'};snippet_EDG_phasetransitions;
```

The output file contains the results from our algorithm along with the results of other comparison algorithms, including ScaledSGD, RieEDG, and ALM, as referenced below.

# Acknowledgements
Our code is adapted from [ScaledSGD](https://github.com/Hong-Ming/ScaledSGD), [ALM](https://github.com/abiy-tasissa/Nonconvex-Euclidean-Distance-Geometry-Problem-Solver) and [MatrixIRLS](https://github.com/ckuemmerle/MatrixIRLS). We sincerely thank the authors for releasing their code.
