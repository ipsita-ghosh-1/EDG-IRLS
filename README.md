# 📌 Sample-Efficient Geometry Reconstruction from Euclidean Distances using Non-Convex Optimization (NeurIPS 2024)  

📄 **[Paper Link.](https://openreview.net/pdf?id=Yu7H8ZOuI2)**  

Welcome to the official repository of **Sample-Efficient Geometry Reconstruction from Euclidean Distances using Non-Convex Optimization**, presented at **NeurIPS 2024**! 

> **<p align="justify"> Abstract:** *The problem of finding suitable point embedding or geometric configurations given only Euclidean distance 
 information of point pairs arises both as a core task and as a sub-problem in a variety of machine learning applications. In this paper, we aim
to solve this problem given a minimal number of distance samples. To this end, we leverage continuous and non-convex rank minimization formulations of the problem and establish a local convergence guarantee for a variant of iteratively reweighted least squares (IRLS), which applies if a minimal random set of observed distances is provided. As a technical tool, we establish a restricted isometry property (RIP)
restricted to a tangent space of the manifold of symmetric rank-r matrices given random Euclidean distance measurements, which might be of independent interest for the analysis of other non-convex approaches. Furthermore, we assess data efficiency, scalability and generalizability of different reconstruction algorithms through numerical experiments with simulated data as well as real-world data, demonstrating the proposed algorithm’s ability to identify the underlying geometry from fewer distance samples compared to the state-of-the-art.*

# Execution 
For executing our algorithm run the following command: <br>

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
