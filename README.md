# SEGRED-NCO
Official Code for the Paper "Sample Efficient Geometry Reconstruction from Euclidean Distances using Non-Convex Optimization"

# Execution 
For executing our algorithm run the following command:
```
problem = struct;problem.type = 'GaussianData';problem.n = 500;N0=500;N0_inner=500;N0_firstorder=2000;tolerance = 1e-10;r_list=[3];oversampling_list=[2];instancesize=2;alg_name={'MatrixIRLS','ScaledSGD','AL_BurerMonteiro','ReiEDG'};snippet_EDG_phasetransitions
```

# Acknowledgements
Our code is adapted from [ScaledSGD](https://github.com/Hong-Ming/ScaledSGD), [ALM](https://github.com/abiy-tasissa/Nonconvex-Euclidean-Distance-Geometry-Problem-Solver) and [MatrixIRLS](https://github.com/ckuemmerle/MatrixIRLS). We sincerely thank the authors for releasing their code.
