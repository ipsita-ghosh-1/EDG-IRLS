rng('shuffle'); % choose random seed
%% Chose option struct 
opts = struct;
opts.symmetricflag = true;
%%% General parameters
opts.N0 = N0; % this is the max. number of outer iterations for IRLS
opts.N0_inner = N0_inner; % this is the max. number of inner iterations for IRLS (max. nr. of iterations of iterative linear system solver)
opts.N0_firstorder = N0_firstorder; % this is the max. number of iterations for first-order algorithms (which includes all other algorithms except from MatrixIRLS currently)
opts.verbose = 1; % for printing status updates
%%% IRLS specific option parameters
opts.tol_CG = 1e-14; % IRLS inner problem tolerance
opts.tol = tolerance; % tolerance parameter for all algorithms
% opts.efficiencymode = 'fast';%
% opts.mode_linsolve = 'rangespace';%
opts.efficiencymode = 'fast';%[];%
opts.mode_linsolve ='tangspace';%'tangspace';%'rangespace';%'tangspace_exact';%rangespace';'tangspace'
opts.tol_CG = 1e-10;
opts.decomposition = 'svd';
%%% Option parameters for AugmentedLagrangian
opts.maxit = N0_firstorder;
%
%%% Option parameters for ScaledSGD
opts.learning_rate = 0.2; % based on the choice in the paper [Zhang et al. 2022] for ScaledSGD for EDM.
%% Create problem instances and run algorithms
[outs_cell,success_matrix_procrustes, success_matrix,success_matrix_rankr,rel_error_PointDist,rel_error_Dist,...
    rel_error_DistObserved,times,prob_cell]= ...
    Phase_transition_parallel(problem,opts,r_list,oversampling_list,...
    instancesize,alg_name);
%% Save results
curdate = datetime;
date = string(curdate);
filename = date+'_alg_name_'+'_type_'+[problem.type]+'_ovrsmp_'+...
[num2str(min(oversampling_list)),'_',...
num2str(max(oversampling_list))];
filename = filename+'_rmax_'+num2str(max(r_list))+'.mat';
save(filename,'success_matrix_procrustes','success_matrix','success_matrix_rankr','rel_error_PointDist','rel_error_Dist',...
    'rel_error_DistObserved','opts','problem',...
    'oversampling_list','r_list','alg_name','times','prob_cell','outs_cell');
