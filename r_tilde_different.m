%% Create/load data
rng('shuffle');


%% Set algorithmic options
lambda = 0;
opts_custom = struct;
opts_custom.symmetricflag = true;
%%% Ground truth rank (potentially not known to algorithm)
opts_custom.r = problem.r ;
%%% Rank estimate used by algorithms
opts_custom.rtilde = 2*problem.r  ; % rank estimate used by algorithms
%%% IRLS options
N0 = 200; 
N0_inner = problem.n;
N0_firstorder = 2000;
tolerance = 1e-14;
opts_custom.verbose = 1;
opts_custom.N0 = N0;
opts_custom.N0_inner = N0_inner;
opts_custom.N0_firstorder = N0_firstorder;
opts_custom.efficiencymode = 'fast';%[];%
opts_custom.mode_linsolve ='tangspace';%'tangspace';%'rangespace';%'tangspace_exact';%rangespace';'tangspace'
opts_custom.tol_CG = 1e-6;
opts_custom.decomposition = 'svd';
opts_custom.tol = tolerance;
% options for BB gradient method with nonmonotone line search
% optional for line search algorithm
opts_custom.xtol = tolerance;
opts_custom.gtol = tolerance; 
opts_custom.ftol = tolerance; 
opts_custom.alph  = 1e-3; 
opts_custom.rho  = 1e-4; 
opts_custom.sigma  = 0.1; 
opts_custom.eta  = 0.8; 
%%
opts_custom.learning_rate = 0.2; % based on the choice in the paper [Zhang et al. 2022] for ScaledSGD for EDM.
%% Create problem instances and run algorithms
[outs_cell,success_matrix_procrustes, success_matrix,success_matrix_rankr,rel_error_PointDist,rel_error_Dist,...
    rel_error_DistObserved,times,prob_cell]= ...
    Phase_transition_parallel(problem,opts_custom,r_list,oversampling_list,...
    instancesize,alg_name);
%% Save results
curdate = datetime;
date = string(curdate);

filename = 'tilde2_'+date+'_alg_name_'+'_type_'+[problem.type]+'_ovrsmp_'+...
[num2str(min(oversampling_list)),'_',...
num2str(max(oversampling_list))];
filename = filename+'_rmax_'+num2str(max(r_list))+'.mat';
save(filename,'success_matrix_procrustes','success_matrix','success_matrix_rankr','rel_error_PointDist','rel_error_Dist',...
    'rel_error_DistObserved','opts_custom','problem',...
    'oversampling_list','r_list','alg_name','times','prob_cell','outs_cell');