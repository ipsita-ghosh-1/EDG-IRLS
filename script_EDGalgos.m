%% Create/load data
problem = struct;
problem.type = 'GaussianData';
% problem.n = 500;
% problem.type = 'GaussianDataIllcond';
% problem.n = 400;
% problem.r = 5;
% problem.modeX0 = 'condition_control_1/x2';%'condition_control_log';%
% problem.cond_nr = sqrt(1e5);
% oversampling = 2.0;
rng('shuffle');
prob = dataloader_EDG(problem,oversampling);



%% Set algorithmic options
lambda = 0;
opts_custom = struct;
opts_custom.symmetricflag = true;
%%% Ground truth rank (potentially not known to algorithm)
opts_custom.r = prob.r ;
%%% Rank estimate used by algorithms
opts_custom.rtilde = prob.r ; % rank estimate used by algorithms
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
opts_custom.tol_CG = 1e-10;
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
%% Run differnet algorithms
alg_name = {'MatrixIRLS','ScaledSGD','AL_BurerMonteiro','ReiEDG'};%,'AL_BurerMonteiro'};%,'ReiEDG'};
% alg_name = {'MatrixIRLS'}
%alg_name = {'MatrixIRLS','AL_BurerMonteiro'}; % <- add here other algorithm to be used as further cells in cell array
% [Xr,outs,alg_name] = run_EDG_algos(prob,alg_name,opts_custom);




[outs_cell,success_matrix_procrustes, success_matrix,success_matrix_rankr,rel_error_PointDist,rel_error_Dist,...
    rel_error_DistObserved,times,prob_cell]= ...
    Phase_transition_parallel(problem,opts_custom,r_list,oversampling_list,...
    instancesize,alg_name);
%% Save results
curdate = datetime;
date = string(curdate);
filename = date+'_alg_name_'+'_type_'+[problem.type]+'_ovrsmp_'+...
[num2str(min(oversampling_list)),'_',...
num2str(max(oversampling_list))];
filename = filename+'_rmax_'+num2str(max(r_list))+'.mat';
save(filename,'success_matrix_procrustes','success_matrix','success_matrix_rankr','rel_error_PointDist','rel_error_Dist',...
    'rel_error_DistObserved',['opts_cu' ...
    'stom'],'problem',...
    'oversampling_list','r_list','alg_name','times','prob_cell','outs_cell');


