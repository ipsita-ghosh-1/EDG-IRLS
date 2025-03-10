function opts = getDefaultOpts_IRLS

opts.p              = 0;  % Schatten-p parameter (p close to 0 for very non-convex objective, p=1 for convex)
opts.N0             = 200;  % max. number of outer iterations (second order method like IRLS)
opts.N0_inner       = 50;  % max. number of inner iterations (second order method like IRLS)
opts.N_SVD          = 20;  % max. number of iterations for power method-type solver for partial SVD (such as bksvd)
opts.tol            = 1e-9;  % stopping criterion, lower bound on relative change of Frobenius norm
opts.tol_CG         = 1e-5;
opts.epsmin         = 1e-15;  % minimal value for epsilon smoothing
opts.use_min        = 1;
opts.type_mean      = 'optimal';
opts.mode_linsolve  = 'tangspace';
opts.preconditioning= 'wo_precond';
opts.qmean_para     = min(opts.p/(opts.p-2),0);
opts.mode_eps       = 'oracle_model_order';% iter_diff'
opts.tracking       = 0;
opts.adaptive_cg    = 0;
opts.verbose        = 1;
opts.decomposition  = 'svd';%'svd';
opts.use_Hermitian  = false;
opts.symmetricflag  = false;
opts.lambda         = 0;
end

    