function [Xr,outs] = ...
    MatrixIRLS(prob,lambda,opts)
% Given partial observations of a matrix, attempts to find a low-rank completion.
% =========================================================================
% Parameters
% ----------
% prob: struct. Contains problem parameters.
% lambda: double. Data fit parameter, choose as 0 for equality data
%           constraint, and larger value for noisy/inaccurate observations. 
%           Large lambda corresponds to a loose data fit, see (1) above.
% opts: struct. Contains algorithmic options. Fields include:
%       - p         double, between 0 and 1. Parameter of (non-)convexity of
%                   objective function.
%                   = 0: sum of log objective /log-det objective.
%                   > 0 and < 1: Schatten-p quasi-norm objective.
%                   = 1: Nuclear norm type objective.
%       - N0        int. Max. number of (outer) IRLS iterations 
%                   (default 200). 
%       - N0_inner  int. Max. number of inner iteration in the conjugate
%                   gradient method (default 200).
%       - tol       double. Stopping criterion, lower bound for relative 
%                   change of Frobenius norm (default 1e-9).
%       - tol_CG    double. Determines the tolerance for stopping criterion
%                   of the conjugate gradient method (bound on the relative
%                   residual norm). (default 1e-5).
%       - epsmin    Minimal value for regularization parameter epsilon
%                   (default 1e-15).
%       - mode_eps   character string. Choice of rule to update the
%                    regularization parameter epsilon (default
%                    'oracle_model_order').
%                   = 'oracle_model_order': Uses knowledge about the
%                      model order/rank "r", epsilon choice as in [1-3].
%                   = 'auto_decay': eps that is automatically
%                      superlinearly decaying. First value eps0 of epsilon 
%                      is based on the Frobenius norm of first iterate
%                      X^{(l)}, choose eps at iteration l such that:
%                      eps(l) = eps0.*(tol_outer/eps0)^(((l-1)/(N0-1))^(2-p));
%                   = 'iter_diff': epsilon choice rule similar to [4], 
%                      based on the norm of the difference of two conse-
%                      cutive iterates.
%       - use_min   Boolean (default 1).
%                   = 1: Force monotonicity of regularization parameter 
%                   epsilon.
%                   = 0: No forced monotonicity, plain application of
%                   epsilon choice rule.
%       - rtilde    int. Rank estimate
%       - R         int. Maximum rank parameter for intermediate iterates.
%                   (default: min(d1,d2)). Can be set to smaller values 
%                   (upper bound for the rank of the solution), used to 
%                   accelerate algorithm for large problems if target rank 
%                   is misspecified or if a non-oracle choice of the 
%                   smoothing parameter epsilon is chosen.
%       - type_mean   character string, determines type of weight operator
%                   underlying the IRLS method, see [3].
%                   = 'optimal': Corresponds to the optimal weight
%                      operator: Geometric mean of p = 0, and
%                      (p/(p-2))-power mean for Schatten-p objective.
%                   = 'geometric': Geometric mean.
%                   = 'harmonic': Harmonic mean.
%                   = 'arithmetic: Arithmetic mean.
%                   = 'min': Based on minimum of right- and
%                      left-weights.
%                   = 'qmean': q-power mean.
%       - qmean_para  double. Only relevant if type_mean=='qmean', 
%                      determines the value of "q" in the q-power mean,
%                      cf. [3].
%       - 0
%                   Boolean (default = 0). Can be set to 1, but that is
%                   only interesting for theoretical reasons.
%                   = 1: In the definition of the weight operator, this 
%                      means that the weights on the "antisymmetric part" 
%                      are increased as in [3]. 
%                   = 0: Weight operator without artificially increased
%                      antisymmettic port.
%       - objective character string, decides precise form of objective
%                       function to be optimized by IRLS.
%                   = 'objective_thesis': log(max(sigma,eps))+0.5.*
%                       (min(sigma,eps)^2/eps^2-1).
%                       (as used in [1-3])
%                   = 'pluseps': Corresponds to objective
%                       of the type:
%                       log(max(sigma,eps)+eps)+0.25.*(min(sigma,eps)^2/eps^2-1)
%                   = 'pluseps_squared_max': Corresponds to objectives
%                       of the type:
%                       log(max(sigma,eps)^2+eps^2)+0.5.*(min(sigma,eps)^2/eps^2-1)
%                   = 'pluseps_squared': Corresponds to objectives
%                       of the type:
%                       log(sigma^2+eps^2)  
%       - verbose   int. Detemines the level of displayed text and
%                   amount of output variables (default 1)
%                   = 0: Very little text output, only output basic
%                        parameters such as regularization parameters
%                        eps_c.
%                   = 1: Return also intermediate iterates and a few 
%                        algorithmic parameters, more text display.
%                   = 2: Return also many further parameters.
%       - saveiterates  Boolean. Save or do not save intermediate iterates
%                   (default 1).
%                   = 0: Only output last iterate.
%                   = 1: Output all intermediate iterates.
%       - mode_linsolve  character string, describes linear system to be solved.
%                   = 'tangspace': Linear system based on tangent space,
%                   see Appendix of [1-2] or Section 2.6.2. of [3].
%                      Solve linear system by
%                      conjugate gradients, using fast operators and a
%                      using a low-rank structure of weight operator.
%                      Works with the full space approach, i.e., it 
%                      solves the linear system.
%                      (lambda.*W+Phi'*Phi)*x=Phi'(y) for x.
%                   = 'rangespace': Solves linear system of the range space
%                      iteratively by CG, i.e., the system
%                      (lambda*Id+Phi*W^(-1)*Phi')*z=y for z (then
%                      the new iterate is x = W^(-1)*Phi'*z).
%       - tangent_para character string, decides what kind of computational 
%                   representation of the tangent space is used, cf. also
%                   [4]. Default: 'extrinsic'.
%                   = 'intrinsic': Uses an intrinsic representation of 
%                      elements of the tangent space (which uses 
%                      r_k*(d1+d2-r_k) parameters). This is used to obtain 
%                      to have a 1-to-1 correspondence with the implemen-
%                      tation detailed in Appendix A of [2], which is
%                      needed to obtain make the well-conditioning result
%                      of Theorem 4.2 of [2] applicable.
%                   = 'extrinsic': Uses an extrinsic representation of 
%                      elements of the tangent space (which uses 
%                      r_k*(d1+d2+r_k) parameters). Compared to the imple-
%                      mentation detailed in Appendix A of [2], this does 
%                      not use Step 6 in Algorithm 3 nor Step 1 of
%                      Algorithm 4 in [2]. This can be used as a default,
%                      as in practice, we typically observed little 
%                      performance gain compared to the intrinsic
%                      representation, and since the extrinsic one is 
%                      conceptually simpler.
%       - adaptive_cg  boolean. 
%                   = 1: Uses varying implentations of the least squares 
%                      solve, depending on the conditioning of
%                      repsective problems (not implemented yet).
%                   = 0: Choose always the linear system described by
%                      'mode_linsolve'.
% Returns
% ----------
% Xr:   struct. Contains the last iterate in space-efficient form, see
%       cofficient matrices in Appendix of [1]. In particular,
%                   Xr.Gam1 = \tilde{Gamma}_1  [(r_k x r_k) matrix]
%                   Xr.Gam2 = \tilde{Gamma}_2  [(r_k x d_2) matrix]
%                   Xr.Gam3 = \tilde{Gamma}_3  [(d_1 x r_k) matrix]
%                   Xr.res_range = r_{k+1}     [(m x 1) vector].
%                   Xr.U    = U^{(k)}          [(d_1 x r_k) matrix w/
%                                                   orthonormal columns]
%                   Xr.V    = V^{(k)}          [(d_2 x r_k) matrix w/
%                                                   orthonormal columns]
%                   Dense matrix X^{(N)} can be recovered from Xr
%                   using the function 'get_full_mat' such that
%                   X^{(N)} = U^{(k)}*(\tilde{Gamma}_1*V^{(k)'}+\tilde{Gamma}_2)
%                   + \tilde{Gamma}_3*V^{(k)'} + P_{\Omega}^*(r_{k+1}).
% outs: struct. Contains algorithmic information such as the smoothing
%                   parameters eps, level of detail depends on 'saveiterates' 
%                   and 'verbose'. Fields include:
%       - X         (1 x N) cell, with structs: As Xr, but for each
%                   iteration (only if opts.saveiterates == 1).
%       - N         int. Number of performed outer iterations. 
%       - eps       (1 x N) vector. Smoothing parameters '\epsilon_k' for 
%                   each iteration.
%       - r_greatereps (1 x N) vector. Current active ranks 'r_k' for each
%                   iteration, such that r_k = |\{i \in [d]:
%                   \sigma_i(X^{(k)}) > \epsilon_k\}|, cf. [1-3].
%       - time      (1 x N) vector. Cumulative timestamps at each
%                   iteration.
% =========================================================================% =========================================================================
%% Obtain algorithmic options

AdjInvVersion = [];
global PhiPhiT_inv_y_c

[N0,N0_inner,N_SVD,p,tol_CG,tol_outer,objective,...
    type_mean,qmean_para,mode_linsolve,adaptive_cg,mode_eps,epsmin,...
    use_min,tracking,verbose,...
    recsys,saveiterates,tangent_para,forcefull]=...
    get_options(opts);
if isfield(opts,'X0')
    X0_oracle = opts.X0;
else
    X0_oracle = [];
end
if adaptive_cg == 1
    mode_linsolve = 'rangespace';
    if verbose > 1
        mode_linsolve_it = cell(1,N0);
    end
end
if isfield(opts,'decomposition')
    decomposition = opts.decomposition;
end

%% Obtain problem information
[n,meas_op,y]=get_problem_parameters(prob,opts);
d=min(n.d1,n.d2);
D=max(n.d1,n.d2);
%% Obtain oracle information, if applicable
oracle_knowledge = 1;
r           = opts.rtilde;       % Rank to be used in the definition of 
                                 % smoothing parameter choice rule
%% Initialize logging information
eps             = zeros(1,N0);
r_greatereps    = zeros(1,N0);
time            = zeros(1,N0);
if saveiterates
    X           = cell(1,N0);
    Z           = cell(1,N0);
    sings       = cell(1,N0);
    U           = cell(1,N0);
    V           = cell(1,N0);
end
if verbose > 1
    N_inner = zeros(1,N0);
    resvec  = cell(1,N0);
    quot_r  = zeros(1,N0);
    cond_Xell = zeros(1,N0);
end 
% =========================================================================
if tracking  % used if objective value F_{\epsilon,k} is to be tracked
    obj_after_eps_update = zeros(1,N0);
    obj_before_eps_update = zeros(1,N0);
    rankobj_before_eps = zeros(1,N0);
    rankobj_after_eps = zeros(1,N0);
end
%% Initialization
[weight_op,X_c,eps_c,r_c,r_upper,first_iterate,X_c_newbase] =...
    MatrixIRLS_initialize_weights(n,opts,oracle_knowledge,r,y,meas_op,prob);
weight_op_previous = weight_op;
z_c = [];
%% Start iterations
if verbose > 0
    fprintf(['Run MatrixIRLS with p=%.2f, data fit parameter lambda=%.2e and smoothing parameter rule "%-11s"...\n'], ...
    p,lambda,mode_eps);
end
k=1;
tic
while k <= N0
    %% Calculate the solution of the minimization problem
    if mod(k,25) == 0 && verbose > 0
        disp(['Begin Iteration ',num2str(k),'...']);
    end
    eps_previous = eps_c;
    if first_iterate
        if weight_op.symmetricflag
            dim_tangent = r_c*((r_c+1)/2+n.d1);
        else
            dim_tangent = r_c*(n.d1+n.d2+r_c) ;
        end
        [X_c_previous] = prepare_iterate_output(zeros(dim_tangent,1),z_c,...
            y,meas_op,weight_op,1,lambda,...%zeros(length(y),1),meas_op,weight_op,1,lambda,...
            mode_linsolve,first_iterate);
        X_c_previous.U = weight_op.U;
        X_c_previous.V = weight_op.V;
    else

        X_c_previous   = X_c;
    end
    X_prev_newbase = X_c_newbase;
    %% Update the iterate: Determine X^{(k)}
    [X_c,X_Tn_c,z_c,N_inner_c,resvec_c,relres_c,cond_c]=update_iterate(weight_op,...
        meas_op,n,y,X_prev_newbase,z_c,...
        lambda,tol_outer,tol_CG,mode_linsolve,...
        first_iterate,N0_inner,...
        type_mean,qmean_para,tangent_para,AdjInvVersion,objective,...
        X_c_previous.res_range.');
    if strcmp(meas_op.problem_type,'MultidimensionalScaling')
        X_c.ytilde = y';
    end
    %% Calculate norm difference of consecutive iterates for stopping criterion (and epsilon rule)   
    diff_c  = compute_iterate_diff_norm(X_c,X_c_previous,meas_op,'proxy');
    norm_Xc = norm_frob_compact(X_c,meas_op,weight_op);
   rel_chg = diff_c/norm_Xc;
    if not(first_iterate)
        if rel_chg < tol_outer
            N0 = k;
        end
        if N_inner_c == 0
            N0 = k-1;
        end
    end
    weight_op_previous = weight_op;
    %% Define handle for obtaining spectral information
    forcefull = false;
    if strcmp(opts.decomposition,'svdexact')
        forcefull = true;
    end
    [X_c_handle,X_c_full] = get_matrix_handle(X_c,meas_op,weight_op,n,'forcefull',forcefull);
    % X_c_full(end,end)
    % tmp = kappainv_EDG(X_c_full);
    % norm(tmp(meas_op.Omega)-prob.Dist_sampled(meas_op.Omega))
    %% Obtain spectral information
    [weight_op.U, sing_c, weight_op.V,singmin,...
        weight_op.Uperp,weight_op.Vperp,weight_op.symTmat,weight_op.symInd] = ...
        get_weight_update(X_c_handle,d,r_c,r,oracle_knowledge,N_SVD,opts,...
    'decomposition',decomposition,'Uinit',weight_op.U,'Vinit',weight_op.V);
    weight_op.sing=diag(sing_c);
    %% Initalize current rank
    eps0 = weight_op.sing(1)/2;
    %norm(X_c_full-X_c_full')./norm(X_c_full)
%     if not(isempty(initialization))
%         eps0 = weight_op.sing(1)/2;
    if first_iterate
        weight_op.sing = weight_op.sing;
        first_iterate = 0;
        r_c = r+1;
    end
    %% Update smoothing parameter \epsilon_k
    if tracking
        X_c.U = weight_op_previous.U;
        X_c.V = weight_op_previous.V;
        X_c_full = get_densemat_from_iterate(X_c,meas_op);
        X_c.Fullmat = X_c_full;
        if k > 1
            X_c_full_old = X_c_full;
        else
            X_c_full_old = eye(n.d1,n.d2);
        end
        [~,singval_mat_c,~]=svd(X_c_full);
        sing_c_all = diag(singval_mat_c);
        eps_c = update_epsilon(eps_c,mode_eps,...
        oracle_knowledge,n,epsmin,use_min,p,sing_c_all,r,...
        diff_c,eps0,tol_outer,k,N0,singmin,X0_oracle,X_c_full);
    else
        if first_iterate
        else
            eps_c = update_epsilon(eps_c,mode_eps,...
            oracle_knowledge,n,epsmin,use_min,p,weight_op.sing,r,...
            diff_c,eps0,tol_outer,k,N0,singmin,X0_oracle);
        end

    end
    %% Update weight operator: W^{(k)}
    weight_op.eps = eps_c;
    if not(strcmp(objective,'pluseps_squared'))
        [r_c,weight_op,eps_c] = update_rankpara(X_c_handle,weight_op,...
            eps_c,r_upper,r_c,N_SVD,0,opts,'decomposition',decomposition);
    end
    if tracking
        weight_op.Uperp = null(weight_op.U(:,1:r)');
        weight_op.Vperp = null(weight_op.V(:,1:r)');
    else
        weight_op.Uperp = [];
        weight_op.Vperp = [];
    end
    
    %%% Calculate quotient quot_r_c between smallest singular value larger 
    %%% than \epsilon_k and \epsilon_k
    if r_c == 0
        if not(strcmp(mode_eps,'auto_decay'))
            N0 = k-1;
        end
        quot_r_c = 1;
    else
        quot_r_c = weight_op.sing(r_c)./eps_c;
    end    
    x_c_newbase_resi=proj_tangspace_from_Oprange(X_c.res_range.',meas_op,...
            weight_op,'tangspace');%% update this for 0
    % tmp2 = meas_op.kappa(meas_op.kappainv(meas_op.kappainv(reshape(meas_op.Samp_opmat'*X_c.res_range.',n.d1,n.d1))));
    % tmp2 = (-2).*meas_op.kappa(meas_op.kappainv(reshape(meas_op.Samp_opmat'*X_c.res_range.',n.d1,n.d1)));
    % tmp2 = applyMeasopBackward(X_c.res_range.',meas_op);
    % x_c_newbase_resi = matrixspace_to_tangspace(reshape(tmp2,n.d1,n.d1),weight_op.U,weight_op.V);

    x_c_newbase_Tm1 = proj_tangspace(X_Tn_c,'MatrixCompletion',...
        weight_op_previous,weight_op);
    X_c_newbase=x_c_newbase_resi+x_c_newbase_Tm1;
    [weight_op.Wl_inv,weight_op.Wl_inv_eps] = ...
        set_weight_infovec(weight_op.sing(1:r_c),eps_c,p,objective);
    if tracking
        %% Track value of objective function
            [~,singval_mat_c,~]=svd(X_c_full_old);
            sing_all = diag(singval_mat_c);
            [Wl_inv_forH,Wl_inv_eps_forH] = ...
            set_weight_infovec(sing_all,eps_c,p,objective);
            [H_full,~,~,dHm1_tracking] = weight_op_lowrank_prec(zeros(1,n.d1),zeros(1,n.d2),...
                [Wl_inv_forH;Wl_inv_eps_forH],1,1,...
            type_mean,0,[],[],qmean_para);
            H_tracking = cell(2,2);
            if not(isempty(r))
                compute_Htracking = 1;
                r_tracking = r;
            else
                r_tracking = 1;
%                 r_tracking = nnz(svd(opts.X0)>1e-13);
                compute_Htracking = 0;
            end
            if compute_Htracking
                H_tracking{1,1} = H_full{1}(1:r,1:r);
                H_tracking{1,2} = H_full{1}(1:r,r+1:end);
                H_tracking{2,1} = H_full{1}(r+1:end,1:r);
                H_tracking{2,2} = H_full{1}(r+1:end,r+1:end);
            end
            outs.H_tracking{k} = H_tracking;
            if isfield(opts,'X0')
%                 if not(isempty(r))
                if compute_Htracking
                    outs.H22_UperpX0Vperp{k} = H_tracking{2,2}.*(weight_op.Uperp'*opts.X0*weight_op.Vperp);
                    outs.H22_UperpX0Vperp_norm(k) = norm(H_tracking{2,2}.*(weight_op.Uperp'*opts.X0*weight_op.Vperp));
                    outs.H_X0_norms{k}(1,1) = norm(H_tracking{1,1}.*(weight_op.U(:,1:r)'*opts.X0*weight_op.V(:,1:r)));
                    outs.H_X0_norms{k}(2,1) = norm(H_tracking{2,1}.*(weight_op.Uperp'*opts.X0*weight_op.V(:,1:r)));
                    outs.H_X0_norms{k}(1,2) = norm(H_tracking{1,2}.*(weight_op.U(:,1:r)'*opts.X0*weight_op.Vperp));
                    outs.H_X0_norms{k}(2,2) = outs.H22_UperpX0Vperp_norm(k);
    %                 end
                    outs.UperpU0{k} = weight_op.Uperp'*opts.U0;
                    outs.VperpV0{k} = weight_op.Vperp'*opts.V0;
    %                 norm(outs.UperpU0{k})
                    outs.UperpU0_norm{k} = norm(outs.UperpU0{k});
                    outs.VperpV0_norm{k} = norm(outs.VperpV0{k});
                end
                outs.XminX0spectral(k) = norm(opts.X0-X_c_full);
                sXminX0 = svd(opts.X0-X_c_full);
                outs.XminX0nuclear(k) = sum(sXminX0);
                %norm(H_tracking{2,2}.*(weight_op.Uperp'*(X_c_full_old-X_c_full)*weight_op.Vperp))
%                 if k > 1
%                     disp(['Scalar product:',num2str(trace((H_tracking{2,2}.*(weight_op.Uperp'*(X_c_full_old-X_c_full)*weight_op.Vperp)'...
%                         *(weight_op.Uperp'*(X_c_full_old-X_c_full)*weight_op.Vperp))))]);
%                 end
                singsX0=svd(opts.X0);
                outs.closeness(k)= norm(opts.X0-X_c_full)./singsX0(r_tracking);
            end
        rankobj_before_eps(k)= eps_previous.^(2-p).*get_rankobjective(sing_all,eps_previous,p,objective);
        rankobj_after_eps(k) = eps_c.^(2-p).*get_rankobjective(sing_all,eps_c,p,objective);
        obj_before_eps_update(k) = rankobj_before_eps(k);
        obj_after_eps_update(k) = rankobj_after_eps(k);
    end
    %% Save output information
    if eps_c <= epsmin || quot_r_c > 1e20
        N0 = k;
    end
    eps(k)=eps_c;
    r_greatereps(k)=r_c;
    time(k) = toc;
    X_c.U = weight_op_previous.U;
    X_c.V = weight_op_previous.V;
    if saveiterates
        %% Save all iterates
        X{k}    = X_c;
        sings{k}= weight_op_previous.sing;
        %%%%%%%%%%%%
        U{k}    = weight_op_previous.U;
        V{k}    = weight_op_previous.V;
    end
    if verbose > 1
        N_inner(k) = N_inner_c;
        resvec{k}  = resvec_c;
        quot_r(k)  = quot_r_c;
        if r_c == 0
            cond_Xell(k) = 0;
        else
            cond_Xell(k) = weight_op.sing(1)./eps_c;
        end
        mode_linsolve_it{k} = mode_linsolve;
    end
    if adaptive_cg == 1
        %%% Update implementation mode to be used at next iteration, based
        %%% value of quot_r_c
        mode_linsolve=use_optimal_cg_mode(quot_r_c);
    end
    if strcmp(type_mean,'Hessianupper_1') || strcmp(type_mean,'Hessianupper_2')...
            || strcmp(type_mean,'Hessianupper_Ribeiro')
        %% Modify weight operator information (only for an older implementation)
        weight_op.Wl_inv = weight_op.sing(1:r_c);
        weight_op.Wl_inv_eps    = eps_c;
    end
    outs.cond(k) = cond_c;
    if verbose > 0
        %%% Print algorithmic information
        if k == 1
%                 fprintf(['%-11s k: %3d  N_inner: %5d   ', ...
%                 'relres: %.3e  eps_k: %.3e r_k: %d qt_rk: %.2f relchg: %.2e, cond: %.2f\n'], ...
%                 'init',k,N_inner_c,relres_c,eps_c,r_c,quot_r_c,rel_chg); 
            fprintf(['%-11s k: %3d  N_inner: %5d   ', ...
            'relres: %.3e  eps_k: %.3e r_k: %d qt_rk: %.2f relchg: %.2e, cond: %.2f\n'], ...
            'init',k,N_inner_c,relres_c,eps_c,r_c,quot_r_c,rel_chg,cond_c); 
        else
%                 fprintf(['%-11s k: %3d  N_inner: %5d   ', ...
%                 'relres: %.3e  eps_k: %.3e r_k: %d qt_rk: %.2f relchg: %.2e\n'], ...
%                 mode_linsolve,k,N_inner_c,relres_c,eps_c,r_c,quot_r_c,rel_chg);  
            fprintf(['%-11s k: %3d  N_inner: %5d   ', ...
            'relres: %.3e  eps_k: %.3e r_k: %d qt_rk: %.2f relchg: %.2e, cond: %.2f\n'], ...
            mode_linsolve,k,N_inner_c,relres_c,eps_c,r_c,quot_r_c,rel_chg,cond_c); 
        end
    end
    k=k+1;
end
%% Tidy up the size of the relevant output arrays and cells
Xr      = X_c;
outs.N      = N0;
outs.eps    = eps(1:N0);
if saveiterates
    outs.X      = X(1:N0);
    outs.sings  = sings(1,1:N0);
    outs.UU     = U(1:N0);
    outs.VV     = V(1:N0);
else
    outs.sings  = weight_op_previous.sing;
end
outs.r_greatereps   = r_greatereps(1:N0);
outs.time   = time(1:N0);
outs.opts   = opts;
outs.meas_op = meas_op;
if verbose > 0
    outs.Xr_handle = X_c_handle;
    if verbose > 1
        outs.N_inner  = N_inner(1:N0);
        outs.resvec   = resvec(1:N0);
        outs.quot_r   = quot_r(1:N0);
        outs.cond_Xell = cond_Xell(1:N0);
        outs.mode_linsolve_it = mode_linsolve_it(1:N0);
    end
end
if tracking
    outs.obj_epsA = obj_before_eps_update(1:N0);
    outs.obj_epsB = obj_after_eps_update(1:N0);
    outs.rankobj_before_eps = rankobj_before_eps(1:N0);
    outs.rankobj_after_eps = rankobj_after_eps(1:N0);
    outs.H_tracking =  outs.H_tracking(1:N0);
    if isfield(opts,'X0')
        if compute_Htracking
            outs.H22_UperpX0Vperp = outs.H22_UperpX0Vperp(1:N0);
            outs.H22_UperpX0Vperp_norm = outs.H22_UperpX0Vperp_norm(1:N0);
            outs.H_X0_norms = outs.H_X0_norms(1:N0);
            outs.UperpU0 = outs.UperpU0(1:N0);
            outs.VperpV0 = outs.VperpV0(1:N0);
            outs.UperpU0_norm = outs.UperpU0_norm(1:N0);
            outs.VperpV0_norm = outs.VperpV0_norm(1:N0);
        end
        outs.XminX0spectral = outs.XminX0spectral(1:N0);
        outs.XminX0nuclear = outs.XminX0nuclear(1:N0);
    end
end
end



function [X_c,X_Tn_c,z_c,N_inner_c,resvec_c,relres_c,cond_c] = update_iterate(weight_op,...
        meas_op,n,y,...
        x_prev_newbase,z_start,...
        lambda,tol_outer,tol_CG,mode_linsolve,...
        first_iterate,N0_inner,...
        type_mean,qmean_para,tangent_para,AdjInvVersion,varargin)
weight_op.alpha = 0;
U = weight_op.U;
%U';
V = weight_op.V;
r_c = size(U,2);
d1 = n.d1;
d2 = n.d2;
Wl_inv      = weight_op.Wl_inv;
Wl_inv_eps  = weight_op.Wl_inv_eps;
eps_c       = weight_op.eps;

if weight_op.symmetricflag
    dim_tangent = r_c*((r_c+1)/2+d1);
else
    dim_tangent = r_c*(d1+d2+r_c);
end

%%% Update tolerance parameter of the conjugate gradient method.
tol_CG=max(tol_CG,1e-16);

%%% Load additional variables
X_c4 = varargin{2};

if first_iterate
%     X_c={zeros(r_c,r_c),zeros(r_c,d2),zeros(d1,r_c),...
%         Wl_inv_eps.*y./(Wl_inv_eps+lambda)};
    X_Tn_c = zeros(dim_tangent,1);
    z_c = [];
    N_inner_c = 0;
    relres_c = 0;
    resvec_c = 0;
    cond_c = 0;
    if isfield(weight_op,'alpha') && weight_op.alpha
        fac = 1/meas_op.mu;%(weight_op.alpha+1)^(-1);
    else
        fac = 1;
    end
    y_used = y;
else
    if contains(mode_linsolve,'fullspace')
        PhWmin1Phstar = 1;
    elseif contains(mode_linsolve,'rangespace')
        y_used = y;
        if strcmp(mode_linsolve,'rangespace_matrix')
            [H,~,~,dH] = weight_op_lowrank_prec(U,V,[Wl_inv;Wl_inv_eps],1,-1,...
            type_mean,0,[],[],qmean_para);
            H{1}=H{1}-Wl_inv_eps;
            dH{1}=dH{1}(1:r_c)-Wl_inv_eps;
            dH{2}=dH{2}(1:r_c)-Wl_inv_eps;

            weight_vec_c = get_weight_vec(d1,d2,H,dH,tangent_para,...
            weight_op,0);
            Wmatinv = kron((V.*Wl_inv.')*V'+Wl_inv_eps,U*(Wl_inv.*U')+Wl_inv_eps);
            PhWmin1Phstar= meas_op.Phi*Wmatinv*meas_op.Phi';
            z_c = PhWmin1Phstar\y_used;
%                 x_c = Wmatinv*(meas_op.Phi'*z_c);
            gam = proj_tangspace_from_Oprange(z_c,meas_op,weight_op,'tangspace');
            gam = apply_weight_vec(gam,weight_vec_c,0);
            X_Tn_c = gam;% -X_Tn_c_delta;
            cond_c  = 0;
        else
            [H,~,~,dH] = weight_op_lowrank_prec(U,V,[Wl_inv;Wl_inv_eps],1,-1,...
            type_mean,0,[],[],qmean_para);
            H{1}=H{1}-Wl_inv_eps;
            dH{1}=dH{1}(1:r_c)-Wl_inv_eps;
            dH{2}=dH{2}(1:r_c)-Wl_inv_eps;
            weight_vec_c = get_weight_vec(d1,d2,H,dH,tangent_para,...
                weight_op,0);
            PhWmin1Phstar_handle = @(z) PhWmin1Phstar_rangespace(z,lambda,weight_vec_c,...
                    Wl_inv_eps,weight_op,meas_op);
            if first_iterate
                z_start = zeros(length(y),1);
            else
                if isempty(z_start)
                    tmp = apply_weight_vec(x_prev_newbase,weight_vec_c.^(-1),...
                         0);
                    z_start_test = proj_Oprange_tangspace(tmp,meas_op,weight_op,...
                        'tangspace',0);
                    z_start = z_start_test + Wl_inv_eps^(-1).*X_c4;
                end
            end
            prec = @(z) z;
            [z_c,~,relres_c,N_inner_c,resvec_c] = pcg_1(PhWmin1Phstar_handle,y,tol_CG,...
                N0_inner,prec,[],z_start);
            gam = proj_tangspace_from_Oprange(z_c,meas_op,weight_op,'tangspace');
            gam = apply_weight_vec(gam,weight_vec_c,0);
    %         X_Tn_c_delta =  proj_tangspace_from_Oprange(Wl_inv_eps.*z_c,'MatrixCompletion',...
    %                         U,V,meas_op.sps_plc,'tangspace');
            X_Tn_c = gam;% -X_Tn_c_delta;
            cond_c  = 0;
        end
    elseif contains(mode_linsolve,'tangspace')
        [H,~,~,dH] = weight_op_lowrank_prec(U,V,[Wl_inv;Wl_inv_eps],Wl_inv_eps,1,...
        type_mean,0,[],[],qmean_para);
        PTPhstarPhPTstar_handle = @(gam) PTPhstarPhPTstar_tangspace(gam,...
            weight_op,meas_op,AdjInvVersion);
%         if first_iterate
%              H{1}=lambda.*H{1};
%              dH{1}=lambda.*dH{1}(1:r_c);
%              dH{2}=lambda.*dH{2}(1:r_c);
%              weight_vec_c = get_weight_vec(d1,d2,H,dH,tangent_para,increase_antisymmetricweights);
%         else
        if contains(mode_linsolve,'tangspace_largeeps')
            if lambda > 0
                error('This implementation does not work yet for lambda > 0.')
            end
            dH{1}=dH{1}(1:r_c);
            dH{2}=dH{2}(1:r_c);
            weight_vec_c = get_weight_vec(d1,d2,H,dH,tangent_para,...
                weight_op,0);
            weight_vec_c_largeeps =  1-weight_vec_c;
            if strcmp(tangent_para,'intrinsic')
                M = @(gam) apply_weight_vec(E2D_rankmanifold(PTPhstarPhPTstar_handle(...
                    D2E_rankmanifold(apply_weight_vec(gam,weight_vec_c_largeeps.^(1/2),...
                0),QU,QV)),QU,QV),...
                weight_vec_c_largeeps.^(1/2),increase_antisymmetricweights,d1,d2) + ...
                apply_weight_vec(gam,weight_vec_c,0);
            else
                M = @(gam) apply_weight_vec(PTPhstarPhPTstar_handle(apply_weight_vec(gam,weight_vec_c_largeeps.^(1/2),...
                0)),weight_vec_c_largeeps.^(1/2),0) + ...
                apply_weight_vec(gam,weight_vec_c,0);
            end
        else
            if strcmp(type_mean,'left-sided') && weight_op.alpha == 0 && lambda == 0
                H{1} = ((H{1}./Wl_inv_eps).^(-1)-Wl_inv_eps).^(1/2);%(H{1}./Wl_inv_eps).^(-1);
                dH{1} = ((dH{1}(1:r_c)./Wl_inv_eps).^(-1)-Wl_inv_eps).^(1/2);
                dH{2} = zeros(r_c,1);
            else
                H{1}=(weight_op.alpha*Wl_inv_eps+((Wl_inv_eps+lambda)/Wl_inv_eps).*H{1})./(1-H{1});
                dH{1}=(weight_op.alpha*Wl_inv_eps+((Wl_inv_eps+lambda)/Wl_inv_eps).*dH{1}(1:r_c))./(1-dH{1}(1:r_c));
                dH{2}=(weight_op.alpha*Wl_inv_eps+((Wl_inv_eps+lambda)/Wl_inv_eps).*dH{2}(1:r_c))./(1-dH{2}(1:r_c));
            end
            weight_vec_c = get_weight_vec(d1,d2,H,dH,tangent_para,...
                weight_op,0);
            if strcmp(type_mean,'left-sided') && weight_op.alpha == 0 && lambda == 0
                    M = @(gam) apply_weight_vec(PTPhstarPhPTstar_handle(apply_weight_vec(gam,weight_vec_c,...
                    0)),weight_vec_c,...
                    0) + ...
                    Wl_inv_eps.*gam;
            else
                if strcmp(tangent_para,'intrinsic')
                    M = @(gam) E2D_rankmanifold(PTPhstarPhPTstar_handle(D2E_rankmanifold(gam,QU,QV)),QU,QV) + ...
                    apply_weight_vec(gam,weight_vec_c,0);
                else
                    M = @(gam) PTPhstarPhPTstar_handle(gam) + ...
                    apply_weight_vec(gam,weight_vec_c,0);
                end
            end
        end
        %%%%%%%%%%%%%%%%%%
        % ytmp = y;
        y_used = y;
        ytmp = applyMeasopForwardAdjointInv(y_used,meas_op,AdjInvVersion);
        rhs = proj_tangspace_from_Oprange(ytmp,meas_op,weight_op,'tangspace',0);
        if strcmp(type_mean,'left-sided') && weight_op.alpha == 0 && lambda == 0
            rhs = apply_weight_vec(rhs,weight_vec_c,0);
        end
%             tmp2 = meas_op.kappa(meas_op.kappainv(meas_op.kappainv(reshape(meas_op.Samp_opmat'*y./4,size(U,1),size(U,1)))));%             tmp2 = applyMeasopBackward(y_used,meas_op);% %             tmp2 = (-2).*meas_op.kappa(meas_op.kappainv(reshape(meas_op.Samp_opmat'*y./4,size(U,1),size(U,1))));%             rhs = matrixspace_to_tangspace(reshape(tmp2,size(U,1),size(U,1)),U,V);        
        if strcmp(tangent_para,'intrinsic')
            rhs = E2D_rankmanifold(rhs,QU,QV);
        end
        if contains(mode_linsolve,'tangspace_largeeps')
            rhs = apply_weight_vec(rhs,weight_vec_c_largeeps.^(1/2),0);
            prec = @(z) z;
        else
            prec = @(z) apply_weight_vec(z,(ones(dim_tangent,1)+(weight_vec_c)).^(-1),0);
        end
        if strcmp(tangent_para,'intrinsic')
            x_prev_newbase = E2D_rankmanifold(x_prev_newbase,QU,QV);
        end
        if strcmp(mode_linsolve,'tangspace_largeeps')
%                 gamma_start = x_prev_newbase;
            gamma_start= apply_weight_vec(x_prev_newbase,...
                weight_vec_c_largeeps.^(-1/2),0);
        else
%           gamma_start = apply_weight_vec(x_prev_newbase,d1,d2,...
%                     (((Wl_inv_eps+lambda)/Wl_inv_eps)^(-1).*weight_vec_c+ones(dim_tangent,1)).^(-1),...
%                     increase_antisymmetricweights,U,V);
            gamma_start = x_prev_newbase;
        end
          % ((Wl_inv_eps.*weight_vec_c)./(Wl_inv_eps+lambda)+ones(dim_tangent,1)).^(-1),...
    %             gamma_start = gamma_old;
    %             gamma_start = proj_tangspace(gamma_old,'MatrixCompletion',...
    %                   U_X_m1,V_X_m1,U_X_c,V_X_c);
    %             gamma_start = proj_tangspace(x_Tn_c_old,'MatrixCompletion',...
    %                 U_X_m1,V_X_m1,U_X_c,V_X_c);
    %             gamma_start = ((weight_vec_c)+ones(dim_tangent,1)).^(-1).*gamma_start;
    %        gamma_start=proj_tangspace(gamma_start,'MatrixCompletion',U_X_c,V_X_c,U_X_c,V_X_c);
           %         gamma_start = x_prev_newbase;
                    %((weight_vec_c)+ones(dim_tangent,1)).^(-1).*x_prev_newbase;
    %             gamma_start = zeros(dim_tangent,1);
            %gamma_old;%zeros(dim_tangent,1); % prec = @(z) apply_weight_vec(z,d1,d2,(Wl_inv_eps+lambda)./weight_vec_c,...
%             increase_antisymmetricweights);
    %     if verbose <= 1
        if contains(mode_linsolve,'exact')
            M_mat = zeros(length(rhs),length(rhs));
            for i=1:length(rhs)
                e_vec = zeros(length(rhs),1);
                e_vec(i) = 1;
                M_mat(:,i)=M(e_vec);
            end
            gamma = M_mat\rhs;
            relres_c = 0;
            N_inner_c = 1;
            resvec_c = 0;
            cond_c  = 0;
        else
            cond_c  = 0;
%             rhs_explained = M(gamma_start);
%             [gamma,~,relres_c,N_inner_c,resvec_c] = pcg_1(M,rhs,tol_CG,N0_inner,prec,[]);
% %             gamma = gamma_start+gamma_delta;
            % prec = @(z) z;
            rhs_explained = M(gamma_start);
            [gamma_delta,~,relres_c,N_inner_c,resvec_c] = pcg_1(M,rhs-rhs_explained,tol_CG,N0_inner,prec,[]);
            gamma = gamma_start+gamma_delta;
            if strcmp(type_mean,'left-sided') && weight_op.alpha == 0 && lambda == 0
                gamma = apply_weight_vec(gamma,weight_vec_c,0);
            end
%             [gamma,~,relres_c,N_inner_c,resvec_c] = pcg_1(M,rhs,tol_CG,N0_inner,prec,[]);
        end
    %     else
    %         [gamma,~,~,N_inner_c,resvec_c] = pcg_1(M,rhs,tol_CG,N0_inner,prec,[],gamma_start);
    %     end
    %      error_gam = check_intangspace(gamma,U,V)
        if not(first_iterate)
            if contains(mode_linsolve,'tangspace_largeeps')
                gamma = apply_weight_vec(gamma,weight_vec_c_largeeps.^(1/2),...
                    0);
            end
            if strcmp(tangent_para,'intrinsic')
                gamma = D2E_rankmanifold(gamma,QU,QV);
            end
        end
        X_Tn_c = gamma;
        z_c= [];
    end
end
[X_c] = prepare_iterate_output(X_Tn_c,z_c,y_used,...
    meas_op,weight_op,Wl_inv_eps,lambda,mode_linsolve,first_iterate);
end

function [X_c] = prepare_iterate_output(X_Tn_c,z_c,y_used,...
    meas_op,weight_op,Wl_inv_eps,lambda,mode_linsolve,...
    first_iterate)

global PhiPhiT_inv_y_c
% X_Tn_c = gamma;
% X_Tn_c = gamma;
[X_c,X_Tn_c_1,X_Tn_c_2,X_Tn_c_3] = get_Tk_matrices(X_Tn_c,weight_op);
% X_c.Gam1 = X_Tn_c_1;
% X_c.Gam2 = X_Tn_c_3;
% X_c.Gam3 = X_Tn_c_2;
if isfield(weight_op,'alpha') && weight_op.alpha
    fac = 1/meas_op.mu;%(weight_op.alpha+1)^(-1);
else
    fac = 1;
end
if first_iterate
    PhiPhiT_inv_y_c = zeros(length(y_used),1);
    X_c.res_range = applyMeasopForwardAdjointInv(y_used,meas_op).';
else
    if contains(mode_linsolve,'tangspace')
        res_c = y_used-proj_Oprange_tangspace(X_Tn_c,meas_op,weight_op,'tangspace',0);
        X_c.res_range = Wl_inv_eps.*applyMeasopForwardAdjointInv(res_c,meas_op).'./(Wl_inv_eps+lambda);
    elseif strcmp(mode_linsolve,'rangespace')
        X_c.res_range = Wl_inv_eps.*z_c.';
    end
end
X_c.symmetricflag = weight_op.symmetricflag;
end

function mode_linsolve=use_optimal_cg_mode(quot_r_c)
%%%% Adaptively choose which linear system to solve (among ones
%%%% that are all equivalent to each other)
quot_r_threshold = 1.1; 
% (those values are subject to change, but they work well already)

if quot_r_c < quot_r_threshold
    mode_linsolve = 'rangespace';
else
    mode_linsolve = 'tangspace';
end
end

function eps_c = update_epsilon(eps_c,mode_eps,...
    oracle_knowledge,n,epsmin,use_min,p,varargin)
sing_c   = varargin{1};
k        = varargin{6};
%%% Keep this below if second epsilon update allows for eps increase
% if k == 2
%    use_min = false; 
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin >= 8
    singmin = varargin{8};
end

if strcmp(mode_eps,'oracle_model_order')
    if oracle_knowledge
        r = varargin{2};
        if p == 1
            eps_rule = sing_c(r+1);
        else
            eps_rule = sing_c(r+1);%/min(n.d1,n.d2)^(p/(2-p));
            if nargin >= 8 && not(isempty(singmin))
                eps_rule = max(eps_rule,max(singmin));
            end
        end
    else
        error('This epsilon rule is only possible if oracle knowledge is available.')
    end
elseif strcmp(mode_eps,'oracle_ellinfty')
    if oracle_knowledge
        r = varargin{2};
        if p == 1
            eps_rule = sing_c(r+1)/min(n.d1,n.d2);
        else
            eps_rule = sing_c(r+1)/min(n.d1,n.d2)^(p/(2-p));
            if nargin >= 8 && not(isempty(singmin))
                eps_rule = max(eps_rule,max(singmin));
            else
            end
        end
    else
        error('This epsilon rule is only possible if oracle knowledge is available.')
    end
elseif contains(mode_eps,'oracle_ell1')
    if oracle_knowledge
        r = varargin{2};
        sX_c_complement = sing_c(r+1:end);
        if strcmp(mode_eps,'oracle_ell1')
            eps_rule = norm(sX_c_complement,1)/((min(n.d1,n.d2))^(1/2));
        elseif strcmp(mode_eps,'oracle_ell1_nodiv')
            eps_rule = norm(sX_c_complement,1)/((min(n.d1,n.d2))^(0));
        elseif strcmp(mode_eps,'oracle_ell1_gooddiv')
            eps_rule = norm(sX_c_complement,1)/((min(n.d1,n.d2))^(1));%/(2-p)));
        end
    else
        error('This epsilon rule is only possible if oracle knowledge is available.')
    end
elseif contains(mode_eps,'oracle_ell2')
    if oracle_knowledge
        r = varargin{2};
        sX_c_complement = sing_c(r+1:end);
        if strcmp(mode_eps,'oracle_ell2')
            eps_rule = norm(sX_c_complement,2)/(min(n.d1,n.d2)^(1/2));
        elseif strcmp(mode_eps,'oracle_ell2_nodiv')
            eps_rule = norm(sX_c_complement,2)/(min(n.d1,n.d2)^(0));
        elseif strcmp(mode_eps,'oracle_ell2_gooddiv')
            eps_rule = norm(sX_c_complement,2)/(min(n.d1,n.d2)^((1/(2-p))-1/2));
        end
    else
        error('This epsilon rule is only possible if oracle knowledge is available.')
    end
elseif strcmp(mode_eps,'auto_decay')
        eps0 = varargin{4};
        tol_outer = varargin{5};
        N0   = varargin{7};
        
        eps_rule = eps0.*(tol_outer/eps0)^(((k-1)/(N0-1))^(2-p));
elseif strcmp(mode_eps,'groundtruth_S1')
        X0 = varargin{9};
        Xk = varargin{10};
        eps_rule = norm(svd(Xk-X0),2)/10;

elseif strcmp(mode_eps,'auto_decay2')
        eps0 = varargin{4};
        eps_rule = min(eps0.*(eps_c)^(2-p),min(eps0.*(eps_c/eps0)^(2-p),eps0));
elseif strcmp(mode_eps,'iter_diff')
        diff_c = varargin{3};
        
        eps_rule = 0.9.*diff_c/sqrt(min(n.d1,n.d2));
end

if use_min
   eps_c=max(min(eps_rule,eps_c),epsmin); 
else
   eps_c=max(eps_rule,epsmin);
end
end


function z_new = PhWmin1Phstar_rangespace(z,lambda,weight_vec_c,...
            Wl_inv_eps,weight_op,...
            meas_op)
% To Do (Apr 14, 2021): Adapt this to meas_op API.


if 0
    U = weight_op.U;
    V = weight_op.V;
    d1 = size(U,1);
    d2= size(V,1);
    gam = proj_tangspace_from_Oprange(z,meas_op,...
                    weight_op,'tangspace',0);
    gam = apply_weight_vec(gam,weight_vec_c,0,d1,d2);
    z1 = proj_Oprange_tangspace(gam,meas_op,weight_op,...
        'tangspace',0);
else
    gam = proj_tangspace_from_Oprange(z,meas_op,weight_op,'tangspace');
    gam = apply_weight_vec(gam,weight_vec_c,0);
    z1 = proj_Oprange_tangspace(gam,meas_op,weight_op,'tangspace');
end
%
%z_new = z1 + (lambda+Wl_inv_eps).*z;
z_new = z1 + (lambda+Wl_inv_eps).*applyMeasopForwardAdjoint(z,meas_op);
end

function [gam_new,tmp_AdjInv] = PTPhstarPhPTstar_tangspace(gam,weight_op,...
            meas_op,AdjInvVersion,tmp_AdjInv_old)
        
% if isfield(meas_op,'sps_plc')
%     sps_plc = meas_op.sps_plc;
%     rowind  = meas_op.rowind;
%     colind  = meas_op.colind;
%     m=length(rowind); 
% else
%     Phi     = meas_op.Phi;
%     m       = size(Phi,1);
% end
% [d1,R] = size(U);
% d2= size(V,1);
if strcmp(meas_op.problem_type,'MultidimensionalScaling')
%     gam_new = PTstarPhKminPhstarPT(gam,meas_op,U,V) - ...
%         PTstarGGstarPT(gam,meas_op,U,V);
%     weight_op = struct;
%     weight_op.U = U;
%     weight_op.V = V;
%     test_proj_tangspace(gam,weight_op,weight_op,meas_op);
%     gam_new = test_proj_tangspace(PTstarPhKminPhstarPT(...
%         test_proj_tangspace(gam,weight_op,weight_op,meas_op),meas_op,U,V),weight_op,weight_op,meas_op) - ...
%         PTstarGGstarPT(gam,meas_op,U,V);
%     gam_new = gam+(PTstarPhKminPhstarPT(gam,meas_op,U,V)) - ...
%         PTstarGGstarPT(gam,meas_op,U,V);

    % gam_new = (PTstarPhKminPhstarPT(gam,meas_op,weight_op)) + ...
    %     PTstarGGstarperpPT(gam,meas_op,weight_op);

    [gam_new_1] = PTstarPhKminPhstarPT(gam,meas_op,weight_op,AdjInvVersion);
    % gam_new = gam_new_1 + PTstarGGstarperpPT(gam,meas_op,weight_op);
    gam_new = gam_new_1;
%     gam_new = 0.5.*(PTstarPhKminPhstarPT(gam,meas_op,U,V) +PTstarPhKminPhstarPT_transp(gam,meas_op,U,V)) - ...
%         PTstarGGstarPT(gam,meas_op,U,V);
else
    tmp = proj_Oprange_tangspace(gam,meas_op,weight_op,...
    'tangspace',0);
    if strcmp(meas_op.problem_type,'CommunityDetection') || strcmp(meas_op.problem_type,'MaxCut fixedval') ...
            || strcmp(meas_op.problem_type,'MultidimensionalScaling')
        tmp = applyMeasopForwardAdjointInv(tmp,meas_op);%meas_op.PhiPhiT_inv*tmp;
    end
    gam_new = proj_tangspace_from_Oprange(tmp,meas_op,...
            weight_op,'tangspace',0);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% only needed for MultidimensionalScaling %%%%
function  [gam_new] = PTstarPhKminPhstarPT(gam,meas_op,weight_op,AdjInvVersion)
 
    tmp = proj_Oprange_tangspace(gam,meas_op,weight_op,'tangspace',0);
    tmp_AdjInv = applyMeasopForwardAdjointInv(tmp,meas_op,AdjInvVersion);%     tmp = meas_op.T\tmp;
    gam_new = proj_tangspace_from_Oprange(tmp_AdjInv,meas_op,...
                weight_op,'tangspace',0);

    % tmp2 = applyMeasopBackward(tmp_AdjInv,meas_op);
    % gam_new = matrixspace_to_tangspace(reshape(tmp2,...
    %     size(weight_op.U,1),size(weight_op.U,1)),weight_op);
    % norm(gam_new-gam_new_test)
end


function gam_new = PTstarGGstarperpPT(gam,meas_op,weight_op)
% X = tangspace_to_matrixspace(gam,U,V);
% % tmp = meas_op.kappa(meas_op.kappainv(meas_op.kappainv(meas_op.kappa(X))));
% % tmp = (-2).*meas_op.kappa(meas_op.kappainv(meas_op.kappa(X)));
% tmp = Proj_GGstarPerp(X);
% gam_new = matrixspace_to_tangspace(tmp,U,V);

U      = weight_op.U;
d1 = size(U,1);
% norm(weight_op.U-weight_op.V)
if ~weight_op.symmetricflag
    V      = weight_op.V;
    d2      =size(V,1);
else
    V = U;
    d2 = d1;
end
[~,M1,M2,M3] = get_Tk_matrices(gam,weight_op);

LeftMat1 = U*M1+M2;
sum1 = sum(LeftMat1,1)*(V.'./d1);
sum1 = sum1 + sum(U,1)*(M3./d1);
sum2 = LeftMat1*(sum(V,1).'./d2);
sum2 = sum2 + U*(sum(M3,2)./d2);
sumsum = sum(sum2)./d1;

gam_new = matrixspace_to_tangspace({ones(d1,1),sum1.'-sumsum.*ones(d2,1)},weight_op);
gam_new = gam_new + matrixspace_to_tangspace({sum2,ones(d2,1)},weight_op);
% fprintf(['Orig norm: %.4e, After norm %.4e: \n'],norm(gam),norm(gam_new))
gam_new = zeros(length(gam),1);
end

function Proj_X = Proj_GGstar(X)
sum2 = sum(X,2)./size(X,1);
sum1 = sum(X,1)./size(X,1);
sumsum = sum(sum2)./size(X,1);
Proj_X = X - sum2*ones(1,size(X,1))- ones(size(X,1),1)*sum1 + sumsum;
end

function Proj_X = Proj_GGstarPerp(X)
sum2 = sum(X,2)./size(X,1);
sum1 = sum(X,1)./size(X,1);
sumsum = sum(sum2)./size(X,1);
% Proj_X = X - sum2*ones(1,size(X,1))- ones(size(X,1),1)*sum1 + sumsum;
Proj_X = + sum2*ones(1,size(X,1))+ ones(size(X,1),1)*sum1 - sumsum;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [N0,N0_inner,N_SVD,p,tol_CG,tol_outer,objective,...
    type_mean,qmean_para,mode_linsolve,adaptive_cg,mode_eps,epsmin,...
    use_min,tracking,verbose,recsys,saveiterates,...
    tangent_para,forcefull] = get_options(opts)

N0     = opts.N0;       % Maximal number of iterations for the IRLS algorithm
N0_inner=opts.N0_inner; % Maximal number of inner iterations (e.g. in the CG method)
N_SVD  = opts.N_SVD;    % max. number of iterations for power method-type solver for partial SVD (such as bksvd)
p      = opts.p;        % Parameter describing the non-convexity of the objective function
                        % (p = 0 for logarithmic objective, 0 < p < 1 for
                        % Schatten-p quasi-norm objective)
tol_CG = opts.tol_CG;   % Tolarance used in the stopping criterion of the inner CG method
tol_outer = opts.tol;   % Tolerance used as stopping criterion of outer IRLS iterations
                        % (min. change of relative Frobenius norm)
if not(isfield(opts,'objective')) || isempty(opts.objective)
    objective = 'objective_thesis';
else
    objective   = opts.objective;% objective = 'objective_thesis': Corresponds to using objectives
end
%                           of the type:
%                           log(max(sigma,eps))+0.5.*(min(sigma,eps)^2/eps^2-1).
%                           (as used in the Ph.D. thesis [Kuemmerle
%                           2019])
                        % mode = 'pluseps': Corresponds to using objectives
%                           of the type:
%                           log(max(sigma,eps)+eps)+0.25.*(min(sigma,eps)^2/eps^2-1)
                        % mode = 'pluseps_squared':Corresponds to using objectives
%                           of the type:
%                           log(max(sigma,eps)^2+eps^2)+0.5.*(min(sigma,eps)^2/eps^2-1)
type_mean = opts.type_mean;% type_mean = 'harmonic': corresponds to usage 
                        % of harmonic mean of weight matrices, as Kuemmerle
if strcmp(type_mean,'qmean')
    qmean_para = opts.qmean_para; 
    % parameter for the q-mean: = 1 for arithmetic mean,
    %                           = 0 for heometric mean,
    %                           = -1 for harmonic mean,
    %                           = -infty for "min" mean.
elseif strcmp(type_mean,'optimal')
    if p == 0
        type_mean = 'geometric';
        qmean_para = [];
    else
        type_mean = 'qmean';
        qmean_para = p/(p-2);
    end
else
    qmean_para = [];
end
mode_linsolve = opts.mode_linsolve;
adaptive_cg = opts.adaptive_cg;

 %(just incorporated for mode = 'fullspace' yet)    % and Sigl (2017).
                        %    type_mean = 'Wasserstein': corresponds to a
                        %    Wasserstein mean of the inverses of the
                        %    one-sided weight matrix inverses, cf. Bhatia,
                        %    Jain, Lim (2018)
mode_eps   = opts.mode_eps;%= 'classical_r': eps choice dependent on true
%                                           rank r, as classical
%                          = 'auto_decay' : eps that is automatically
%                             superlinearly decaying, rule such that eps0=
%                             sing_val_1(x^(1)) and 
%                             eps(iter) = 10^(-(iter/N0)^(2-p)*10)
%                          = 'iter_diff'  : eps similar to Voronin &
%                                           Daubechies (2016), i.e.,
%                             eps(iter)= min(eps(iter-1),\|x{iter}-x{iter-1}\|_2
%                                      + eps0*10^(-(iter/10)*(2-p)*10))
epsmin     = opts.epsmin;  
use_min    = opts.use_min; % == 1: In the epsilon choice, use the minimum 
                           %        of the previous one and the one calculated
                           %        by the rule.
                           % == 0: Use the one calculated by the epsilon rule.
if isfield(opts,'saveiterates') && opts.saveiterates == 1
    saveiterates = 1;
else
    saveiterates = 0;
end

if isfield(opts,'recsys') && opts.recsys == 1
    recsys = 1;
else
    recsys = 0;
end
if isfield(opts,'efficiencymode') && strcmp(opts.efficiencymode,'fast')
    forcefull = false;
else
    forcefull = true;
end

if isfield(opts,'tracking') &&  ~isempty(opts.tracking)
    tracking=opts.tracking;
else
    tracking = 0;
end

verbose     = opts.verbose; % verbose = 0: No output except of last iterate
                            % verbose = 1: Return also intermediate
                            % iterates and a few algorithmic parameters 
                            % verbose = 2: Return also many further
                            % parameters
                            
if isfield(opts,'tangent_para') &&  ~isempty(opts.tangent_para)
    tangent_para = opts.tangent_para;
else
    tangent_para = 'extrinsic';
end

end

function  [n,meas_op,Yval]=get_problem_parameters(prob,opts)

n.d1   = prob.d1;       % First dimension of matrix X to be recovered
n.d2   = prob.d2;       % Second dimension of matrix X to be recovered

meas_op.problem_type = prob.problem_type;
if isfield(opts,'efficiencymode') && strcmp(opts.efficiencymode,'fast')
    meas_op.efficiencymode = 'fast';
elseif ~isfield(opts,'efficiencymode')
    meas_op.efficiencymode = [];
else
    meas_op.efficiencymode = opts.efficiencymode;
end
Yval   = prob.y;

% [W_sampled,W_sampled_long,meas_op_out,ind_diag] = create_Wmatrix_MDS(prob.Phi,n);
% meas_op_out.T = W_sampled'*W_sampled;
if isfield(opts,'efficiencymode') && strcmp(opts.efficiencymode,'fast')
    meas_op_out.sps_meas = prob.Phi;
    meas_op_out.Omega = find(meas_op_out.sps_meas);
end
meas_op_out.sps_upper = prob.Phi;
meas_op_out.problem_type = meas_op.problem_type;
meas_op_out.efficiencymode = meas_op.efficiencymode;   
meas_op_out.mode_linsolve = opts.mode_linsolve;
meas_op = meas_op_out;
[rowind,colind] = find(prob.Phi); 
Yval   = prob.y;   
Omega  = find(prob.Phi);
meas_op.Omega  = Omega;
meas_op.m = length(Omega);
meas_op.rowind = rowind;
meas_op.colind = colind;
meas_op.OmegaT = sub2ind([n.d2,n.d1],meas_op.colind,meas_op.rowind) ;
meas_op.sps_plc = sparse(rowind,colind,ones(length(Omega),1),n.d1,n.d2);
if not(strcmp(opts.efficiencymode,'fast'))
    meas_op.Phi = sparse(1:length(Omega),Omega,ones(length(Omega),1),length(Omega),n.d1*n.d2);
end
meas_op.rowind_diag_set = cell(1,n.d1);
meas_op.colind_diag_set = cell(1,n.d1);
for l = 1:n.d1
    meas_op.rowind_diag_set{l} = find(meas_op.rowind==l)';
    meas_op.colind_diag_set{l} = find(meas_op.colind==l)';
    meas_op.unionind_diag_set{l} = [meas_op.rowind_diag_set{l},meas_op.colind_diag_set{l}];
end
maxx_ind = 1;
for l = 1:n.d1
    maxx_ind = max(length(meas_op.unionind_diag_set{l}),maxx_ind);
end
meas_op.IndMat = (meas_op.m+1).*ones(maxx_ind,n.d1);
for l = 1:n.d1
    meas_op.IndMat(1:length(meas_op.unionind_diag_set{l}),l) = meas_op.unionind_diag_set{l};
end
%meas_op.T_handle = @(y) apply_T_handle(y,meas_op,n);
if not(strcmp(opts.efficiencymode,'fast'))
    [W,W_long,meas_op_W,ind_diag,V] = create_Wmatrix_MDS(prob.Phi,n);
    meas_op.Phi = meas_op_W.Phi;
    meas_op.PhiPhiT = meas_op.Phi*meas_op.Phi';
    diag_indices = find(speye(n.d1));
    Omega_ext = [Omega;diag_indices];
    VE_Omega= V(:,Omega_ext);
    meas_op.PhiVPhiVT = VE_Omega'*VE_Omega;
end
% diag_indices = find(speye(n.d1));
% Omega_ext = [Omega;diag_indices];
meas_op.T_handle = @(y) apply_T_handle(y,meas_op);
% WE_Omega= W_long;
% Identity_test = WE_Omega'*VE_Omega;
% PhiPhiT_short=W'*W;
% dim=length(Omega_ext);
% PhiVPhiVT_short = meas_op.PhiVPhiVT(1:dim,1:dim)-1/n.d1^2;
meas_op.PhiVPhiVT_handle = @(y) apply_PhiVPhiVT_handle(y,meas_op);
% MatT=extract_matrix_from_handle(meas_op.PhiVPhiVT_handle,dim);
end

function outp=apply_PhiVPhiVT_handle(y,meas_op)
n = size(meas_op.sps_plc,1);
% entry_sums = zeros(n.d1,1);
% for l=1:n.d1
%     entry_sums(l) = sum(y(meas_op.unionind_diag_set{l}));
% end
yupper= y(1:meas_op.m);
ylower= y(meas_op.m+1:end);
avg1=sum(yupper);
avg2=sum(ylower);
yt = y;
yt(meas_op.m+1) = 0;
% yt=sparse(yt);
entry_sums = sum(yt(meas_op.IndMat),1)';


outp_u = -(entry_sums(meas_op.rowind)+entry_sums(meas_op.colind))./(2*n)+(1/2).*yupper+avg1/n^2;
outp_l = (2/n).*ylower-avg2/n^2;
outp = [outp_u;outp_l];
% for l=1:meas_op.m
%    T_y(l) = entry_sums(meas_op.rowind(l))+entry_sums(meas_op.colind(l))+2.*y(l);
% end
end

function [T_y]=apply_T_handle(y,meas_op)
% entry_sums = zeros(n.d1,1);
% for l=1:n.d1
%     entry_sums(l) = sum(y(meas_op.unionind_diag_set{l}));
% end
yt = y;
yt(meas_op.m+1) = 0;
% yt=sparse(yt);
entry_sums = sum(yt(meas_op.IndMat),1)';

T_y = entry_sums(meas_op.rowind)+entry_sums(meas_op.colind)+2.*y;
% for l=1:meas_op.m
%    T_y(l) = entry_sums(meas_op.rowind(l))+entry_sums(meas_op.colind(l))+2.*y(l);
% end
end

function [U,sing,V,singmin,Uperp,Vperp,symTmat,symInd] = get_weight_update(X_c_handle,d,r_c,r,...
    oracle_knowledge,N_SVD,opts,varargin)
decomposition = 'svd';
Uinit = [];
Vinit = [];
if ~isempty(varargin)
    for tt = 1:2:length(varargin)
        switch lower(varargin{tt})
            case 'decomposition'
                decomposition = varargin{tt+1};
        end
        switch varargin{tt}
            case 'Uinit' 
                Uinit = varargin{tt+1};
            case 'Vinit'
                Vinit = varargin{tt+1};
        end
    end
end

if oracle_knowledge
    [U,sing,V,singmin] = update_matrixdec_from_handle(X_c_handle,...
        min(d,max(r_c,r+1)),N_SVD,decomposition,opts,'Uinit',Uinit,'Vinit',Vinit);
else
    [U,sing,V,singmin] = update_matrixdec_from_handle(X_c_handle,...
        min(d,r_c),N_SVD,decomposition,opts,'Uinit',Uinit,'Vinit',Vinit);
end
if isfield(opts,'tracking') && opts.tracking
    Uperp = null(U');
    Vperp = null(V');
else
    Uperp = [];
    Vperp = [];
end
if opts.symmetricflag
    symTmat = sparse(ones(r_c,r_c));
    symTmat = triu(symTmat,0);
    symInd = find(symTmat);
else
    symTmat = [];
    symInd = [];
end
end


function [weight_op,X_c,eps_c,r_c,r_upper,first_iterate,X_c_newbase] = MatrixIRLS_initialize_weights(n,...
    opts,oracle_knowledge,r,y,meas_op,prob)
if isfield(opts,'initialization')
    initialization = opts.initialization; % = []: No prior knowledge, use uniform weights.
                                      % otherwise use:
                                      % weight_op.U = initialization.U 
                                      % (d1 x r_c) matrix with orthonormal
                                      % columns, left singular vector
                                      % matrix.
                                      % weight_op.V = initialization.V 
                                      % (d2 x r_c) matrix with orthonormal
                                      % columns, right singular vector
                                      % matrix.
                                      % weight_op.sing =
                                      % initialization.sing
                                      % (r_c x 1) vector with
                                      % non-increasing entries. Singular
                                      % values.
                                      % eps_0 = initialization.eps
                                      % initial smoothing parameter epsilon
else
   initialization = []; 
end
weight_op.symmetricflag = opts.symmetricflag;
if isfield(opts,'alpha')
    weight_op.alpha = opts.alpha;
end
objective = 'objective_thesis';
if isfield(opts,'R')
    r_upper = opts.R; % Upper bound for number of singular values to be calculated
else
    r_upper = Inf;
end
d = min(n.d1,n.d2);
r_upper = min(r_upper,min(n.d1,n.d2));
if not(isempty(initialization))
    r_c = length(initialization.sing);
else
    if oracle_knowledge
        r_c = r;
    else
        r_c = r_upper;
    end
end
if opts.symmetricflag
    weight_op.symTmat = sparse(ones(r_c,r_c));
    weight_op.symTmat = triu(weight_op.symTmat,0);
    weight_op.symInd = find(weight_op.symTmat);
else
    weight_op.symTmat = [];
end
weight_op.p = opts.p;
if not(isempty(initialization))
    weight_op.U = initialization.U;
    weight_op.V = initialization.V;
    weight_op.sing = initialization.sing;
    eps_c = initialization.eps;
    weight_op.eps = eps_c;
    first_iterate = 0;
    X_c.U = weight_op.U;
    X_c.V = weight_op.V;
    X_c.Gam1 = diag(weight_op.sing);
    X_c.Gam2 = zeros(r_c,n.d2);
    X_c.Gam3 = zeros(n.d1,r_c);
    X_c.res_range = zeros(1,length(y));
    
else %% use uniform weights in first iteration
    weight_op.U = eye(n.d1,r_c);
    weight_op.V = eye(n.d2,r_c);
    weight_op.sing = zeros(r_c,1);
    eps_c=Inf;
    weight_op.eps = eps_c;
    first_iterate = 1;
    X_c = [];
end
X_c_newbase = zeros(r_c*(r_c+n.d1+n.d2),1);
weight_op.Wl_inv        = max(weight_op.sing,eps_c).^(2-weight_op.p);
weight_op.Wl_inv_eps    = eps_c^(2-weight_op.p); 
end