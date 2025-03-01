function [Xr,outs,alg_name] = run_EDG_algos(prob,alg_name,opts_custom)
% This function runs different algorithms (indicated by 'alg_name')
% for a given Euclidean distance geometry problem with problem definition
% encoded by structure array 'prob'.

% =========================================================================
% Parameters
% ----------
% prob: structure array. Encodes problem definition model and relevant
%       parameters
% r:    integer. Target rank used by reconstruction algorithms (if
%       applicable)
% alg_name:     (1 x nr_algos) cell of character strings. Indicating which 
%               algorithm to use.
% opts_custom:  struct. Options to override standard options of any 
%               algorithm.
% Returns
% ----------
% Xr:   (1 x nr_algos_total) cell. The l-th cell contains algorithmic
%       iterates of l-th algorithm (if opts.saveiterates == 1) or
%       the last iterate (if opts.saveiterates == 0).
% outs: (1 x nr_algos_total) cell. The l-th cell contains additional
%       information about the progress of the l-th algorithm.
% alg_name: (1 x nr_algos_total) cell of character strings. Contains name
%       of used algorithms, including reflecting potentially additional 
%       parameter choices (and thus, might differ from input 'alg_name').
% =========================================================================

%%% prepare multiple algorithm instances of IRLS variants if multiple
%%% values for quasi-norm parameter p are provided
opts_extendnames_1 = {'p'};
alg_name_extend = {'MatrixIRLS'};
[opts_new,alg_name1,alg_name_out1] = extend_algopts(opts_custom,opts_extendnames_1,...
    alg_name,alg_name_extend);

opts_extendnames_2 = {'type_mean','mode_eps'};
[opts_new,alg_name2,alg_name_out] = extend_algopts(opts_custom,opts_extendnames_2,...
    alg_name_out1,alg_name_extend,opts_new);
alg_name = alg_name_out;

nr_algos = length(alg_name);
outs  = cell(1,nr_algos);
Xr    = cell(1,nr_algos);
opts  = cell(1,nr_algos);
lambda = cell(1,nr_algos);
%params = cell(1,nr_algos);

%%% run algorithm (with default options, priority on option in opts_new if
%%% applicable)
for l=1:nr_algos
    % nr_algos
    opts{l} = getDefaultOpts_IRLS;
    opts{l} = setExtraOpts(opts{l},opts_new{l});
    opts
    if contains(alg_name{l},'MatrixIRLS')        
        prob_c = prob;
        prob_c.problem_type = 'MultidimensionalScaling';
        [X_c,outs{l}] = ...
        MatrixIRLS(prob_c,opts{l}.lambda,opts{l});
        % outs{l}.P = get_pointmatrix_from_iterate(X_c,opts{l}.r);
        % %%
        % Xr{l} = get_densemat_from_iterate(X_c,outs{l}.meas_op);
        rankr_threshold = 1e-6;
        if X_c.res_range < rankr_threshold
            if (norm(X_c.Gam1) < rankr_threshold) && (norm(X_c.Gam2) < rankr_threshold)
                outs{l}.rank_out = outs{l}.r_greatereps(end);
            else
                outs{l}.rank_out = 2*outs{l}.r_greatereps(end);
            end
        else
            outs{l}.rank_out = size(prob.Dist,1);
        end
        %outs{l}.rank_out = outs{l}.r_greatereps(end);
%%
        if isfield(opts{l},'saveiterates') && opts{l}.saveiterates == 1
            outs{l}.P = cell(outs{l}.N,1);
            Xr{l} = cell(outs{l}.N,1);
            for h=1:outs{l}.N 
                % Xr{l,h}=outs{l}.X{h};
                outs{l}.P{h} = get_pointmatrix_from_iterate(outs{l}.X{h},opts{l}.r);
                % 
                Xr{l}{h} = get_densemat_from_iterate(outs{l}.X{h},outs{l}.meas_op);
            end
        else
        %      Xr{l}{1} = X_c;
                outs{l}.P = get_pointmatrix_from_iterate(X_c,opts{l}.r);
        %%
                Xr{l} = get_densemat_from_iterate(X_c,outs{l}.meas_op);
        end
    elseif contains(alg_name{l},'AL_BurerMonteiro')
        opts{l}.maxit = opts{l}.N0_firstorder;
        opts{l}.printenergy = opts{l}.verbose; 
        lsopts.maxit = 20;
        lsopts.xtol = 1e-8;
        lsopts.gtol = 1e-8; 
        lsopts.ftol = 1e-10; 
        lsopts.alpha  = 1e-3;
        lsopts.rho  = 1e-4; 
        lsopts.sigma  = 0.1; 
        lsopts.eta  = 0.8; 
        % [iterate_alignment, ipm, outs{l}] = alternating_completion(prob.Dist,prob.Phi_all,opts{l},lsopts);
        [iterate_alignment, ipm, outs{l}] = alternating_completion(prob.Dist_sampled,prob.Phi_all,opts{l},lsopts);

        % Xr{l} = GCor;
    %     outs{l}.DistRec = ipm;
        outs{l}.rank_out = opts{l}.rtilde;
        if isfield(opts{l}, 'saveiterates') && opts{l}.saveiterates == 1
            % Initialize cell array to store P iterates
            % P_iterates = cell(opts{l}.maxit, 1);
            outs{l}.P = cell(size(iterate_alignment,2),1);
            Xr{l} = cell(size(iterate_alignment,2),1);
            for h = 1:size(iterate_alignment,2)
                % Save P in each iteration
                outs{l}.P{h} = iterate_alignment{h}; % Assuming iterate_alignment{h} contains the P matrix for iteration h
                % Save IPM (inner-product matrix)
                Xr{l}{h} = ipm{h}; % Assuming ipm{h} contains the inner-product matrix for iteration h
            end
        else
            % Save the final P matrix
            outs{l}.P = iterate_alignment;
            % Save the final inner-product matrix
            Xr{l} = ipm;
        end
            
   elseif contains(alg_name{l},'ReiEDG')
        G_kappa = kappa_EDG(prob.Dist);
        disp('R^*_Omega R_Omega Algorithm Running') %you can comment this line out, i use it just for sanity purposes
        G_hat = RstarR_omega(G_kappa,prob.OmegaTuple); %standard initialization
        result_print = 1;
        diff_threshing = 0;
        % opts{l}.saveiterates = 1;
        norm_thresh = opts{l}.tol;
        rel_thresh = opts{l}.tol;
        N0 = opts{l}.N0_firstorder;
        G = []; % Ground truth gram matrix
        [G_approx,P_approx,iterate_alignment,time] = distgeo_rgrad_RstarR(G_hat,prob.OmegaTuple,opts{l}.rtilde,N0, ...
            prob.Dist_sampled,G, ...
            result_print,diff_threshing,norm_thresh,rel_thresh);
        % Xr{l} =  G_approx;
        % outs{l}.P_approx = P_approx;
        outs{l}.rank_out = size(P_approx,2);
        outs{l}.time = time;
        % outs{l}.P = iterate_alignment
        if isfield(opts{l},'saveiterates') && opts{l}.saveiterates == 1
            outs{l}.P = cell(size(iterate_alignment,2),1);
            Xr{l} = cell(size(iterate_alignment,2),1); %opts{l}.N0_firstorder
            for h=1:size(iterate_alignment,2)
                % Xr{l,h}=outs{l}.X{h};
                outs{l}.P{h} = iterate_alignment{h};
                % 
                Xr{l}{h} = G_approx{h};
            end
        else
            outs{l}.P = iterate_alignment;
            Xr{l} =  G_approx;
        end
    elseif contains(alg_name{l},'ScaledSGD')
        epochs = opts{l}.N0_firstorder; %5*opts{l}.N0_firstorder;
        doscale = true;
        lossfun = 'EDM';
        [P_c, outs{l}] = scaledsgd(prob.Dist_sampled, opts{l}.rtilde, epochs, opts{l}.learning_rate, lossfun, doscale,1);
        % outs{l}.P = P_c;
        outs{l}.rank_out = size(P_c,2);
        % Xr{l} = P_c*P_c';
        if isfield(opts{l},'saveiterates') && opts{l}.saveiterates == 1
            outs{l}.P = cell(size(P_c,2),1);
            Xr{l} = cell(size(P_c,2),1); %epochs
            for h=1:size(P_c,2)
                % Xr{l,h}=outs{l}.X{h};
                outs{l}.P{h} = P_c{h};
                % 
                Xr{l}{h} = outs{l}.P{h} *outs{l}.P{h}';
            end
        else
            outs{l}.P = P_c;
            Xr{l} =  P_c*P_c';
        end
    end
    outs{l}.opts=opts{l};
    if opts{l}.verbose > 0
        disp(['Time elapsed for ',alg_name{l},': ',num2str(outs{l}.time(end))])
    end
end
    

end

function opts = setExtraOpts(opts,opts_new)
    if ~isempty(opts_new)
        optNames = fieldnames(opts_new);
        for i = 1:length(optNames)
            optName = optNames{i};
            opts.(optName) = opts_new.(optName);
        end
    end

end

