function [outs_cell, success_matrix_procrustes,success_matrix,success_matrix_rankr,rel_error_PointDist,rel_error_Dist,...
    rel_error_DistObserved,times,prob_cell]= ...
    Phase_transition_parallel(problem,opts,r_list,oversampling_list, ...
    instancesize,alg_name)
success_threshold = 1e-4;
nr_algos = length(alg_name);
Xr_cell= cell(1,instancesize);
alg_name_cell = cell(1, instancesize);
prob_cell = cell(length(r_list),size(oversampling_list,2),instancesize);
outs_cell = cell(length(r_list),size(oversampling_list,2),instancesize);
% Define 4-tensor that contains all relative Frobenius norm errors with respect to the ground truth distance matrix 
rel_error_Dist = zeros(length(r_list),size(oversampling_list,2),nr_algos,instancesize);
rel_error_gram = zeros(length(r_list),size(oversampling_list,2),nr_algos,instancesize);

%Define 4-tensor that contains all relative frob error between ground truth
%and procastes alignment 
rel_error_PointDist = zeros(length(r_list),size(oversampling_list,2),nr_algos,instancesize);

% Define 4-tensor that contains all relative l2-norm errors on observed
% entries of the distance matrix
rel_error_DistObserved = zeros(length(r_list),size(oversampling_list,2),nr_algos,instancesize);
% Define matrix of success rates for ground truth recovery
success_matrix = zeros(length(r_list),size(oversampling_list,2),nr_algos);
% Define matrix of success rates for rank-r completion of Gram matrix
% complatible with measurements
success_matrix_rankr = zeros(length(r_list),size(oversampling_list,2),nr_algos);
success_matrix_procrustes = zeros(length(r_list),size(oversampling_list,2),nr_algos);
times = zeros(length(r_list),size(oversampling_list,2),nr_algos,instancesize);
for i=1: length(r_list)
    problem.r = r_list(i);
    % opts.rtilde;
    opts.r = problem.r; % ground truth rank
    opts.rtilde = problem.r; % rank estimate used by algorithms
    for j = 1:size(oversampling_list,2)
        success = zeros(instancesize,nr_algos);
        success_rankr_completion = zeros(instancesize,nr_algos);
        success_procrustes_dist = zeros(instancesize,nr_algos);
        oversampling = oversampling_list(j);
        parfor t=1:instancesize
            [prob] = dataloader_EDG(problem,oversampling);
            [Xr,outs,alg_name_out] = run_EDG_algos(prob,alg_name,opts);
            % Save results in cell arrays
            alg_name_cell{t} = alg_name_out;
            Xr_cell{t} = Xr;
            
            proxy_prob = prob;
            %%
            %remove prob.dist for large data 
            % remove X fro protein and US city
            % fields_to_keep = {'Omega','y','d1','d2','X'};
            fields_to_Remove = {'Phi','sampled_ind','Phi_all'...
             'OmegaTuple','Omega_all','y_all', 'problem_type'...
             'Dist','Dist_sampled','n'}
            proxy_prob = rmfield(proxy_prob,fields_to_Remove);
            prob_cell{i,j,t} = proxy_prob;
            % if strcmp(alg_name_out,'MatrixIRLS')
            %     outs_proxy = outs;
            %     field_to_Remove = {'meas_op'};
            %     outs_proxy = rmfield(outs_proxy,field_to_Remove);
            %     outs_cell{i,j,t} = outs_proxy;
           
            outs_cell{i,j,t} = outs;
            % end 
        %% Show statistics
            for k =1:nr_algos
                opts.rtilde
                Dist_estimated = kappainv_EDG(Xr_cell{t}{k});
                times(i,j,t,k) = outs{k}.time(end);
                rel_error_Dist(i,j,k,t) = (norm(Dist_estimated-prob.Dist,'fro')/norm(prob.Dist,'fro'));
                % rel_error_gram(i,j,k,t) = norm((Xr_cell{t}{k})-kappa_EDG(prob.Dist),'fro')./norm(kappa_EDG(prob.Dist),'fro')
                if rel_error_Dist(i,j,k,t) < success_threshold
                    success(t,k) = 1;
                end 
                rank_output = outs_cell{i,j,t}{k}.rank_out;
                % rel_error_PointDist(i,j,k,t) = sqrt(procrustes(prob.X, outs_cell{i,j,t}{k}.P))/norm(prob.X,'fro');
                % rel_error_PointDist(i,j,k,t) = sqrt(procrustes(prob.X, outs_cell{i,j,t}{k}.P));

                rel_error_DistObserved(i,j,k,t) = (norm(Dist_estimated(prob.Omega_all)-prob.Dist(prob.Omega_all),'fro')./norm(prob.Dist(prob.Omega_all),'fro'));
                if rel_error_DistObserved(i,j,k,t) < success_threshold && rank_output == prob.r
                    success_rankr_completion(t,k) = 1;
                end
                if rel_error_PointDist(i,j,k,t) < success_threshold
                    success_procrustes_dist(t,k) = 1;
                end
                if opts.verbose
                    disp(['Rank ',num2str(r_list(i)),', oversampling ',num2str(oversampling_list(j)),...
                        ', rel. rec. error for instance ',num2str(t),' for ',...
                        alg_name{k},': ',num2str(rel_error_Dist(i,j,k,t))])
                    disp(['Rank ',num2str(r_list(i)),', oversampling ',num2str(oversampling_list(j)),...
                        ', rel. rec. error gram for instance ',num2str(t),' for ',...
                        alg_name{k},': ',num2str(rel_error_gram(i,j,k,t))])
                    disp(['Rank ',num2str(r_list(i)),', oversampling ',num2str(oversampling_list(j)),...
                        ', rel. rec. error on observed entries for instance ',num2str(t),' for ',...
                        alg_name{k},': ',num2str(rel_error_DistObserved(i,j,k,t))])
                    disp(['Rank ',num2str(r_list(i)),', oversampling ',num2str(oversampling_list(j)),...
                        ', rel. error for pracastes alginment',num2str(t),' for ',...
                        alg_name{k},': ',num2str(rel_error_PointDist(i,j,k,t))])
                    % disp(['Relative Frob. error for ',alg_name{l},' :', num2str(norm(kappainv_EDG(Xr{l})-prob.Dist,'fro')./norm(prob.Dist,'fro'))])
                end
                % success_matrix(i, j, k) = sum(success(:,k)) / instancesize;
            end
        end
        % Xr_matrix = cell2mat(Xr_cell);
        for k=1:nr_algos
            success_matrix(i, j, k) = sum(success(:,k)) / instancesize;
            success_matrix_rankr(i, j, k) = sum(success_rankr_completion(:,k)) / instancesize;
            success_matrix_procrustes(i,j,k) = sum(success_procrustes_dist(:,k))/ instancesize;
        end
           
    end
end 
delete(gcp);
