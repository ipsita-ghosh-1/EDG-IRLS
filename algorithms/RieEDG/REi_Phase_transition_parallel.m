function [outs_cell, success, success_matrix]= ...
    REi_Phase_transition_parallel(problem,opts,r_list,oversampling_list, ...
    instancesize,alg_name)
success_threshold = 1e-4;
nr_algos = length(alg_name);
Xr_cell= cell(1,instancesize);
alg_name_cell = cell(1, instancesize);

success = zeros(instancesize,nr_algos);

success_matrix = zeros(length(r_list),size(oversampling_list,2),nr_algos);

for i=1: length(r_list)
    problem.r = r_list(i);
    opts.r = problem.r; % ground truth rank
    opts.rtilde = problem.r; % rank estimate used by algorithms
    for j = 1:size(oversampling_list,2)
        oversampling = oversampling_list(j);
        parfor t=1:instancesize
            [prob] = dataloader_EDG(problem,oversampling);
            [Xr,outs,alg_name_out] = run_EDG_algos(prob,alg_name,opts);
            % Save results in cell arrays
            alg_name_cell{t} = alg_name_out
            Xr_cell{t} = Xr;

            outs_cell{t} = outs;
                
        %% Show statistics
            for k =1:nr_algos

                if contains(alg_name{k},'ReiEDG') 

                    error = norm(Xr_cell{t}{k} - G,'fro')/norm(G,'fro'); 
                else error = (norm(kappainv_EDG(Xr_cell{t}{k})-prob.Dist,'fro')./norm(prob.Dist,'fro'))
                if error < success_threshold
                    success(t,k) = 1;

                end 
                % success_matrix(i, j, k) = sum(success(:,k)) / instancesize;
            end
        end
        % Xr_matrix = cell2mat(Xr_cell);
        for k=1:nr_algos
            success_matrix(i, j, k) = sum(success(:,k)) / instancesize;
        end
           
    end
end 
delete(gcp);
