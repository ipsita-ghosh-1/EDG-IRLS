function [prob] = dataloader_EDG(problem,oversampling)
%dataloader_EDG This function prepares a structure array 'prob' containing
%all necessary data and problem parameters (provided via 'problem') for 
% any EDG algorithm to be run.

if strcmp(problem.type,'GaussianData')
    X = rand(problem.n,problem.r);
elseif strcmp(problem.type,'GaussianDataIllcond')
    complex = false;
    [X,~] = sample_X0_lowrank(problem.n,problem.n,problem.r,...
        problem.modeX0,complex,problem.cond_nr);
    X=100.*X./norm(X);
    [problem.n,problem.r] = size(X);
% elseif strcmp(problem.type,'ScaledSGD_IllCond')
%     if not(problem.r == 3)
%         error('Example need to have rank 3.')
%     end
    [~,~,~, ~, ~, XI] = generate_EDM_ScaledSGD_mod(problem.n,problem.nr_outliers);
    X = XI;
elseif strcmp(problem.type,'USCities')
    load('./data/UScities.mat');
    X = spt(:,1:2);
    [problem.n,problem.r] = size(X);
    % disp( size(X))
elseif strcmp(problem.type,'1BPM')
    load('./data/1BPM.mat'); 
    X = prob.A';
    [problem.n,problem.r] = size(X);
    disp( size(X))
elseif strcmp(problem.type,'1AX8')
    load('./data/1AX8.mat'); 
    X = prob.A';
    [problem.n,problem.r] = size(X);
    disp( size(X))
end

P = X';
Dist = bsxfun(@plus,dot(P,P,1)',dot(P,P,1))-2*(P'*P);



if strcmp(problem.type,'SNL')

    % dist = squareform(pdist(P'));
    % D = dist.*dist;
    
    % deg_freedom = problem.n*problem.r-problem.r*(problem.r-1)/2;
    m = floor(problem.n /2)
    % E = D(1:m,1:m);
    gamma = 0.01;
    Weight=rand(m-1,m-1);
    Weight(Weight>1-gamma)=1;
    Weight(Weight<1)=0;
    Weight(Weight>0)=1;
    for i=1:m-1
        Weight(i,i)=1;
            for j=i+1:m-1
                Weight(i,j)=Weight(j,i);
            end
    end
    OmegaLin = find(Weight == 1);
    size(OmegaLin')
 
else 
    deg_freedom = problem.n*problem.r-problem.r*(problem.r-1)/2;

    m_sym = floor(oversampling*deg_freedom)
   
    if m_sym <= problem.n*(problem.n-1)/2
        OmegaLin = sort(randperm(problem.n*(problem.n-1)/2,m_sym));
        size(OmegaLin)
    end
end 



[OmegaUpperTrig,OmegaMatrix] = get_OmegaUpperTrig(problem.n,OmegaLin);

prob.n = problem.n;
prob.Phi = OmegaUpperTrig;
prob.sampled_ind = OmegaLin;
prob.Phi_all = sparse(OmegaMatrix);
Omega = find(prob.Phi);
[OmegaI,OmegaJ] = find(prob.Phi_all);
prob.OmegaTuple = [OmegaI,OmegaJ];
prob.y = [Dist(Omega);zeros(problem.n,1)];
prob.Omega = Omega;
prob.Omega_all = find(prob.Phi_all);
prob.y_all = Dist(prob.Omega_all);
prob.problem_type = 'MultidimensionalScaling';
prob.d1 = size(prob.Phi,1);
prob.d2 = size(prob.Phi,2);
prob.r = problem.r; % 
prob.Dist = Dist;
prob.Dist_sampled = prob.Dist.*prob.Phi_all;
prob.X = X;

end