
% This determines which algorithm you're using. Set it to
% 0 for P_omega, 2 for RstarR_omega 
descent_mode = 2;

%there used to be another initialization, not included due to poor performance
initialization_type = 0;

%set sampling_mode = 0 for Bernoulli sampling of the entries, sampling_mode
%= 1 for a fixed value m number of samples
sampling_mode = 1;

% Select the names of the data types
% data_array = ["./data/1k.off","./data/cow.off","swiss","cities","ellipse"];
data_array = ["cities"];
% Choose the sampling rates
% rate_array = [0.15,0.10,.07,.05,.03,.02,.01]
rate_array = [0.75];%,0.15,0.10,.07,.05,.03,.02,.01];
%rate_array = [.9,.7,.5,.3]
% Number of times to run each dataset (with a new randomized subset of
% points
num_trials = 1;
% Max iterations stopping condition for the algorithms
num_iter = 300;
% relative difference threshold 
norm_thresh = 10^-7;
rel_thresh = 10^-7;
% 3d result array
results = zeros(size(data_array,2),size(rate_array,2),num_trials);
i =1;
%load data
% for i=1:length(data_array)
%     data_loc = data_array(i)
%     % if data_loc == "ellipse"
%     %     t = linspace(0,2*pi,600);
%     %     x = 100*cos(t);
%     %     y = 0.01*sin(t);
%     %     P = [x' y'];
%     % elseif data_loc == "./data/cow.off" %datapoints P and mesh trg
%     %     [P, trg] = ReadOFF(data_loc,'1');
%     % elseif data_loc == "./data/1k.off"
%     %     [P, trg] = ReadOFF(data_loc,'1'); %datapoints P and mesh trg
%     % elseif data_loc == "swiss"
%     %     load('./data/ptswiss.mat');
%     %     P = pt;
%     % elseif data_loc == "cities"
%     %     load('./data/UScities.mat');
%     %     P = spt(:,1:2); %can make the last arg : to get altitude component
%     % else
%     %     disp 'You should write something here to load your dataset!'
%     %     break
%     if data_loc == "cities"
load('./data/UScities.mat');
P = spt(:,1:2); 

    % end

%number of datapoints
%n_used = 100;
%P = P(1:n_used,:);
n = size(P,1);
%dimension of datapoints
d = size(P,2);
%number of squared distances (useful for one of the sampling modes)
L = n*(n-1)/2;
%make sure that the datapoints have zero mean for reconstruction
P = P - sum(P,1)/n; %centering the data

%build the true gram and distance matrices
G = P*P';%gram matrix
D = ones(n,1)*diag(G)'+diag(G)*ones(1,n)-2*G;%distance matrix
for q=1:n
    D(q,q) = 0.0;
end
    
L = n*(n-1)/2;
G = -1/2*(eye(n) - 1/n*ones(n))*D*(eye(n)-1/n*ones(n));

%cycle through rates
for j=1:size(rate_array,2)
rate = rate_array(j)
    for k=1:num_trials
        k
        [samples,Weight] = Construct_Samples(n,n,L,rate,sampling_mode);

        if descent_mode == 0
            if initialization_type == 0
                  D_hat = P_omega(D,samples);% this is the standard initialization
            else
                disp 'Choose 0 for initialization' 
                break
            end
            [D_approx,~,~] = distgeo_rgrad_Pomega(D_hat,samples,d,num_iter,D,rate,1,1,norm_thresh);
            results(i,j,k) = norm(D_approx-D,'fro')/norm(D,'fro');

        elseif descent_mode == 2

            disp('R^*_Omega R_Omega Algorithm Running') %you can comment this line out, i use it just for sanity purposes
            if initialization_type == 0
                G_hat = RstarR_omega(G,samples); %standard initialization
            else
                disp 'Choose either 0 or 1 for initialization'
                break
            end
            [G_approx,P_approx] = distgeo_rgrad_RstarR(G_hat,samples,d,num_iter,G,rate,1,0,norm_thresh,rel_thresh);
            results(i,j,k) = norm(G_approx-G,'fro')/norm(G,'fro'); 
            %norm(G_approx-G,'fro')/norm(G,'fro')
        else
            disp('Choose either 0 or 2')
        end
    end
   
end
disp(results)
% end 

%ViewMesh(P_approx,trg) this is the function you'd use to visualize the
%sphere or the cow. The other 3 datasets you can just scatter plot          disp('P_Omega Algorithm Running')
                % if initialization_type == 0
    
