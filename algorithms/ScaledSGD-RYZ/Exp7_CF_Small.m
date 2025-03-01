%% Experiment 7: Small-scale Item-item Collaborative Filtering
clear; addpath('Functions')
loader = load('Data/CF_1M.mat');
spdata = loader.spdata;
d = loader.n_movie;
r = 3;
epochs = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ScaledSGD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(1); learning_rate = 1e3;
[~, fscsgd, aucscsgd] = bpr_scaledsgd(spdata, d, r, epochs, learning_rate, true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SGD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(1); learning_rate = 5e-2;
[~, fsgd, aucsgd]     = bpr_scaledsgd(spdata, d, r, epochs, learning_rate, false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NP-Maximum %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(1); learning_rate = 1e-1; 
[~, np_maximum] = bpr_npmaximum(spdata, d, 100, learning_rate);

%%%%%%%%%%%%%%%%%%%%%%%% Plot ScaledSGD vs SGD %%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotCF(fscsgd,fsgd,aucscsgd,aucsgd,np_maximum)
save('Data/Results/EXP7.mat','fscsgd','fsgd','aucscsgd','aucsgd','np_maximum')