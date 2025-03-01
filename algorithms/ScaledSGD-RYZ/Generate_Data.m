%% Experiment 1 and Experiment 4: Synthetic Data (30 x 30 Symmetric Matrix with Rank 3)
clear; addpath('Functions')
rng(1)   % Random seed
n = 30;  % Size of matrix
r = 3;   % Rank
[MW, MI] = generate_MAT(n,r);
filename = ['Data/MAT_', num2str(n),'.mat'];
save(filename,'MW','MI','r')

%% Experiment 2: Euclidean Distance Matrix (30 Sample Points)
clear; addpath('Functions')
rng(1)    % Random seed
n = 30;   % Number of points
[DW, MW, XW, DI, MI, XI] = generate_EDM(n);
filename = ['Data/EDM_', num2str(n),'.mat'];
save(filename,'DW','MW','XW','DI','MI','XI')

%% Experiment 3: Huge-scale Collaborative Filtering (MovieLens25M, 100M sample)
clear; addpath('Functions')
rng(1)    % Random seed
% Warning: to speed up the sampling process, this code requires 32GB of RAM
% and 10 core of CPU, as well as the MATLAB Parallel Computing Toolbox.
% If the MATLAB Parallel Computing Toolbox is not available on your
% platform, please set UseParall = false.
UseParall = true;
million = 1e6;
if ~isfile('Data/ml-25m/ratings.csv')
    url = 'https://files.grouplens.org/datasets/movielens/ml-25m.zip';
    unzip(url, 'Data');
end
data=csvread('Data/ml-25m/ratings.csv',1,0);
train_size = 100*million;
test_size = 10*million;
[spdata, n_movie] = generate_CF_LR(data, train_size, test_size, UseParall);
filename = 'Data/CF_100M.mat';
save(filename, 'spdata', 'n_movie', '-v7.3')

%% Experiment 5 and Experiment 6: Synthetic Data Noisy Case (30 x 30 Symmetric Matrix with Rank 3)
clear; addpath('Functions')
addpath('Functions')
rng(1)    % Random seed
n = 30;   % Size of matrix
r = 3;    % Rank
SNR = 15; % Signal to noise ratio
[MW, MI] = generate_MAT_Noise(n, r, SNR);
filename = ['Data/MAT_Noise_', num2str(n),'.mat'];
save(filename,'MW','MI','r','SNR')

%% Experiment 7: Small-scale Collaborative Filtering (MovieLens Latest Small, 1M sample)
clear; addpath('Functions')
rng(1)    % Random seed
million = 1e6;
if ~isfile('Data/ml-latest-small/ratings.csv')
    url = 'https://files.grouplens.org/datasets/movielens/ml-latest-small.zip';
    unzip(url, 'Data');
end
data=csvread('Data/ml-latest-small/ratings.csv',1,0);
train_size = 1*million;
test_size = 0.1*million;
[spdata, n_movie] = generate_CF_SM(data, train_size, test_size);
filename = 'Data/CF_1M.mat';
save(filename, 'spdata', 'n_movie')

%% Experiment 8: Medium-scale Collaborative Filtering (MovieLens Latest Full, 10M sample)
clear; addpath('Functions')
rng(1)    % Random seed
% Warning: this code requires the MATLAB Parallel Computing Toolbox.
% If the MATLAB Parallel Computing Toolbox is not available on your
% platform, please set UseParall = false.
UseParall = true;
million = 1e6;
if ~isfile('Data/ml-latest/ratings.csv')
    url = 'https://files.grouplens.org/datasets/movielens/ml-latest.zip';
    unzip(url, 'Data');
end
data=csvread('Data/ml-latest/ratings.csv',1,0);
train_size = 10*million;
test_size = 1*million;
[spdata, n_movie] = generate_CF_LR(data, train_size, test_size, UseParall);
filename = 'Data/CF_10M.mat';
save(filename, 'spdata', 'n_movie')

%% Experiment 9: Large-scale Collaborative Filtering (MovieLens Latest Large-scale, 30M sample)
clear; addpath('Functions')
rng(1)    % Random seed
% Warning: this code requires the MATLAB Parallel Computing Toolbox.
% If the MATLAB Parallel Computing Toolbox is not available on your
% platform, please set UseParall = false.
UseParall = true;
million = 1e6;
if ~isfile('Data/ml-latest/ratings.csv')
    url = 'https://files.grouplens.org/datasets/movielens/ml-latest.zip';
    unzip(url, 'Data');
end
data=csvread('Data/ml-latest/ratings.csv',1,0);
train_size = 30*million;
test_size = 3*million;
[spdata, n_movie] = generate_CF_LR(data, train_size, test_size, UseParall);
filename = 'Data/CF_30M.mat';
save(filename, 'spdata', 'n_movie')
