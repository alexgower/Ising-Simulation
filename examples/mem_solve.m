clear;
addpath(genpath('..'));

% Tiling Parameters
sz = [4 4 4]; % 6x6x6 3D lattice
flist = [1 1 1 1 1]; % tiling cube ratios (f21, f22, f41, f42, f6)
% flist = [1 1 1 1 1]; % tiling cube ratios (f21, f22, f41, f42, f6)
fRBM = false; % non-RBM structure

% Algorithm
% [bond-based memory, double precision]
% i.e. bond-based memory != 2.5, plaque-based memory = 2.5 - MAYBE
% i.e. double precision = 0/1 = single/double precision for W
falgo = [2 1];

alpha = 0.799; 
beta = 1.040; 
gamma = 0.836; 
delta = 7.068; 
zeta = 2.201; % time-scale parameters

xini = 0.683; % memory initial value
conf = []; % random spin initialization
t = 6.192; % time before restart - as a power of 2 i.e. will be 2^t
dt = -3.184; % maximum stepsize - as a power of 2 i.e. will be 2^dt
vars = [alpha beta gamma delta zeta xini t dt];

% Simulation
T = 2^12; % total simulation time
monitor = [0 1 1]; % [quiet record save] 

% Generate W matrix instance and solution energy Esol given above tiling parameters
% e.g. for 3d, W is a (n,m,k,12) matrix i.e. left,right,front,back,down,up couplings per vertex v in 3D cubic lattice)
[W,Esol] = tiling(sz,flist);

% Run memory algorithm to get the best solution energy Ebest for above instance
% Ebest = mem(vars,falgo,Esol,W,fRBM,T,[],monitor);
[Ebest,tt,step,conf,state] = mem(vars,falgo,Esol,W,fRBM,T,[],monitor);

% Print results
disp("Energy of the instance: " + Esol)
disp("Energy of the best solution: " + Ebest)

% ALEX - save generated W to file
% h5create('W_matrix.h5', '/W', size(W)); % Create the HDF5 file and dataset
% h5write('W_matrix.h5', '/W', W); % Write the matrix W to the dataset
