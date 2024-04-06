%%%%% tune memory parameters based on simplex descent in parameter space
%%%%% note: uses parallel processes
clear;
addpath(genpath('..')); 

%%%%% PARAMETERS %%%%%
% varsmem = [0.799 1.040 0.836 7.068 2.201 0.683 6.192 -3.184]; % original best parameters for memory
varsmem = [0.80489      1.0503     0.83947      7.0614       2.227     0.69678      6.4122     -3.2073] % Alex
% varlb = [0 0 0 0 0 0 4 -5]; % parameter lower bounds
varlb = [0 0 0 0 0 0 40 -50]; % parameter lower bounds
% varub = [2 2 2 1 5 1 10 -2]; % parameter upper bounds
varub = [20 20 20 10 50 10 100 -20]; % parameter upper bounds
indices = 1:8; % tune over all 8 indices

% Algorithm
% [bond-based memory, double precision]
% i.e. bond-based memory != 2.5, plaquette-based memory = 2.5
% i.e. double precision = 0/1 = single/double precision for W
falgo = [2 2]; % first 2: memory; second 2: double precision % CHECK maybe 2=1=double precision

sz = [6 6 6];
flist = [0 0 0 0 1]; % code 11 = gaussian distribution % CHECK actually 11 = fully frustated?
fRBM = false; % RBM instances
oim = 1; % Use oim version of memory algorithm

steps = 100; % simplex steps
runs = 40; % random instances per step
T = 2^10; % simulation time

conf = [];
monitor = [0 0 0];

fSA = 0; % method of parameter searching; false=simplex descent; true=simulated annealing
fAPX = 0; % leave 0 if solution unknown
gap = 0; % leave 0 if solution unknown

var = tuning(varsmem,varlb,varub,indices,falgo,sz,flist,fRBM,steps,runs,T,conf,monitor,fSA,fAPX,gap,oim);