clear;
addpath(genpath('..'));

% Tiling Parameters
sz = [8 8 8]; % nxmxk 3D lattice
flist = [1 1 1 1 1]; % tiling cube ratios (f21, f22, f41, f42, f6)
fRBM = false; % non-RBM structure

% Algorithm
% [bond-based memory, double precision]
% i.e. bond-based memory != 2.5, plaquette-based memory = 2.5
% i.e. double precision = 0/1 = single/double precision for W
falgo = [2 1];
gauge_transform_W = true;

% ORIGINAL VERSION OPTIMAL TIMESCALE PARAMETERS FOR PEI
alpha = 0.799; 
beta = 1.040; 
gamma = 0.836; 
delta = 7.068;
zeta = 2.201;
xini = 0.683; % memory initial value
t = 6.192; % time before restart - as a power of 2 i.e. will be 2^t
dt = -3.184; % maximum stepsize - as a power of 2 i.e. will be 2^dt

% vars = [alpha beta gamma delta zeta xini t dt]; % Pei optimal vars
vars = [0.80489      1.0503     0.83947      7.0614       2.227     0.69678      6.4122     -3.2073] % alex optimal vars

% Use same starting conf for both algorithms so know not initial condition difference
% (Use cos(conf) for Pei algorithm)
% phase_init_oim = acos(double(-1+2*round(rand(sz,'single')))); % SAME AS PEI CASE BUT WITH PHASES - CHANGE BACK
phase_init_oim = acos(-1+2*rand(sz,'single'));
sigma_init = cos(phase_init_oim) % COULD ALSO LET IT START AT UP/DOWN IF BETTER FOR PEI
d = length(sz);
X = xini*ones([sz 2*d]);
Y = ones([sz 2*d]);

conf{1} = sigma_init;
conf{2} = X;
conf{3} = Y;
conf_oim{1} = phase_init_oim;
conf_oim{2} = X;
conf_oim{3} = Y;


% Simulation
T = 2^12; % total simulation time (he did 2^12)
monitor = [0 1 0]; % [quiet record save] 

% Generate W matrix instance and solution energy Esol given above tiling parameters
% e.g. for 3d, W is a (n,m,k,6) matrix i.e. left,right,front,back,down,up couplings per vertex v in 3D cubic lattice)
[W,Esol] = tiling(sz,flist);

% If gauge transform W
if gauge_transform_W
    % Generate random solution vector s
    s = -1+2*round(rand(sz));
    % Gauge transform W = W.*s.*get_v(s)
    W = W .* s .* get_V(s);
end


% Run memory algorithm to get te best solution energy Ebest for above instance
% Ebest = mem(vars,falgo,Esol,W,fRBM,T,[],monitor);
[Ebest_oim,tt_oim,step_oim,conf_oim,state_oim] = mem_oim(vars,falgo,Esol,W,fRBM,T,conf_oim,monitor);

[Ebest,tt,step,conf,state] = mem(vars,falgo,Esol,W,fRBM,T,conf,monitor);

% Print results
disp("Energy of the solution: " + Esol)
disp("Energy of the best instance found from Pei: " + Ebest)
disp("Energy of the best instance found from OIM: " + Ebest_oim)



% UNCLEAR WHETHER MY ROUNDING PROCEDURES MAKE A DIFFERENCE MUCH
best_conf = acos(state.best_conf);
best_conf_oim = state_oim.best_conf;
[optimal_ising_configuration, optimal_ising_energy] = optimal_phase_rounding(best_conf,W);
disp("Optimal Plane Cut Energy Pei: " + optimal_ising_energy)
[optimal_ising_configuration_oim, optimal_ising_energy_oim] = optimal_phase_rounding(best_conf_oim,W);
disp("Optimal Plane Cut Energy OIM: " + optimal_ising_energy_oim)

if ~gauge_transform_W
    found_ferromagnetic_state = isequal(optimal_ising_configuration,ones(size(optimal_ising_configuration)));
    disp("Optimal Pei Ising Configuration is ferromagnetic: " + found_ferromagnetic_state)
    found_ferromagnetic_state_oim = isequal(optimal_ising_configuration_oim,ones(size(optimal_ising_configuration_oim)));
    disp("Optimal OIM Ising Configuration is ferromagnetic: " + found_ferromagnetic_state_oim)
end

% ALEX - save generated W to file
filename = 'W_matrix.h5';
if exist(filename, 'file') == 2 % The second argument 'file' ensures we're looking for files specifically
    delete(filename);
end
h5create(filename, '/W', size(W)); % Create the HDF5 file and dataset
h5write(filename, '/W', W); % Write the matrix W to the dataset


%%% PHASE GRAPH PLOTTING %%%
oim_phase_plotter(state_oim.confs,state_oim.confs_tlist);