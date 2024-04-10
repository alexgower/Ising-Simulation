function single_instance_Es_by_algorithm = mem_compare(sz,flist)

     %%% 0.A PROBLEM GENERATION %%%

     % Tiling Parameters
     % sz = [8 8 8]; % nxmxk 3D lattice
     % flist = [0 0 0 0 1]; % tiling cube ratios (f21, f22, f41, f42, f6)
     fRBM = false; % non-RBM structure
     gauge_transform_W = true;

     % 0.A i) Generate W matrix instance and solution energy Esol given above tiling parameters
     % e.g. for 3d, W is a (n,m,k,6) matrix i.e. left,right,front,back,down,up couplings per vertex v in 3D cubic lattice)
     [W,Esol] = tiling(sz,flist);

     % 0.A ii) If gauge transform W then gauge transform so planted solution isn't ferromagnetic
     if gauge_transform_W
     % Generate random solution vector s
     s = -1+2*round(rand(sz));
     % Gauge transform W = W.*s.*get_v(s)
     W = W .* s .* get_V(s);
     end

     % 0.A iii) Generate initial configuration
     % Use same starting conf for both algorithms so know not initial condition difference
     % (Use cos(conf) for Pei algorithm)
     % phase_init_oim = acos(double(-1+2*round(rand(sz,'single')))); % SAME AS PEI CASE BUT WITH PHASES - CHANGE BACK
     phase_init_oim = acos(-1+2*rand(sz,'single'));
     sigma_init = cos(phase_init_oim); % COULD ALSO LET IT START AT UP/DOWN IF BETTER FOR PEI - THIS IS ONLY MODIFCATION TO PEI




     %%% 0.B ALGORITHM SETUP %%%
     % [bond-based memory, double precision]
     % i.e. bond-based memory != 2.5, plaquette-based memory = 2.5
     % i.e. double precision = 0/1 = single/double precision for W
     falgo = [2 1];
     T = 2^12; % total simulation time (he did 2^12)
     monitor = [0 1 0]; % [quiet record save] 


     % Memory algorithm parameters
     % vars = [alpha beta gamma delta zeta xini t dt];
     vars_pei = [0.799 1.040 0.836 7.068 2.201 0.683 6.192 -3.184]; % Pei optimal vars
     vars_oim_memory = [0.80489 1.0503 0.83947 7.0614 2.227 0.69678 6.4122 -3.2073]; % Alex optimal vars
     xini = 6.

     % Set up initial parameters for memory algorithms
     d = length(sz);

     conf{1} = sigma_init; % [-1,1] config
     conf{2} = vars_pei(6)*ones([sz 2*d]);; % X
     conf{3} = ones([sz 2*d]); % Y

     conf_oim{1} = phase_init_oim; % Phase config
     conf_oim{2} = vars_oim_memory(6)*ones([sz 2*d]); % X
     conf_oim{3} = ones([sz 2*d]); % Y





     %%% RUN ALGORITHMS %%%

     % 1. Run classic OIM algorithm to get the best solution energy Ebest for above instance
     [Efinal_oim_classic,final_phases_oim_classic] = classic_oim(W,phase_init_oim);

     % 2. Run OIM adapted memory algorithm to get te best solution energy Ebest for above instance
     [Ebest_oim_memory,tt_oim,step_oim,conf_oim,state_oim] = mem_oim(vars_oim_memory,falgo,Esol,W,fRBM,T,conf_oim,monitor);
     Efinal_oim_memory = cube_ising_energy(sign(cos(conf_oim{1})),W);

     % 3. Run usual Pei algorithm to get the best solution energy Ebest for above instance
     [Ebest_pei,tt,step,conf,state_pei] = mem_pei(vars_pei,falgo,Esol,W,fRBM,T,conf,monitor);
     Efinal_pei = cube_ising_energy(sign(conf{1}),W);


     %%% PRINT RAW RESULTS %%%
     disp("Energy of the solution: " + Esol)
     disp(" ")
     disp("Energy of the final instance found from OIM-Classic: " + Efinal_oim_classic)
     disp("Energy of the final | best instances found from OIM-Memory: " + Efinal_oim_memory + " | " + Ebest_oim_memory)
     disp("Energy of the final | best instances found from Pei-Memory: " + Efinal_pei + " | " + Ebest_pei)


     %% ROUNDING PROCEDURES %%
     final_conf_oim_classic = final_phases_oim_classic;
     best_conf_oim_memory = state_oim.best_conf;
     final_conf_oim_memory = conf_oim{1};
     best_conf_pei = acos(state_pei.best_conf);
     final_conf_pei = acos(conf{1});
     [~, Ebest_rounded_oim_memory] = optimal_phase_rounding(best_conf_oim_memory,W, true);
     [~, Efinal_rounded_oim_memory] = optimal_phase_rounding(final_conf_oim_memory,W, true);
     [~, Ebest_rounded_pei] = optimal_phase_rounding(best_conf_pei,W, true);
     [~, Efinal_rounded_pei] = optimal_phase_rounding(final_conf_pei,W,true);
     [~, Efinal_rounded_oim_classic] = optimal_phase_rounding(final_conf_oim_classic,convert_W_to_J(W), false);


     %%% PRINT ROUNDED RESULTS %%%
     disp(" ")
     difference = "";
     if Efinal_rounded_oim_classic~= Efinal_oim_classic
          difference = " (IMPROVEMENT)";
     end
     disp("Energy of the (rounded) final instance found from OIM-Classic: " + Efinal_rounded_oim_classic + difference)
     difference = "";
     if Ebest_rounded_oim_memory~= Ebest_oim_memory || Efinal_rounded_oim_memory~= Efinal_oim_memory
          difference = " (IMPROVEMENT)";
     end
     disp("Energy of the (rounded) final | best instances found from OIM-Memory: " + Efinal_rounded_oim_memory + " | " + Ebest_rounded_oim_memory)
     difference = "";
     if Ebest_rounded_pei~= Ebest_pei || Efinal_rounded_pei~= Efinal_pei
          difference = " (IMPROVEMENT)";
     end
     disp("Energy of the (rounded) final | best instances found from Pei-Memory: " + Efinal_rounded_pei + " | " + Ebest_rounded_pei)


     %%% Create and return single_instance_dEs_by_algorithm array %%%
     single_instance_Es_by_algorithm = [Esol, Efinal_oim_classic, Efinal_oim_memory, Ebest_oim_memory, Efinal_pei, Ebest_pei, Efinal_rounded_oim_classic, Efinal_rounded_oim_memory, Ebest_rounded_oim_memory, Efinal_rounded_pei, Ebest_rounded_pei];

end
















     %%%%% OPTIONAL LATER %%%%%

     %%% EXTRA PLOTTING %%%

     % % OIM-Memory Phase Plotter
     % oim_phase_plotter(state_oim.confs,state_oim.confs_tlist);

     % % Ferromagnetic check if applicable
     % if ~gauge_transform_W
     % found_ferromagnetic_state = isequal(optimal_ising_configuration,ones(size(optimal_ising_configuration)));
     % disp("Optimal Pei Ising Configuration is ferromagnetic: " + found_ferromagnetic_state)
     % found_ferromagnetic_state_oim = isequal(optimal_ising_configuration_oim,ones(size(optimal_ising_configuration_oim)));
     % disp("Optimal OIM Ising Configuration is ferromagnetic: " + found_ferromagnetic_state_oim)
     % end


     %%% Save generated W to file %%%
     % filename = 'W_matrix.h5';
     % if exist(filename, 'file') == 2 % The second argument 'file' ensures we're looking for files specifically
     %     delete(filename);
     % end
     % h5create(filename, '/W', size(W)); % Create the HDF5 file and dataset
     % h5write(filename, '/W', W); % Write the matrix W to the dataset
