clear;
addpath(genpath('..')); 


%%% PROBLEM SET GENERATION %%%
sz = [8 8 8]; % nxmxk 3D lattice
flist = [0 0 0 0 1]; % tiling cube ratios (f21, f22, f41, f42, f6)
number_of_problems = 50;
parallel_processors = 50;

Es_by_algorithm = zeros(number_of_problems, 11);

if parallel_processors
     parpool(parallel_processors);
     parfor i=1:number_of_problems
          % single_instance_Es_by_algorithm = [Esol, Efinal_oim_classic, Efinal_oim_memory, Ebest_oim_memory, Efinal_pei, Ebest_pei, Efinal_rounded_oim_classic, Efinal_rounded_oim_memory, Ebest_rounded_oim_memory, Efinal_rounded_pei, Ebest_rounded_pei];
          single_instance_Es_by_algorithm = mem_compare(sz, flist);
          Es_by_algorithm(i, :) = single_instance_Es_by_algorithm;
     end
else
     for i=1:number_of_problems
          % single_instance_Es_by_algorithm = [Esol, Efinal_oim_classic, Efinal_oim_memory, Ebest_oim_memory, Efinal_pei, Ebest_pei, Efinal_rounded_oim_classic, Efinal_rounded_oim_memory, Ebest_rounded_oim_memory, Efinal_rounded_pei, Ebest_rounded_pei];
          single_instance_Es_by_algorithm = mem_compare(sz, flist);
          Es_by_algorithm(i, :) = single_instance_Es_by_algorithm;
     end
end

%%% SAVE RESULTS %%% 

% Save the results as .mat
save('results.mat', 'Es_by_algorithm');

% Save the results as table in .csv
Problem = (1:number_of_problems)';
Esol = Es_by_algorithm(:, 1);
Efinal_oim_classic = Es_by_algorithm(:, 2);
Efinal_oim_memory = Es_by_algorithm(:, 3);
Ebest_oim_memory = Es_by_algorithm(:, 4);
Efinal_pei = Es_by_algorithm(:, 5);
Ebest_pei = Es_by_algorithm(:, 6);
Efinal_rounded_oim_classic = Es_by_algorithm(:, 7);
Efinal_rounded_oim_memory = Es_by_algorithm(:, 8);
Ebest_rounded_oim_memory = Es_by_algorithm(:, 9);
Efinal_rounded_pei = Es_by_algorithm(:, 10);
Ebest_rounded_pei = Es_by_algorithm(:, 11);
T = table(Esol, Efinal_oim_classic, Efinal_oim_memory, Ebest_oim_memory, Efinal_pei, Ebest_pei, Efinal_rounded_oim_classic, Efinal_rounded_oim_memory, Ebest_rounded_oim_memory, Efinal_rounded_pei, Ebest_rounded_pei);
writetable(T,'results.csv','Delimiter',',')  



%%% PLOT RESULTS ON A HISTOGRAM %%%
dEs_oim_classic = -Es_by_algorithm(:, 2) + Es_by_algorithm(:, 1);
dEs_oim_memory = -Es_by_algorithm(:, 3) + Es_by_algorithm(:, 1);


% Make histogram with dE bins for each algorithm
% Create the histogram for the baseline
binEdges = linspace(min([dEs_oim_classic; dEs_oim_memory]), max([dEs_oim_classic; dEs_oim_memory]), 30);
histogram(dEs_oim_classic, 'BinEdges', binEdges, 'FaceColor', 'b', 'FaceAlpha', 0.5);
hold on;
histogram(dEs_oim_memory, 'BinEdges', binEdges, 'FaceColor', 'r', 'FaceAlpha', 0.5);

xlabel('dE from Ground State Energy');
ylabel('Count');
title('OIM-Classic v OIM-Memory for Cubic Planted Solution Problems');
legend('OIM-Classic', 'OIM-Memory');
hold off;


