function [optimal_ising_configuration, optimal_ising_energy] = optimal_phase_rounding(phases, W)
    d = length(size(phases));
    
    % Generate all cutting planes
    epsilon = 0; % 0.001; % Small offset
    planes_normal_phases = phases + epsilon + pi/2; % Vectorized addition

    optimal_ising_configuration = sign(cos(phases)); % Initialize with z-rounded configuration
    % Calculate z-rounded energy
    sigma_0 = optimal_ising_configuration;
    Sigma_0 = get_V(sigma_0);
    optimal_ising_energy = double(sum(sigma_0.*sum(Sigma_0.*W,d+1),'all')/2); 

    for i = 1:length(planes_normal_phases)
        plane_normal_phase = planes_normal_phases(i);
        
        % Redefine 'z' direction based on the plane's normal
        new_phases = phases - plane_normal_phase;
        this_ising_configuration = sign(cos(new_phases));

        sigma_0 = this_ising_configuration;
        Sigma_0 = get_V(sigma_0);
        
        this_ising_energy = double(sum(sigma_0.*sum(Sigma_0.*W,d+1),'all')/2);

        % Update if this configuration has the best energy so far
        if this_ising_energy > optimal_ising_energy
            optimal_ising_configuration = this_ising_configuration;
            optimal_ising_energy = this_ising_energy;
        end
    end

    % DO THE SAME WITH PLANES FROM PARALLEL PHASES
    planes_parallel_phases = phases + epsilon; % Vectorized addition

    for i = 1:length(planes_parallel_phases)
        plane_parallel_phase = planes_parallel_phases(i);
        
        % Redefine 'z' direction based on the plane's normal
        new_phases = phases - plane_parallel_phase;
        this_ising_configuration = sign(cos(new_phases));

        sigma_0 = this_ising_configuration;
        Sigma_0 = get_V(sigma_0);
        
        this_ising_energy = double(sum(sigma_0.*sum(Sigma_0.*W,d+1),'all')/2);

        % Update if this configuration has the best energy so far
        if this_ising_energy > optimal_ising_energy
            optimal_ising_configuration = this_ising_configuration;
            optimal_ising_energy = this_ising_energy;
        end
    end
end