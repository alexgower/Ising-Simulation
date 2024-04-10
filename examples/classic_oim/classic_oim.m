function [E, final_phases] = classic_oim(W,conf_init)

    %%% CONVERSIONS %%%
    J = convert_W_to_J(W);
    n = size(J,1);
    u0 = reshape(conf_init, [1,n]);

    %%%% Define ODE parameters %%%
    p{1} = J;
    tSpan = [0 40]; 

    %%%% SOLVE ODE %%%
    [t, u] = ode45(@(t, u) oim_dynamics(p,t,u), tSpan, u0);

    %%% Do energy calculation on final phases %%%
    final_phases = u(end,:);
    sigmas = sign(cos(final_phases));
    E = general_ising_energy(sigmas,J);

    %%% PLOT %%%
    % plot(t,u);
    % xlabel('Time');
    % ylabel('Oscillator Phases');
    % title('Classic OIM Dynamics');


end


%%%%

function K = K(t)
    K = 2*(1.0 + 0.15*t);

end

function K_s = K_s(t)
    if t < 4.0
        K_s = 0.0;
    else
        K_s = 1.0 + 2.0*tanh(10*cos(pi*t)); % 0.5+2.5*tanh(10*cos(pi*(t)))
    end
end

function K_n = K_n(t)
    K_n = 0.8;
end

%%%%%%


function du = oim_dynamics(p, t, u)
    J = p{1};

    n = size(J,1);

    disp(u);

    % Add oscillator coupling component to dphi ---
    dot_product_matrix = -K(t) * J .* sin(u - u.');

    % Sum each row in this matrix to add oscillator coupling component to dphi
    du = sum(dot_product_matrix,2); 

    % Add SHIL SYNC component to dphi ---
    du = du - K_s(t).* sin(2 * u);

end