function [E,C,G] = get_L_oim(v,X,W,fp,check,alpha,beta,fRBM) % THIS IS THE OIM-ADAPTED VERSION OF THE ORIGINAL get_L FUNCTION

    % FOR OIM DEFINITION OF v
    phi = v;
    % FOR PEI DEFINITION OF v
    % phi = acos(v); % THIS IS THE ANALOGOUS SOFT PHASES FOR THE SOFT SPIN VALUES OF v



    % Non-RBM case
    if ~fRBM

        % d = dimensions of lattice
        d = length(size(v));

        % V = 6x(nxmxk) (for 3D) matrix of the 6 nearest neighbors of each site in v
        % Remember v elements are continuous soft-spin PHASES between 0 and pi
        % Therefore v0 and V0 are the hard-spin +1 -1 versions of v and V
        % BUT NEED SIGMA=COS(PHI) CONVERSION FIRST BEFORE SIGN NOW FOR OIM
        % V = get_V(v); sigma_0 = sign(v); Sigma_0 = sign(V);
        Phi = get_V(phi); sigma_0 = sign(cos(phi)); Sigma_0 = sign(cos(Phi)); % WE HAVE CONFIRMED PHI CREATES SAME HARD SPINS
        

        % E calculates the ISING energy of the rounded hard-spin version of v
        % V0 .* W turns all s_j nearest neighbour values to W_{ij}s_{j}
        % sum(V0 .* W,d+1) sums over all nearest neighbours of each site to get \Sigma_{j \in N} W_{ij}s_{j}
        % v0 .* sum(V0 .* W,d+1) turns to \Sigma_{j \in N} s_{i}*W_{ij}*s_{j}
        % sum(v0.*sum(V0.*W,d+1),'all')/2 turns to \frac{1}{2} \Sigma_{i} \Sigma_{j \in N} s_{i}*W_{ij}*s_{j}
        % i.e. because W is symmetric, we need to renormalise by factor of 2
        % REMEMBER PEI TREATS GROUND STATE AS MOST POSITIVE ENERGY SO W=+1 IS FERROMAGNETIC BOND STILL
        % ALL THIS HARD SPIN STUFF STAYS THE SAME FOR OIM VERSION
        E = double(sum(sigma_0.*sum(Sigma_0.*W,d+1),'all')/2);


        % % lfield =  W_{ij}v_{j} but with the continuised soft-spin values of v
        % % is the LOCAL FIELD contributions from each nearest neighbour direction for each site v
        % % i.e. V = 6x(nxmxk) (for 3D) matrix of the 6 nearest neighbors of each site in v
        % % so lfield = W.*V; is the local field contribution from each nearest neighbour direction for each site v
        % WE WILL USE THIS FOR THE MEMORY VARIABLES CODE BUT NOT THE PHASE CODE
        % BUT WE NEED THE SIGMA = COS(PHI) CONVERSION FIRST FOR OIM
        % lfield = W.*cos(V); 

        % PEI DEFINITION OF V
        % lfield = W.*V;
        % lfield_phi = W.*cos(Phi);

        % OIM DEFINITION OF V
        % lfield = W.*cos(phi);


        % Non plaquette-based memory case (i.e. bond memory case)
        if ~fp
            % C_{ij} = (W_{ij}*\sigma_{i}*\sigma_{j} + 1) as in Pei paper (except factor of 1/2 absorbed into gamma)
            % i.e. lfield = 6x(nxmxk) (for 3D) matrix of the 6 nearest neighbours bond contributions of each site in v
            % i.e. v = nxmxk matrix of the soft-spin values of each site in v
            % Therefore v.*lfield = 6xnxmxk matrix of the local energy contributions from each nearest neighbour direction for each site v
            % (Remember C_ij lives on the bonds)
            % WE JUST NEED THE SIGMA=COS(PHI) CONVERSION FIRST FOR OIM
            % C = cos(v).*lfield+1;

            % C = v .* lfield + 1;
            % C = cos(phi) .* W .* cos(Phi) + 1;

            C = 0.5*W.*(cos(phi-Phi) + cos(phi+Phi)) + 1;
            
        
        % Plaquette-based memory case
        else
            % The checks just ensure that C and X are only non-zero in the black part of the checkerboard
            % C = wp_to_w(get_ww(v.*lfield),check,1)+1;
            % X = get_W(wp_to_w(X,check,0));
        end

        % ORIGINAL VERSION
        % G = \dot{\sigma_i} in the Pei paper or \dot{v} in the code nomenclature
        % i.e. the derivative of the continuised soft-spin v values with respect to time
        % THIS IS NOW G = \DOT{\phi} FOR OIM

        % i.e. the first term is gradient descent on the Hamiltonian due to the local field
        % i.e. the second term is gradient descent on the Hamiltonian due to the (memory controlled) softness constraint
        % Note that we sum over dimension d+1 i.e. over all nearest neighbours for each element in v (or equivalently indexed in G)
        
        % i.e. v - V = nxmxkx6 matrix of the differences in phases \phi between each site \phi_i and all of its 6 nearest neighbours \phi_j

        % G_sigma = alpha*sum(lfield,d+1) - 2*beta*sum(X,d+1).*v;
        % G_phi_from_sigma = -sin(phi).*(alpha*sum(lfield,d+1) - 2*beta*sum(X,d+1).*cos(phi));
        % G_sigma_from_phi = (-0.5*alpha*sum(W.*sin(phi-Phi),d+1) -0.5*alpha*sum(W.*sin(phi+Phi),d+1) + beta*sum(X,d+1).*sin(2*phi))./-sin(phi);

        % dH/d(phi_i) definition of dynamcis
        G = -0.5*alpha*sum(W.*sin(phi-Phi),d+1) -0.5*alpha*sum(W.*sin(phi+Phi),d+1) + beta*sum(X,d+1).*sin(2*phi);
        % dH/d(sigma_i) definition of dynamics (but with bad numerical precision with -1/sin(phi) jacobian factor)
        % G = (-0.5*alpha*sum(W.*sin(phi-Phi),d+1) -0.5*alpha*sum(W.*sin(phi+Phi),d+1) + beta*sum(X,d+1).*sin(2*phi))./(-sin(phi));




    % RBM case - TODO understand later
    else
        % [~,~,lfield,E,C] = get_E(v,W,true); 
        % G = alpha*lfield - 2*beta*[sum(X,2).' sum(X,1)].*v;
    end

end


%%%%% MATRIX NON-EQUALITY CHECKER CODE %%%%%
% A = lfield;
% B = lfield_phi;
% diffIndex = abs(A - B) > 1e-6;
% if any(diffIndex(:))
%     [row, col] = find(diffIndex);
%     for i = 1:length(row)
%         fprintf('A(%d,%d) = %g differs from B(%d,%d) = %g\n', ...
%                 row(i), col(i), A(row(i), col(i)), ...
%                 row(i), col(i), B(row(i), col(i)));
%     end
% else
%     disp('All elements are equal within the specified tolerance.');
% end