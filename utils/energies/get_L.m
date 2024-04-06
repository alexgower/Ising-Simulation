function [E,C,G] = get_L(v,X,W,fp,check,alpha,beta,fRBM)

    % Non-RBM case
    if ~fRBM

        % d = dimensions of lattice
        d = length(size(v));

        % V = 6x(nxmxk) (for 3D) matrix of the 6 nearest neighbors of each site in v
        % Remember v elements are continuous soft-spins between -1 and 1
        % Therefore v0 and V0 are the hard-spin versions of v and V
        V = get_V(v); v0 = sign(v); V0 = sign(V);

        % E calculates the ISING energy of the rounded hard-spin version of v
        % V0 .* W turns all s_j nearest neighbour values to W_{ij}s_{j}
        % sum(V0 .* W,d+1) sums over all nearest neighbours of each site to get \Sigma_{j \in N} W_{ij}s_{j}
        % v0 .* sum(V0 .* W,d+1) turns to \Sigma_{j \in N} s_{i}*W_{ij}*s_{j}
        % sum(v0.*sum(V0.*W,d+1),'all')/2 turns to \frac{1}{2} \Sigma_{i} \Sigma_{j \in N} s_{i}*W_{ij}*s_{j}
        % i.e. because W is symmetric, we need to renormalise by factor of 2
        E = double(sum(v0.*sum(V0.*W,d+1),'all')/2);


        % lfield =  W_{ij}v_{j} but with the continuised soft-spin values of v
        % is the LOCAL FIELD contributions from each nearest neighbour direction for each site v
        % i.e. V = 6x(nxmxk) (for 3D) matrix of the 6 nearest neighbors of each site in v
        % so lfield = W.*V; is the local field contribution from each nearest neighbour direction for each site v
        lfield = W.*V; 

        % Non plaquette-based memory case (i.e. bond memory case)
        if ~fp
            % C_{ij} = (W_{ij}*v_{i}*v_{j} + 1) as in Pei paper (except factor of 1/2 absorbed into gamma)
            % i.e. lfield = 6x(nxmxk) (for 3D) matrix of the 6 nearest neighbours bond contributions of each site in v
            % i.e. v = nxmxk matrix of the soft-spin values of each site in v
            % Therefore v.*lfield = 6xnxmxk matrix of the local energy contributions from each nearest neighbour direction for each site v
            % (Remember C_ij lives on the bonds)
            C = v.*lfield+1;
        
        % Plaquette-based memory case
        else
            % The checks just ensure that C and X are only non-zero in the black part of the checkerboard
            C = wp_to_w(get_ww(v.*lfield),check,1)+1;
            X = get_W(wp_to_w(X,check,0));
        end


        % G = \dot{\sigma_i} in the Pei paper or \dot{v} in the code nomenclature
        % i.e. the derivative of the continuised soft-spin v values with respect to time
        % i.e. the first term is gradient descent on the Hamiltonian due to the local field
        % i.e. the second term is gradient descent on the Hamiltonian due to the (memory controlled) softness constraint
        % Note that we sum over dimension d+1 i.e. over all nearest neighbours for each element in v (or equivalently indexed in G)
        G = alpha*sum(lfield,d+1) - 2*beta*sum(X,d+1).*v;






    % RBM case - TODO understand later
    else
        [~,~,lfield,E,C] = get_E(v,W,true); 
        G = alpha*lfield - 2*beta*[sum(X,2).' sum(X,1)].*v;
    end

end