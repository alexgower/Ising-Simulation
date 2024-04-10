% Function to convert (n,m,k,6) matrix W which contains nearest neighbour couplings for all nxmxk spins in a cubic lattice
% to a (n*m*k, n*m*k) matrix J which contains the couplings between all spins in the lattice
% W is a symmetric matrix and so J will be symmetric as well and so we still need the factor of 1/2 in the energy function as usual
% Similarly W and J both use the +J convention for favourable bonds so we still need the factor of -1 in the energy function as usual

function J = convert_W_to_J(W)
    % Assuming n, m, k are the dimensions of the 3D lattice and W is your 4D matrix
    [n, m, k, ~] = size(W);

    % Initialize J
    N = n * m * k;
    J = zeros(N, N);

    % Function to map 3D indices to a single 1D index, now correctly including 'n', 'm', and 'k'
    map_3d_to_1d = @(n_0, m_0, k_0, n, m, k) (k_0-1)*n*m + (m_0-1)*n + n_0;

    % Populate J
    for n_0 = 1:n
        for m_0 = 1:m
            for k_0 = 1:k
                i = map_3d_to_1d(n_0, m_0, k_0, n, m, k);
                % Neighbors: left, right, front, back, down, up
                neighbors = [
                    map_3d_to_1d(mod(n_0-2, n)+1, m_0, k_0, n, m, k), % Left
                    map_3d_to_1d(mod(n_0, n)+1, m_0, k_0, n, m, k),   % Right
                    map_3d_to_1d(n_0, mod(m_0-2, m)+1, k_0, n, m, k), % Front
                    map_3d_to_1d(n_0, mod(m_0, m)+1, k_0, n, m, k),   % Back
                    map_3d_to_1d(n_0, m_0, mod(k_0-2, k)+1, n, m, k), % Down
                    map_3d_to_1d(n_0, m_0, mod(k_0, k)+1, n, m, k)    % Up
                ];
                
                for index = 1:6
                    j = neighbors(index);
                    J(i, j) = W(n_0, m_0, k_0, index);
                end
            end
        end
    end

    if ~issymmetric(J)
        error('J is not symmetric. Please check the population logic.');
    end
end