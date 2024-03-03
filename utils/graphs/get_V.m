function V = get_V(v)

% GET_V - Returns the 3D array V of the 4/6 nearest neighbors of the 2D/3D array of spins v 

% sz = (n,m,k) etc for 3D, d = dimension of lattice
sz = size(v); d = length(sz);

    % 2D Case - TODO think later
    if d == 2

        % V = cat(3,circshift(v,-1,1),circshift(v,-1,2)); % DELETE?
            
        V  = cat(3, circshift(v,1,1), circshift(v,-1,1), ...
                    circshift(v,1,2), circshift(v,-1,2));
        

    % 3D Case
    elseif d == 3
            
        % V = cat(4,circshift(v,-1,1),circshift(v,-1,2),circshift(v,-1,3));  % DELETE?
            

        % cat(4,... stores in some new overarching '4th dimension' (with 6 rows corresponding to 6 rows corresponding to 6 nearest neighbor types)
        % 6 different sub-matrices of the size of the 3D array v, but with values corresponding to the 6 nearest neighbors of the 3D array v
        % Which we get from a +-1 circshift (periodic BCs) in dimensions 1/2/3
        V  = cat(4, circshift(v,1,1), circshift(v,-1,1), ...
                    circshift(v,1,2), circshift(v,-1,2), ...
                    circshift(v,1,3), circshift(v,-1,3));


    else
        
        fprintf('Error');
        return;
        
    end

end