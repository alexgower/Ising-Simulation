function w = get_ww(W)

% This function turns a W matrix [of size (n,m,k,6) for 3D etc] which includes all nearest neighbour couplings per site in v
% into a w matrix [of size (n,m,k,3) for 3D etc] which includes just couplings in the right,back,up directions etc


% Remember W is a matrix of size (n,m,k,6) for 3D etc
sz = size(W);
% d = dimension of lattice
d = length(sz)-1;


if d == 1 % 1D case
    w = W; % Trivial case
elseif d == 2 % 2D Case
    w = W(:,:,[2 4]); % Just include right and up couplings
elseif d == 3 % 3D Case
    w = W(:,:,:,[2 4 6]); % Just include right, back and up couplings
else   
    fprintf('Error');
    return;    
end

end