function W = get_W(w)
%%% CONVERT w MATRIX of size (n,m,k,3) (for 3D) i.e. right,back,up couplings per vertex v in 3D cubic lattice 
%%% to W MATRIX of size (n,m,k,6) i.e. left,right,front,back,down,up couplings per vertex v in 3D cubic lattice

% Get full size of w e.g. (n,m,k,3) for 3D
sz = size(w); 

% get class of w e.g. double, single etc
fd = class(w);

% dimension of lattice 'dim' = length(sz) - 1 as we ignore final parameter for 3 couplings (for 3D) etc
dim = length(sz)-1;

% 2D CASE
if dim == 2

    % (Extract n,m from sz)
    n = sz(1); m = sz(2);

    % Create W matrix of size (n,m,4) (initialise to zeros) describing couplings of each vertex v in 2D lattice to its 4 neighbours 
    W = zeros([n m 4],fd);

    % Firstly the couplings 2 and 4 are copied from w to W
    W(:,:,[2 4]) = w;
    % TODO UNDERSTAND THIS BIT
    W(:,:,1) = circshift(w(:,:,1),1,1);
    W(:,:,3) = circshift(w(:,:,2),1,2);


% 3D CASE
elseif dim == 3
   
    % Extract n,m,k from sz
    n = sz(1); m = sz(2); k = sz(3);

    % Create W matrix of size (n,m,k,6) (initialise to zeros) describing couplings of each vertex v in 3D lattice to its 6 neighbours
    W  = zeros([n m k 6],fd);

    % Firstly the couplings 2,4,6 (right,back,up) are copied from couplings 1,2,3 of w to W
    W(:,:,:,[2 4 6]) = w;


    % n1 = circshift(1:n,1); m1 = circshift(1:m,1); k1 = circshift(1:k,1); - DELETE?
    % W(:,:,:,[1 3 5]) = cat(4,w(n1,:,:,1),w(:,m1,:,2),w(:,:,k1,3)); - DELETE?

    % Note circshfit(M,1,1) shifts by 1 UNIT in dimension 1
    % All copulings 1 are taken from couplings 1 from w (for each v), but just shifted by 1 unit back in the dimension 1 (n)
    W(:,:,:,1) = circshift(w(:,:,:,1),1,1);
    % All copulings 3 are taken from couplings 2 from w (for each v), but just shifted by 1 unit back in the dimension 2 (m)
    W(:,:,:,3) = circshift(w(:,:,:,2),1,2);
    % All copulings 5 are taken from couplings 3 from w (for each v), but just shifted by 1 unit back in the dimension 3 (k)
    W(:,:,:,5) = circshift(w(:,:,:,3),1,3);
    
else
    
    fprintf('Error');
    return;
    
end

end