function w = wp_3d(wp,check,map)

% ULTIAMTELY THIS FUNCTION TAKES A 'problem based' wpcube MATRIX (DEFINED AS 12 EDGE VALUES PER UNIT CUBE C DEFINED AT EACH V=(n,m,k) LATTICE SITE)
% (BUT WITH 1-12 CONVENTION LABELS RE-ORDERED TO BE CONVENIENT FOR VERTEX COUPLINGS NOT SUBPROBLEM DEFINITIONS)
% AND CONVERTS IT TO A w MATRIX (DEFINED AS 3 EDGE VALUES PER LATTICE SITE V=(n,m,k) IN THE RIGHT, BACK, UP DIRECTIONS (ETC))    

% sz is full lattice size i.e. (n,m,k,12) for 3D 
% Ignore the 12/4 at the end so now sz = (n,m,k) for 3D
% So now sz = (n,m,k)
sz = size(wp); sz = sz(1:end-1); 

% fd is the precision level of wp (e.g. 'single' or 'double')
fd = class(wp);

% Explicitly extract the lattice size (n,m,k) from sz
n = sz(1); m = sz(2); k = sz(3);

% n1/m1/k1 gets all the lattice vectors shifted by 1 in the NEGATIVE direction (and circshift handles periodic BCs)
% n2/m2/k2 gets all the lattice vectors shifted by 1 in the POSITIVE direction (and circshift handles periodic BCs)
n1 = circshift(1:n,1); m1 = circshift(1:m,1); k1 = circshift(1:k,1);
n2 = circshift(1:n,-1); m2 = circshift(1:m,-1); k2 = circshift(1:k,-1);


if ~map % i.e. if map = 0

    % Multiply wp matrix element-wise by checkerboard matrix (checkerboard defined only on U like cube subproblems)
    % (Note that since cubes are EDGE DISJOINT then checkerboard matrix doesn't affect validity of planted solution)
    % I think this is just here to delete all subcubes that are not in U (for check=0) or offset-U (for check=1) - DOUBLE-CHECK
    % (e.g. for convenience if we want to define over all V cubes initially)
    wp = wp.*check;


    % NOW WE CONVERT wpcube - defined as 12 edge values per unit cube C defined at each v=(n,m,k) lattice site
    % TO w - defined as 3 edge values per lattice site v=(n,m,k) in the right, back, up directions (etc)

    % But note since each edge value is shared by two RAW (in V) unit cubes, (but obviously only in 1 in U unit cube) 
    % There are 4 RAW (in V) unit cubes where these 3 edge values could come from (although in reality only 1 contribution will be non-zero for each edge value)

    % My ''Pei Convention' picture makes the specific numbers used here obvious
    % For each v=(n,m,k) lattice site, firstly just keep the [1,2,3] edge values of the unit cube defined at v if v is in U
    w = wp(:,:,:,[1 2 3]) + ...
        cat(4,wp(:,m1,k1,4),wp(n1,:,:,[5 6])) + ... % Then add the [5,6] edge values from unit cube 1 behind in n direction, and [4] value 1 behind in m and k direction
        cat(4,wp(:,m1,:,7),wp(n1,:,k1,8),wp(:,m1,:,9)) + ... % Then add the [7] value 1 behind in m direction, [8] value 1 behind in n and k direction, and [9] value 1 behind in m direction
        cat(4,wp(:,:,k1,[10 11]),wp(n1,m1,:,12)); % Then add the [10,11] edge values from unit cube 1 behind in k direction, and [12] value 1 behind in n and m direction



else % i.e. if map = 1 - TODO CHECK what this case does

    w = zeros([sz 12],fd);
    w(:,:,:,[1 2 3]) = wp;
    w(:,:,:,[5 6]) = wp(n2,:,:,[2 3]);
    w(:,:,:,[7 9]) = wp(:,m2,:,[1 3]);
    w(:,:,:,[10 11]) = wp(:,:,k2,[1 2]);
    w(:,:,:,4) = wp(:,m2,k2,1);
    w(:,:,:,8) = wp(n2,:,k2,2);
    w(:,:,:,12) = wp(n2,m2,:,3);
    w = w.*check;

end

end