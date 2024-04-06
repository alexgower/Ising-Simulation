function check = checkerboard(sz,shade)

% d = dimension of lattice
d = length(sz);

% 2D CASE
if d == 2

    % We have an (n,m) 2D lattice
    n = sz(1); m = sz(2);

    check = reshape([ones(1,n/2); zeros(1,n/2)],[n 1]); % This just creates a column 1 0 1 0... of length n
    check = repmat([check ~check],[1 m/2]); % This repeats the two-column batch [1 0 1 0.. then 0 1 0 1...] m/2 times to fully checker the 2d lattice 
    
    % Shade is a boolean which determines whether the checkerboard is 'black or white' (1 or 0) first
    if shade
        check = ~check;
    end



% 3D CASE
elseif d == 3
    
    % We have an (n,m,k) 3D lattice
    n = sz(1); m = sz(2); k = sz(3);

    % This just creates a column 1 0 1 0... of length n
    check = reshape([true([1 n/2]); false([1 n/2])],[n 1]); 

    % This repeats the two-column batch [1 0 1 0.. then 0 0 0 0...] m/2 times to make one slice (U LATTICE checkerboard) tiling of the 3D lattice
    check1 = repmat([check false(n,1)],[1 m/2]);  

    % This repeats the two-column batch [0 0 0 0.. then 0 1 0 1...] m/2 times to make one opposite slice (U LATTICE checkerboard) tiling of the 3D lattice
    check2 = repmat([false(n,1) ~check],[1 m/2]);

    if shade % If shade is true, we want to swap the two slices order
        check = check1; check1 = check2; check2 = check;
    end

    % cat(3,check1,check2) creates an alternating TWO-SLICE (U LATTICE checkerboard) pattern
    % repmat(...,[1 1 k/2]) repeats the two-slice pattern k/2 times to fully checker the 3D (U LATTICE) lattice
    check = repmat(cat(3,check1,check2),[1 1 k/2]);

end

end