function w = wp_to_w(wp,check,map)

% sz is full lattice size i.e. (n,m,k,12) for 3D and (n,m,4) for 2D
sz = size(wp); 
% Ignore the 12/4 at the end so now sz = (n,m,k) for 3D and (n,m) for 2D
sz = sz(1:end-1); 
% d is the dimension of the lattice
d = length(sz);

% If check is not already a checkboard matrix, take 'check' as the 'shade' value and create a checkerboard matrix with it
% TODO UNDERSTAND LATER
if numel(check) == 1
    check = checkerboard(sz,check);
end


% If the lattice is 2D, call the 2D function, else call the 3D function
% Map seems to be a bool (0 or 1) value which determines the way to reduce wp to w - CHECK
if d == 2  
    w = wp_2d(wp,check,map);  
elseif d == 3
    w = wp_3d(wp,check,map);
end

end
