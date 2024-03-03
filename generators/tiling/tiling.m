function [W,Esol] = tiling(sz,flist)

d = length(sz); 

if flist == 'g' % Gaussian
    w = normrnd(0,1,[sz d]); % Sets right,back,up couplings per vertex v in 3D cubic lattice as Gaussian random in w
    Esol = 0; 
elseif flist == 'r' % Random {-1,1}
    w = -1 + 2*round(rand([sz d],'single')); % Sets right,back,up couplings per vertex v in 3D cubic lattice as random {-1,1} in w
    w = double(w); % Convert to double
    Esol = 0;
elseif d == 2 % 2D
    [w,Esol] = tiling_2d(sz(1),sz(2),flist);
elseif d == 3 % 3D F21,F22,F41,F42,F6 subproblem planted solution
    [w,Esol] = tiling_3d(sz(1),sz(2),sz(3),flist); 
end

% w = gauge_lattice(w); % TODO DELETE?

%%% CONVERT w MATRIX of size (n,m,k,3) (for 3D) i.e. right,back,up couplings per vertex v in 3D cubic lattice 
%%% to W MATRIX of size (n,m,k,12) i.e. left,right,front,back,down,up couplings per vertex v in 3D cubic lattice)
W = get_W(w);

end