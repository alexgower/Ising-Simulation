function cube = f21()

% The F21 class of subproblem can be considered as just having
% 1 antiferromagnetic edges somewhere in one of the 12 edges of the cube.

cube = ones(1,12); % i.e. set ferromagnetic edges by default
i = ceil(rand('single')*12); % pick a random edge to make antiferromagnetic
cube(i) = cube(i) - 2; % make it antiferromagnetic

end