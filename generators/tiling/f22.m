function cube = f22()

% The F22 class of subproblem can be considered as just having
% 2 antiferromagnetic edges on OPPOSITE edges of the same face
% (this will in turn create 2 frustrated facets)

faces = get_faces(); % Get 6x12 matrix of 12 edges associated with each of 6 faces
cube = ones(1,12); % Set all edges to ferromagnetic by default

i = ceil(rand('single')*6); % Pick a random face
j = ceil(rand('single')*2); % Pick a random pair of edges on that face

index = [faces(i,j) faces(i,j+2)]; % Get the indices of the OPPOSITE (hence +2) edges to make antiferromagnetic

cube(index) = cube(index) - 2; % Set the edges to antiferromagnetic

end