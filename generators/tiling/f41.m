function cube = f41()

% The F41 subproblem class is best dispalyed in Fig 2 of Hameze2018

cube = ones(1,12); % Set all edges to ferromagnetic by default
faces = get_faces(); % Get a 6x12 matrix of the 12 edges for each of the 6 faces

fi = 2*ceil(rand('single')*3)-1; % Pick a random first face out of 1, 3, 5 (the second face will be fi+1)
i = ceil(rand('single')*2); 

j1 = ceil(rand('single')*4); % Pick a random edge on the first face to make antiferromagnetic
% Pick the edge on the second face to make antiferromagnetic (can't be same or opposite edge to j1)
% (note this either adds 1 or 3 (mod 4) to the relative (1,2,3,4) edge index on the face)
j2 = mod4(j1+2*i-1); 

i1 = faces(fi,j1); % Get the index of the edge on the first face to make antiferromagnetic
i2 = faces(fi+1,j2); % Get the index of the edge on the second face to make antiferromagnetic

cube([i1 i2]) = cube([i1 i2]) - 2; % Set the edges to antiferromagnetic

end