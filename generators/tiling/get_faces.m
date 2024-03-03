function faces = get_faces() %% Returns the 4 edges associated with each face of the cube

faces = zeros(6,4); % 6 faces, 4 edges per face

% Each row represents a face, and the 4 columns represent the 4 edges 
% associated with that face.
faces(1,:) = [1 2 3 4]; 
faces(2,:) = [5 6 7 8];
faces(3,:) = [1 10 5 9];
faces(4,:) = [3 11 7 12];
faces(5,:) = [2 11 6 10];
faces(6,:) = [4 12 8 9];

end