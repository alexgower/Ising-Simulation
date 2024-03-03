function cube = f42()

% The F42 subproblem class consists of all ferromagnetic edges except for 2 antiferromagnetic edges
% which are DIAGONALLY OPPOSITE each other. There are 6 possible subproblems in this class.

cube = ones(1,12); % Set all edges to ferromagnetic by default

% Choose one of the 6 possible subproblems
i = ceil(rand('single')*6); 

% This messy bit of code just uses the Pei convention to choose these diagonally opposite edges indices
if i <= 4
    i1 = i;
    i2 = mod4(i+2) + 4;

    % Note mod4(x) = mod(x-1,4)+1 i.e. it is a 1 based modulo 4 function
    % i.e. mod4(1) = 1, mod4(2) = 2, mod4(3) = 3, mod4(4) = 4, mod4(5) = 1, mod4(6) = 2 etc

else % i = 5 or 6
    i1 = i + 4;
    i2 = i1 + 2;
end

cube([i1 i2]) = cube([i1 i2]) - 2; % Set the selected edges to antiferromagnetic

end