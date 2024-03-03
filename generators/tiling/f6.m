function cube = f6()

% The F6 subproblem class has all 6 facets frustrated and best displayed in Fig2 of Hamze2018
% Therefore in general we can choose a random edge on face 1 (edges 1,2,3,4) to be antiferromagnetic
% And then just make sure that the other two edges on the other two antiferromagnetic edges are not diagonally opposite


cube = ones(1,12); % Set all edges as ferromagnetic by default

i = ceil(rand('single')*4); % Choose an edge on face 1 to be antiferromagnetic
j = ceil(rand('single')*2); % TODO

% Note mod4(x) = mod(x-1,4)+1 i.e. it is a 1 based modulo 4 function
% i.e. mod4(1) = 1, mod4(2) = 2, mod4(3) = 3, mod4(4) = 4, mod4(5) = 1, mod4(6) = 2 etc

% This messy bit of code uses the Pei convention to find the other two edges on other two faces
% that are antiferromagnetic (but none diagonally opposite to each other)
if j == 1
    i1 = mod4(i+1) + 4; 
    i2 = mod4(i+3) + 8;
else
    i1 = mod4(i+3) + 4;
    i2 = mod4(i+2) + 8;
end

cube([i i1 i2]) = cube([i i1 i2]) - 2; % Set the chosen edges to be antiferromagnetic

end