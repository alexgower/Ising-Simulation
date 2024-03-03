function [w,E,cost] = tiling_3d(n,m,k,flist)

if flist == 10 % Special case of ferromagnetic Ising model (CHECK)
    w = ones(n,m,k,3);
    E = sum(w(:));
    cost = (sum(abs(w(:)))-E)/2;
    return;
end
if flist == 11 % Special case of fully frustrated Ising model MAYBE (CHECK)
    cube = f6(); cube = permute(cube,[1 4 3 2]);
    wcube = repmat(cube,[n m k 1]);
    w = wp_to_w(wcube);
    E = sum(w(:));
    cost = (sum(abs(w(:)))-E)/2;
    return;
end

% Convert flist of subproblem ratios to cdf (cumulative distribution function)
cdf = to_cdf(flist); 

% Initialize the 3D array of unitcubes
% (i.e. for each vertex v in the 3D lattice, assign 12 values corresponding to the 12 edge values 
% of the 3D unit cube C with v at bottom-left-front (etc CHECK))
% !!! IN PRINCIPLE THIS INCLUDES SPACES FOR 'GAP UNIT CUBES' (i.e. cubes that are not actually used as subproblems)
wcube = zeros(n,m,k,12); 

% n = round(n/2); m = round(m/2); k = round(k/2); % (DELETE?)

% Loop over all vertices v in the 3D lattice
for i = 1:n 
for j = 1:m 
for l = 1:k
   
    % Now only consider vertices u that have same parity co-ordinates
    if (mod(i,2) && mod(j,2) && mod(l,2)) || (~mod(i,2) && ~mod(j,2) && ~mod(l,2))

        % Sample a subproblem for the unit cube C defined by this vertex v/u
        x = cdf_sample(cdf);
        if x == 1  
            cube = f21();
        elseif x == 2 
            cube = f22();
        elseif x == 3 
            cube = f41();
        elseif x == 4 
            cube = f42();
        elseif x == 5 
            cube = f6();
        end

        % Assign the subproblem weights to the 3D array of unitcubes
        wcube(i,j,l,:) = cube;
    end
    
end
end
end

%%% THIS REARRANGES THE EDGES DEFINITION FROM A WAY CONVENIENT TO DEFINE SUBPROBLEMS TO A WAY CONVENIENT TO CALCULATE THE 3 COUPLINGS PER VERTEX
wcube(:,:,:,:) = wcube(:,:,:,[1 4 9 7 2 10 3 6 12 5 8 11]);

%%% THEN CONVERT wpcube to w
% ULTIAMTELY THIS FUNCTION TAKES A wpcube MATRIX (DEFINED AS 12 EDGE VALUES PER UNIT CUBE C DEFINED AT EACH V=(n,m,k) LATTICE SITE)
% (BUT WITH 1-12 CONVENTION LABELS RE-ORDERED TO BE CONVENIENT FOR VERTEX COUPLINGS NOT SUBPROBLEM DEFINITIONS)
% AND CONVERTS IT TO A w MATRIX (DEFINED AS 3 EDGE VALUES PER LATTICE SITE V=(n,m,k) IN THE RIGHT, BACK, UP DIRECTIONS (ETC)) 
w = wp_to_w(wcube,0,0);

% The planted ground state energy is simply the sum of all the weights as all are satisfied by the ground state
% Note Pei measures ground state energy as MOST POSITIVE
E = sum(w(:));

% The cost measures the number of unsatisfied constraints in the planted ground state
cost = (sum(abs(w(:)))-E)/2;

end