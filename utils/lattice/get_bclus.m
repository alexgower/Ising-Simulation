function [list,b,bg,cn] = get_bclus(w,sw)

% Note that bclus uses a (n,m,k,3) w matrix not a (n,m,k,6) W matrix (i.e only right,back,up not left,right,front,back,down,up)
% so that we don't double count the same bond

% Also note bclus stands for 'bond clusters' i.e. the w matrix in question is an x (from X) matrix of bond memory values
% (i.e. larger x values indicate a measured stronger bond between spins - but then simultaneously larger x will soften the spins)

% Note that in 'record' function, we actually don't care about 'list' of spins in some cluster
% Note than in general we don't care about 'cn' list of spins in first cluster


sz = size(w); 
if length(sz) == 2 % 1D case
    N = sum(sz); % MAYBE - understand this later (why they include 4th dimension in sum)
else
    N = prod(sz(1:end-1)); % Number of spins N = n*m*k (etc for 3D)
end


cn = zeros([1 N]); % cn later on will be a logical list describing the spins in cluster 1 - CHECK why

% Turn bond memory values x into a logical matrix of active (1) and inactive (0) bonds
% Also add random noise to the bond memory values first (probably to get new clusters even if the bond memory values are the same as before? - CHECK)
w = logical(floor( w + rand(sz,'single') )); 

% TODO - understand this later
G = get_G(w);

% Get the connected components of the graph G (i.e. the clusters of spins)
% idx is a list of size (1,N) with a cluster index for each spin
% bsz is a list of size (1,N) with the size of each bond cluster
[idx,bsz] = conncomp(G);

% If there are no clusters then return empty list (and empty list of cluster sizes, and empty list of number of clusters of each size)
if isempty(bsz)
    list = []; b = 0; bg = [];
    return;
end

% bg = a unique list of the bond cluster sizes
% b = the corresponding number of bond clusters of each size
[b,bg] = my_gp(bsz.');

% Set cn to 1 for all spins in the first cluster in the list of cluster indices
cn(idx == idx(1)) = 1;

% Non-Swedensen-Wang Case
if ~sw 

    % Randomly select a cluster index from the list of cluster indices
    % (with probability proportional to the size of the cluster - i.e. the number of spins in the cluster)
    % i.e. pick a random spin and then pick the cluster index of that spin
    ind = cdf_sample(to_cdf(bsz));
    % Return a list of 1s and 0s indicating whether the cluster index of each spin is the randomly selected cluster index
    list = find(idx == ind);

% Swedensen-Wang Case
else
    [idx,ind] = discard(idx); % idx is a modified version of idx with single element cluster indexes all set to 0, ind is the list of non-zero indices of idx
    
    % Randomly keep only ~half (as rand 50/50 rounds to 0 or 1) of elements in ind
    % (I think this is analogous to flipping clusters with probability of 0.5 in SW algorithm)
    ind = ind(logical(round(rand([1,length(ind)],'single')))); 

    % If there are no clusters to flip (either randomly not selected, or all single element clusters of both) then return empty list
    if isempty(ind)
        list = 0;
        return;
    end

    % Create a matrix of size (length(ind),N) where we repeat length(ind) rows of idx 
    % (remember idx is a list of size (1,N) with a cluster index for each spin (or 0 if single element cluster))
    idx = repmat(idx,[length(ind) 1]); 
    
    % Create a matrix of size (length(ind),N) where we repeat length(ind) columns of ind
    % (remember ind is a list of size (1,length(ind)) with the UNIQUE cluster indices of non-zero clusters)
    ind = repmat(ind.',[1 N]);
    
    % e.g. the first column of idx now will just be the cluster index of spin 1 repeated length(ind) times
    % e.g. the first column of ind now will just be ind itself i.e. the list of cluster indices to flip
    % Therefore the value of sum(idx==ind,1) for this column will be 1 if the cluster index of spin 1 is in ind, 0 otherwise
    % And sum(idx==ind,1) will be a list of length N with 1s and 0s indicating whether the cluster index of each spin is in ind
    % i.e. whether it should be flipped or not
    list = sum(idx==ind,1);

end

end