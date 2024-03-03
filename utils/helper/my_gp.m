function [b,bg] = my_gp(list)

% TODO description here
% Note list will be bsz i.e. sizes of each bond cluster

% ll is length of bsz i.e. number of clusters
ll = length(list);

% Inititialize b and bg as empty
% bg is a list of cluster sizes
% b is the corresponding NUMBER of clusters of each size
b = []; bg = []; 

% Loop over all clusters, start at the the value at the end of the bsz list
while ll>0
   
    % Find the cluster size of this cluster
    li = list(ll);
    % Find all the indices in bsz of this size
    ind = find(list==li);
    % l is the number of clusters of this size
    l = length(ind);
    
    % Append l to b (i.e. the number of clusters of this size)
    b = [b l]; 
    % Append li to bg (i.e. the size of the clusters of this size)
    bg = [bg li];
    
    % Remove all indices of this size from bsz
    list(ind) = [];

    % Decrement to a previous index in the bsz list (until you hit 0 then exit the loop)
    ll = ll-l;
    
end

end