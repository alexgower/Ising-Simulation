function [idx,ind] = discard(idx)

% This function takes idx and discards single element clusters,
% Modifying idx such that the indexes are all 0 for these single element clusters
% and also returning ind which is the list of all the UNIQUE cluster indices that are not single element clusters

% Find maximum cluster index
m = max(idx); 

% Create empty indices list
ind = [];

% Iterate through each cluster INDEX value, not the values in idx itself
for i = 1:m
    % Find all elements in idx that are in this cluster
    ii = find(idx == i);

    % If it is a single element cluster, rename it's index value as 0 in idx and don't add it's index value to the indices ind
    if length(ii) == 1
        idx(ii) = 0;
    
    % If it is not a single element cluster, add it's index value to the indices ind (and don't change it's index value in idx)
    elseif ~isempty(ii)
        ind = [ind i]; 
    end
end

end