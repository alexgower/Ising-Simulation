function cdf = to_cdf(list)

% This function takes a list of ratio values and returns the cumulative
% distribution function (CDF) of the list.
% It does this by integrating over the list of ratio values (which are essentially unnormalised PDF values,
% but unnormalised by the same factor) and then normalising the whole CDF afterwards

l = length(list); % length of the list (number of ratio values)
listsum = zeros(1,l); % initialise unnormalised CDF as list of zeros
listsum(1) = list(1); % first value of unnormalised CDF is the first value of the list

for i = 2:l
   listsum(i) = sum(list(1:i)); % i'th value of unnormalised CDF is the sum of the first i values of the original list of ratio values
end

cdf = listsum/listsum(l); % normalise the unnormalised CDF by dividing by the last value of the unnormalised CDF

end