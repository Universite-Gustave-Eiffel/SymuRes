function revlist = reverselist(list)
% revlist = reverselist(list)
% Return the element of a list in reverse order
%
% INPUTS
%---- list : vector, list of elements
%
% OUTPUTS
%---- revlist : vector, same size as list, elements of list in reverse order

[nrow, ncol] = size(list); % list may be column or row vector
revlist = zeros(nrow,ncol);
N = max([nrow ncol]);

for i = 1:N
    revlist(i) = list(N-i+1);
end

end