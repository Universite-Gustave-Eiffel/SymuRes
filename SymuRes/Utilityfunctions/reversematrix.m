function revM = reversematrix(M,dim)
% revlist = reverselist(list)
% Return the element of a list in reverse order
%
% INPUTS
%---- M   : 2d matrix, matrix of elements
%---- dim : integer, 1: reverse the rows, 2: reverse the columns 
%
% OUTPUTS
%---- revM : matrix, same size as M, elements in reverse order along
%            dimension 1 or 2

[nrow, ncol] = size(M); % list may be column or row vector
revM = zeros(nrow,ncol);

if dim == 1
    for i = 1:nrow
        revM(i,:) = M(nrow-i+1,:);
    end
elseif dim == 2
    for j = 1:ncol
        revM(:,j) = M(:,ncol-j+1);
    end
end

end