function extarray = arrayextension(array,Next,type)
% extarray = listextension(array,Next,type)
% Extend a 2d matrix, cell or array along one of its dimensions (columns or
% rows) by repeating its content until the corresponding dimension has the
% required size Next
%
% INPUTS
%---- array : 1d or 2d matrix, cell or array
%---- Next  : integer, desired size for the chosen dimension
%---- type  : string, 'column' or 'row', chosen dimension to extend
%
% OUTPUTS
%---- extarray : extended matrix, cell or array

if strcmp(type,'column')
    extarray = array;
    while size(extarray,2) < Next
        extarray = [extarray array];
    end
    extarray = extarray(:,1:Next);
elseif strcmp(type,'row')
    extarray = array;
    while size(extarray,1) < Next
        extarray = [extarray; array];
    end
    extarray = extarray(1:Next,:);
else
    warning('Bad extension type, must be one of the following: column, row')
    extarray = [];
end

end