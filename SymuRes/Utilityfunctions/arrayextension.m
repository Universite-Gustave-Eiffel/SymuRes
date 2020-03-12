function extarray = arrayextension(array,Nmin,type)
% extarray = listextension(array,Nmin,type)
% Extend a 2d matrix, cell or array along one of its dimensions (columns or
% rows) by repeating its content until the corresponding dimension has at
% least a minimum size
%
% INPUTS
%---- array : 1d or 2d matrix, cell or array
%---- Nmin  : integer, minimum size for the chosen dimension
%---- type  : string, 'column' or 'row', chosen dimension to extend
%
% OUTPUTS
%---- extarray : extended matrix, cell or array

if strcmp(type,'column')
    extarray = array;
    while size(extarray,2) < Nmin
        extarray = [extarray array];
    end
elseif strcmp(type,'row')
    extarray = array;
    while size(extarray,1) < Nmin
        extarray = [extarray; array];
    end
else
    warning('Bad extension type, must be one of the following: column, row')
    extarray = [];
end

end