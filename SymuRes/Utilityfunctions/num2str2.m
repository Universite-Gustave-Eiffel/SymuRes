function s = num2str2(x)
% s = num2str2(x)
% Delete the point when converting a float to a string
% Ex: 12.36 converted into 1236; 0.102 converted into 0102
%
% INPUTS
%---- x : float number
%
% OUTPUTS
%---- s : string

sx = num2str(x);
Ns = length(sx);
s = '';
for i = 1:Ns
    if ~strcmp(sx(i),'.')
        s = [s sx(i)];
    end
end

end
    