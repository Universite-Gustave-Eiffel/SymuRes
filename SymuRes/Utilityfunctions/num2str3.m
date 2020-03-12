function s = num2str3(x)
% s = num2str3(x)
% Delete the point when converting a float to a string
% Ex: 12.36 converted into 1236; 0.102 converted into 0102
% Limitation to 5 digits for the decimal part
% Ex: 3.141596 converted into 314159.6
%
% INPUTS
%---- x : float number
%
% OUTPUTS
%---- s : string

intx = floor(x);
maxiter = 5;
if x == intx
    s = num2str(x);
else
    n = 1;
    decx = (x - intx)*10^n;
    while n <= maxiter && decx ~= floor(decx)
        n = n + 1;
        decx = (x - intx)*10^n;
    end
    s = [num2str(intx) num2str(decx)];
end

end
    