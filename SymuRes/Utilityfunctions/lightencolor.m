function newcolor = lightencolor(color0,factor)
% newcolor = lightencolor(color0,factor)
% Lighten or brighten one or more RGB colors
%
% INPUTS
%---- color0 : n-by-3 matrix, the columns are the RGB components in [0,1]
%---- factor : scalar between 0 and 1. 0: no change, 1: white
%
% OUTPUTS
%---- newcolor : n-by-3 matrix, the modified colors

Ncolor = size(color0,1);
Ncomponent = size(color0,2);

newcolor = zeros(Ncolor,Ncomponent);

for i = 1:Ncolor
    for j = 1:Ncomponent
        newcolor(i,j) = color0(i,j) + factor*(1 - color0(i,j));
    end
end

end