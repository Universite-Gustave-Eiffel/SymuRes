function cmap = buildmap(color1,color2,nbColor)
%
% Build a colormap that corresponds to the transition between two colors
%
% INPUTS
%---- color1  : 3-size row vector, first color in RGB format
%---- color2  : 3-size row vector, second color in RGB format
%---- nbColor : integer, number of colors in the colormap
%
% OUTPUTS
%---- cmap : nbColor-by-3 matrix, colormap built

diffcolor = color2 - color1;
nbhalf = floor((nbColor+1)/2);
cmap = ones(nbColor,1)*(color1 + color2);
cmap = (cmap <= 1).*cmap + (cmap > 1);

addcanal = [];
delcanal = [];
for canal = 1:3
    if diffcolor(canal) > 0
        addcanal = [addcanal canal];
    elseif diffcolor(canal) < 0
        delcanal = [delcanal canal];
    end
end

if ~isempty(addcanal) && ~isempty(delcanal)
    for i = 1:length(addcanal)
        canal = addcanal(i);
        cmap(1:nbhalf,canal) = linspace(color1(canal),color2(canal),nbhalf);
    end
    for i = 1:length(delcanal)
        canal = delcanal(i);
        cmap(nbhalf:nbColor,canal) = linspace(color1(canal),color2(canal),nbColor-nbhalf+1);
    end
elseif isempty(delcanal)
    for i = 1:length(addcanal)
        canal = addcanal(i);
        cmap(1:nbColor,canal) = linspace(color1(canal),color2(canal),nbColor);
    end
elseif isempty(addcanal)
    for i = 1:length(delcanal)
        canal = delcanal(i);
        cmap(1:nbColor,canal) = linspace(color1(canal),color2(canal),nbColor);
    end
end

end
    