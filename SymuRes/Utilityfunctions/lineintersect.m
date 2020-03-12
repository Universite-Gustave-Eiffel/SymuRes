function [isinter, xinter, yinter] = lineintersect(pt11,pt12,pt21,pt22)
% [isinter, xinter, yinter] = lineintersect(pt11,pt12,pt21,pt22)
% Return the intersection point between two lines defined by the points
% (pt11,pt12) and (pt21,pt22).
% 
% INPUTS
%---- pt11, pt12 : points [x y] defining the first line
%---- pt21, pt22 : points [x y] defining the second line
%
% OUTPUTS
%---- isinter : boolean, 1: [pt11,pt12] crosses [pt21,pt22], 0: otherwise
%---- xinter  : scalar, abscissa of the intersection point
%---- yinter  : scalar, ordinate of the intersection point

x11 = pt11(1);
y11 = pt11(2);
x12 = pt12(1);
y12 = pt12(2);

x21 = pt21(1);
y21 = pt21(2);
x22 = pt22(1);
y22 = pt22(2);

isinter = 0;

% Test if [pt11 pt12] intersects [pt21 pt22]
det0 = (x11 - x12)*(y21 - y22) - (y11 - y12)*(x21 - x22);
if det0 ~= 0
    % intersection point of the lines (pt11,pt12) and (pt21,pt22)
    xinter = ((x11*y12 - y11*x12)*(x21 - x22) - (x11 - x12)*(x21*y22 - y21*x22))/det0;
    yinter = ((x11*y12 - y11*x12)*(y21 - y22) - (y11 - y12)*(x21*y22 - y21*x22))/det0;
    its1 = 0;
    if x11 < x12
        if x11 <= xinter && xinter <= x12
            its1 = 1;
        end
    elseif x12 < x11
        if x12 <= xinter && xinter <= x11
            its1 = 1;
        end
    else
        if y11 < y12
            if y11 <= yinter && yinter <= y12
                its1 = 1;
            end
        elseif y12 < y11
            if y12 <= yinter && yinter <= y11
                its1 = 1;
            end
        end
    end

    its2 = 0;
    if x21 < x22
        if x21 <= xinter && xinter <= x22
            its2 = 1;
        end
    elseif x22 < x21
        if x22 <= xinter && xinter <= x21
            its2 = 1;
        end
    else
        if y21 < y22
            if y21 <= yinter && yinter <= y22
                its2 = 1;
            end
        elseif y22 < y21
            if y22 <= yinter && yinter <= y21
                its2 = 1;
            end
        end
    end

    if its1 == 1 && its2 == 1
        isinter = 1;
    end
    
else
    % The two lines are parallel
    xinter = Inf; % default
    yinter = Inf; % default
end

end