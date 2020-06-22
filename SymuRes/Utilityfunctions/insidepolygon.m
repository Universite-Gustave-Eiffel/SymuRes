function isinside = insidepolygon(pt,poly)
% isinside = insidepolygon(pt,poly)
% Test if a point (x,y) is inside a polygon or not
%
% INPUT
%---- pt   : [x y] vector, point to test
%---- poly : 2-by-N or N-by-2 matrix, list of the points defining the
%            polygon (vertices)
%
% OUTPUT
% isinside : boolean, 1 if inside the polygon and 0 otherwise

isinside = 0;

[n, m] = size(poly);
[Npts, idim] = max([n m]); 

if Npts > 2 % the polygon must be at least a triangle
    
    % Get the x-values and y-values of the polygon vertices
    if idim == 1
        xpoly = poly(:,1)';
        ypoly = poly(:,2)';
    else
        xpoly = poly(1,:);
        ypoly = poly(2,:);
    end
    
    % Test if the point is inside the rectangle that includes the polygon
    if min(xpoly) <= pt(1) && pt(1) <= max(xpoly) && min(ypoly) <= pt(2) && pt(2) <= max(ypoly)
        j = Npts;
        for i = 1:Npts
            if (ypoly(i)>pt(2)) ~= (ypoly(j)>pt(2)) && pt(1) < xpoly(i)+(xpoly(j)-xpoly(i))/(ypoly(j)-ypoly(i))*(pt(2)-ypoly(i))
                isinside = ~isinside;
            end
            j = i;
        end
    end
end

end