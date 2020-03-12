function [xpath, ypath] = smoothroute(xroute,yroute,Nptsbezier,alphastart,alphaend)
% [xpath, ypath] = smoothroute(xroute,yroute,Nptsbezier,alphaturn,alphaedge)
% Return a smooth path representing a given route of points.
% The algo consists in rounding the corners of the broken line defined by
% the route points, using bezier curves with 3 points [start middle end].
% i.e. the lines [pt1 pt2 pt3], [pt3 pt4 pt5] and so on are represented by
% cubic bezier curves with two middle control points
%
% INPUTS
%---- xroute, yroute : vectors of the same length, x and y coordinates of
%                      the points defining the route
%---- Nptsbezier     : integer, number of points for each bezier curve
%---- alphastart     : scalar in [0,1], define the position of the start
%                      control point regarding the middle point position.
%                      0 means control point = start point, while 1 means
%                      control point = symmetrical of the middle point
%---- alphaend       : scalar in [0,1], define the position of the end
%                      control point regarding the middle point position.
%                      0 means control point = end point, while 1 means
%                      control point = middle point
%
% OUTPUTS
%---- xpath, ypath : row vectors, define the new smooth line representing
%                    the route

Nroute = length(xroute); % Nroute must be > 1

Nbezier = max([floor(Nroute/2) 1]); % nb of bezier curves

% Smooth path
xpath = zeros(1,Nbezier*Nptsbezier);
ypath = zeros(1,Nbezier*Nptsbezier);

% Smooth parameters
alpha1 = alphastart; % to define the control point associated with the end point
alpha2 = alphaend; % to define the control point associated with the start point

if Nroute == 2
    
    % Build an arc between the 2 points with low curvature
    istart = 1; % index of the start point
    iend = 2; % index of the end point
    d = sqrt((xroute(iend) - xroute(istart))^2 + (yroute(iend) - yroute(istart))^2);
    th = atan2(yroute(iend)-yroute(istart),xroute(iend)-xroute(istart));
    
    % Start and end middle control points: in the middle of the 2 points
    x1 = d/2;
    if 0 <= th && th < pi/2
        y1 = d/6;
    elseif pi/2 <= th && th < pi
        y1 = -d/6;
    elseif -pi <= th && th < -pi/2
        y1 = d/6;
    else
        y1 = -d/6;
    end
    xmid1 = xroute(istart) + x1*cos(th) - y1*sin(th);
    ymid1 = yroute(istart) + x1*sin(th) + y1*cos(th);
    xmid2 = xmid1;
    ymid2 = ymid1;
    
    % Set of four control points for the bezier curve
    x = [xroute(istart) xmid1 xmid2 xroute(iend)]';
    y = [yroute(istart) ymid1 ymid2 yroute(iend)]';
    
    % Build the bezier curve
    bc = beziercurve([x y],Nptsbezier);
    xpath = bc(:,1);
    ypath = bc(:,2);
    
else
    
    % Build a succesion of Bezier curves
    i = 1;
    for j = 1:Nbezier
        istart = i; % index of the start point
        imid = min([i+1 Nroute]); % index of the middle point
        iend = min([i+2 Nroute]); % index of the end point
        
        % Start middle control point
        if istart == 1
            xmid1 = xroute(istart) + alpha2*(xroute(imid) - xroute(istart));
            ymid1 = yroute(istart) + alpha2*(yroute(imid) - yroute(istart));
        else
            xmid1 = xroute(istart) + alpha1*(xroute(istart) - xroute(istart-1));
            ymid1 = yroute(istart) + alpha1*(yroute(istart) - yroute(istart-1));
        end
        
        % End middle control point
        xmid2 = xroute(iend) + alpha2*(xroute(imid) - xroute(iend));
        ymid2 = yroute(iend) + alpha2*(yroute(imid) - yroute(iend));
        
        % Set of four control points for the bezier curve
        x = [xroute(istart) xmid1 xmid2 xroute(iend)]';
        y = [yroute(istart) ymid1 ymid2 yroute(iend)]';
        
        % Build the bezier curve
        bc = beziercurve([x y],Nptsbezier);
        xpath(((j-1)*Nptsbezier+1):(j*Nptsbezier)) = bc(:,1);
        ypath(((j-1)*Nptsbezier+1):(j*Nptsbezier)) = bc(:,2);
        i = i + 2;
    end
end

end