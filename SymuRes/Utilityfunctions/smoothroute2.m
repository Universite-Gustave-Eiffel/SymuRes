function [xpath, ypath] = smoothroute2(xroute,yroute,Nptsbezier,tension)
% [xpath, ypath] = smoothroute2(xroute,yroute,Nptsbezier,tension)
% Return a smooth path representing a given route of points using connected
% Bezier curves. The definition of the Bezier control points comes from the
% code of Rob Spencer.
% http://scaledinnovation.com/analytics/splines/aboutSplines.html
%
% INPUTS
%---- xroute, yroute : vectors of the same length, x and y coordinates of
%                      the points defining the route
%---- Nptsbezier     : integer, number of points for each bezier curve
%---- tension        : scalar, define the smoothness of the curve, can be
%                      negative, tension = 0 means straight lines between
%                      points
%
% OUTPUTS
%---- xpath, ypath : row vectors, define the new smooth line representing
%                    the route

Nroute = length(xroute); % Nroute must be > 1

Nbezier = Nroute - 1; % nb of bezier curves

% Smooth path
xpath = zeros(1,Nbezier*Nptsbezier);
ypath = zeros(1,Nbezier*Nptsbezier);

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
    xmid1 = xroute(1);
    ymid1 = yroute(1);
    for i = 1:Nbezier
        istart = i; % index of the start point
        iend = i + 1; % index of the end point
        
        if i < Nbezier
            inext = i + 2;
            d2 = sqrt((xroute(iend) - xroute(istart))^2 + (yroute(iend) - yroute(istart))^2);
            d3 = sqrt((xroute(inext) - xroute(iend))^2 + (yroute(inext) - yroute(iend))^2);
            alpha2 = tension*d2/(d2 + d3);
            alpha3 = tension*d3/(d2 + d3);
            % End control point
            xmid2 = xroute(iend) - alpha2*(xroute(inext) - xroute(istart));
            ymid2 = yroute(iend) - alpha2*(yroute(inext) - yroute(istart));
            % Next start control point
            xmid3 = xroute(iend) + alpha3*(xroute(inext) - xroute(istart));
            ymid3 = yroute(iend) + alpha3*(yroute(inext) - yroute(istart));
        else
            % End control point
            xmid2 = xroute(iend);
            ymid2 = yroute(iend);
        end
        
        % Set of four control points for the bezier curve
        x = [xroute(istart) xmid1 xmid2 xroute(iend)]';
        y = [yroute(istart) ymid1 ymid2 yroute(iend)]';
        
        % Build the bezier curve
        bc = beziercurve([x y],Nptsbezier);
        xpath(((i-1)*Nptsbezier+1):(i*Nptsbezier)) = bc(:,1);
        ypath(((i-1)*Nptsbezier+1):(i*Nptsbezier)) = bc(:,2);
        
        % Next start control point
        xmid1 = xmid3;
        ymid1 = ymid3;
    end
end

end