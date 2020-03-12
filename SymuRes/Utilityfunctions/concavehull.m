function [xcontour, ycontour] = concavehull(xdata,ydata,NumNeighbors)
% [xcontour, ycontour] = concavehull(xdata,ydata,NumNeighbors)
% Return the concave hull of a set of points
% Based on the k-nearest points algorithm of Moreira & Santos (2007)
% http://repositorium.sdum.uminho.pt/bitstream/1822/6429/1/ConcaveHull_ACM_MYS.pdf
%
% INPUTS
%---- xdata        : vector, x-values of the set of points
%---- ydata        : vector, same size as xdata, y-values
%---- NumNeighbors : integer >= 3, indicate the number of nearest neighbors
%                    to explore the next point in the hull
%
% OUTPUTS
%---- xcontour : vector, x-values of the points that define the hull
%---- ycontour : vector, same size as xcontour, y-values

NumNeighbors = max([NumNeighbors 3]);

Ndata = length(xdata);

indexhull = []; % hull set as indexes of points

if Ndata > 3
    NumNeighbors = min([NumNeighbors Ndata-1]);
    indexdata = 1:Ndata; % indexes of points
    
    [y, minindex] = min(ydata);
    ifirst = minindex; % first point of the hull
    indexhull = [indexhull ifirst];
    icurrent = ifirst;
    indexdata(ifirst) = []; % remove the first point of the original set
    prevang = pi;
    k = 2;
    while (icurrent ~= ifirst || k == 2) && ~isempty(indexdata)
        if k == 5
            indexdata = [indexdata ifirst];
        end
        listdist = sqrt((xdata(icurrent) - xdata(indexdata)).^2 + (ydata(icurrent) - ydata(indexdata)).^2);
        [listdist, sortindex] = sort(listdist);
        num = min([NumNeighbors length(indexdata)]);
        indexneighbors = indexdata(sortindex(1:num));
        currang = atan2(ydata(indexneighbors)-ydata(icurrent),xdata(indexneighbors)-xdata(icurrent));
        currang = currang + (sign(currang) == -1).*2*pi;
        prevang = prevang + (sign(prevang) == -1).*2*pi;
        listang = prevang - currang;
        ind = 1:length(listang);
        ind1 = ind(listang >= 0);
        ind2 = ind(listang < 0);
        listang1 = listang(listang >= 0);
        listang2 = listang(listang < 0);
        [listang1, sortindex1] = sort(listang1,'descend');
        [listang2, sortindex2] = sort(listang2,'descend');
        sortindex = [ind2(sortindex2) ind1(sortindex1)];
        %[listang, sortindex] = sort(listang,'descend');
        indexcandidates = indexneighbors(sortindex);
        
%         clf
%         plot(xdata(indexdata),ydata(indexdata),'ok',xdata(indexneighbors),ydata(indexneighbors),'ob',...
%             xdata(indexcandidates(1)),ydata(indexcandidates(1)),'og',xdata(indexhull),ydata(indexhull),'-xr')
%         pause
        
        its = 1;
        i = 0;
        while its == 1 && i < length(indexcandidates)
            i = i + 1;
            if indexcandidates(i) == ifirst
                ilast = 1;
            else
                ilast = 0;
            end
            its = 0;
            j = 2;
            while its == 0 && j < length(indexhull) - ilast
                x11 = xdata(indexhull(k-1));
                y11 = ydata(indexhull(k-1));
                x12 = xdata(indexcandidates(i));
                y12 = ydata(indexcandidates(i));
                
                x21 = xdata(indexhull(k-1-j));
                y21 = ydata(indexhull(k-1-j));
                x22 = xdata(indexhull(k-j));
                y22 = ydata(indexhull(k-j));
                
                % Test if [pt11 pt12] intersects [pt21 pt22]
                its = lineintersect([x11 y11],[x12 y12],[x21 y21],[x22 y22]);
                
                j = j + 1;
            end
        end
        if its == 1
            warning('Since all candidates intersect at least one edge, try again with a higher number of neighbors')
            %concavehull(xdata,ydata,NumNeighbors+1)
        end
        
        %i = 1;
        icurrent = indexcandidates(i);
        indexhull = [indexhull icurrent];
        %prevang = atan2(ydata(indexhull(k))-ydata(indexhull(k-1)),xdata(indexhull(k))-xdata(indexhull(k-1)));
        prevang = atan2(ydata(indexhull(k-1))-ydata(indexhull(k)),xdata(indexhull(k-1))-xdata(indexhull(k)));
        indexdata(indexdata == icurrent) = [];
        
        k = k + 1;
    end
    allinside = 1;
    i = length(indexdata);
    while allinside == 1 && i > 0
        xi = xdata(indexdata(i));
        yi = ydata(indexdata(i));
        % Test if point i is inside the concave hull
        allinside = 1;
        for j = 1:(length(indexhull)-1)
            xj1 = xdata(indexhull(j));
            yj1 = ydata(indexhull(j));
            xj2 = xdata(indexhull(j+1));
            yj2 = ydata(indexhull(j+1));
            prodvec = (xi-xj1)*(yj2-yj1) - (yi-yj1)*(xj2-xj1);
            if prodvec >= 0
                % If at least one cross product is > 0, the point is
                % outside the hull
                allinside = 0;
            end
        end
        i = i - 1;
    end
    if allinside == 0
        warning('Since at least one point is out of the computed polygon, try again with a higher number of neighbors')
        %concavehull(xdata,ydata,NumNeighbors+1)
    end
    xcontour = xdata(indexhull);
    ycontour = ydata(indexhull);
else
    xcontour = xdata;
    ycontour = ydata;
end

end