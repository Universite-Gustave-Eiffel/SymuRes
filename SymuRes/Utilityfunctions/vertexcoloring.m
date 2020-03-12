function colorlist = vertexcoloring(Vertex,NumColors)
% colorlist = vertexcoloring(Vertex,NumColors)
% Return the list of color IDs assigned to the vertex of a graph, with the
% constraint that two adjacent vertex have not the same color if possible
%
% INPUTS
%---- Vertex    : cell array, list of vertex, Vertex{i} contains the
%                 adjacent Vertex of i
%---- NunColors : integer, total number of different colors, color IDs are
%                 then from 1 to NumColors
%
% OUTPUTS
%---- colorlist : vector, same size as Vertex, list of color IDs per vertex

Nvertex = length(Vertex);

colorlist = ones(1,Nvertex);

for i = 1:Nvertex
    colorID = mod(i,NumColors);
    if colorID == 0
        colorID = NumColors;
    end
    colorlist(i) = colorID;
end

if NumColors < Nvertex
    adjacentcolors = 1; % boolean which indicate that two adjacent vertex have the same assigned color
else
    adjacentcolors = 0; % if the number of colors is enough, no risk for same assigned colors
end
maxiter = 100; % maximum number of iterations allowed
iter = 0;

while adjacentcolors == 1 && iter < maxiter
    adjacentcolors = 0;
    for i = 1:Nvertex
        for j = Vertex{i} % loop on all adjacent vertex of vertex i
            if colorlist(j) == colorlist(i)
                adjacentcolors = 1;
                
                % Randomly sort another vertex
                vertexlist = 1:Nvertex;
                [i_, vertexlist] = find(vertexlist ~= i);
                [i_, vertexlist] = find(vertexlist ~= j);
                Nvertexlist = length(vertexlist);
                listindex = ceil(Nvertexlist*rand);
                k = vertexlist(listindex); % another vertex ID
                
                % Switch the colors assigned to vertex i and vertex k
                colori = colorlist(i);
                colork = colorlist(k);
                colorlist(i) = colork;
                colorlist(k) = colori;
            end
        end
    end
    iter = iter + 1;
end

if iter == maxiter
    warning('No perfect vertex coloring found')
end

end