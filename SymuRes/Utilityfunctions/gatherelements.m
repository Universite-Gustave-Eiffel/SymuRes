function indexlist = gatherelements(similarityMat)
% indexlist = gatherelements(similarityMat)
% Gathers N elements in lists by similarity, the similarity being given by
% a N-by-N symmetric matrix indicating M(i,j) = 1 when elements i and j are
% considered similar, and 0 elsewhere
%
% INPUTS
%---- similarityMat : N-by-N matrix, the similarity matrix
%
% OUTPUTS
%---- indexlist : Ni-size cell, gathering the Ni groups of similar elements
%                 (Ni lists of indexes)

Nelem = size(similarityMat,1);
indexlist = cell(1,Nelem);
ilist = 0;
isnewlist = 1;
for i = 1:Nelem
    if ilist > 0
        isnewlist = 1;
        for i2 = 1:ilist
            if ismember(i,indexlist{i2})
                % if not already in a list, open a new list
                isnewlist = 0;
            end
        end
    end
    if isnewlist == 1
        % if a new list has been opened, search for similar paths
        ilist = ilist + 1;
        indexlist{ilist} = [indexlist{ilist} i];
        if i < Nelem
            for i2 = (i+1):Nelem
                if similarityMat(i,i2) == 1
                    % add to the list
                    indexlist{ilist} = [indexlist{ilist} i2];
                end
            end
        end
    end
end

indexlist = {indexlist{1:ilist}};

end