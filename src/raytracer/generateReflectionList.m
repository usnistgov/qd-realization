function triangIdxs = generateReflectionList(previousReflectionList,cadOp)

nTriang = size(cadOp,1);
prevListLen = size(previousReflectionList,1);
% If starting from Tx (1st reflection), add all possible triangles to the
% list
if isempty(previousReflectionList)
    triangIdxs = (1:nTriang).';
    return
end

% else, start from last previous triangle
triangIdxs = cell(prevListLen,1);
for i = 1:prevListLen
    % index of the last triangle re
    lastTriangIdx = previousReflectionList(i,end);
    newListTriang = [1:lastTriangIdx-1, lastTriangIdx+1:nTriang].';
    newListLen = size(newListTriang,1);
    triangIdxs{i} = [repmat(previousReflectionList(i,:), newListLen, 1),...
        newListTriang];
end

% unwrap cell array to matrix
triangIdxs = cell2mat(triangIdxs);

end