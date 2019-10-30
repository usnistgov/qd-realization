function triangIdxs = generateReflectionList(previousReflectionList,...
    cadOp, visibilityMatrix)

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
    % index of the last reflecting triangle
    lastTriangIdx = previousReflectionList(i,end);
    % based on last reflection, list of next possible reflections
    newListTriang = find(visibilityMatrix(:, lastTriangIdx));
    % add new possible reflections to previousReflectionList(i,:)
    % using cell arrays for convenience, later merged into a single matrix
    newListLen = size(newListTriang,1);
    triangIdxs{i} = [repmat(previousReflectionList(i,:), newListLen, 1),...
        newListTriang];
end

% unwrap cell array to matrix
triangIdxs = cell2mat(triangIdxs);

end