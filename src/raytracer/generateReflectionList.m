function triangIdxs = generateReflectionList(previousReflectionList,...
    cadOp, visibilityMatrix)
%GENERATEREFLECTIONLIST Given the the reflection list of the lower order,
%generates the next order reflection list.


% Copyright (c) 2019, University of Padova, Department of Information
% Engineering, SIGNET lab.
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%    http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

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