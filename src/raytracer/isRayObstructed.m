function isObstructed = isRayObstructed(p1, p2,...
    cadData, visibilityMatrix, triangIdxs)
%ISRAYOBSTRUCTED Check if the segment [p1,p2] is obstructed by any triangle
%contained in the CAD file. Partially robust to numeric approximations.


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

if isempty(triangIdxs)
    % No reflections from triangles, thus check all possible obstructions
    possibleObstructors = 1:size(cadData, 1);
else
    % Check obstructions only from triangles visible by all those involved
    % in the path segment
    possibleObstructors = find(all(visibilityMatrix(triangIdxs, :) > 0, 1));
end

for i = possibleObstructors
    
    planeEq = cadData(i, 10:13);
    vertex1 = cadData(i, 1:3);
    vertex2 = cadData(i, 4:6);
    vertex3 = cadData(i, 7:9);
    
    [intersection, t] = planeIntersectsSegment(planeEq, p1, p2, false);
    if ~isempty(t) && ... % valid intersection with plane
            t > 1e-9 && t < 1-1e-9 &&... % intersect, not just "touch"
            pointInTriangle(intersection, vertex1, vertex2, vertex3)
        isObstructed = true;
        return
    end
    
end

isObstructed = false;

end