function isObstructed = isRayObstructed(p1, p2, cadData, triangIdxs)
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

for i = 1:size(cadData, 1)
    if any(i == triangIdxs)
        % testing a triangle on which the bounce happened: not an
        % obstruction
        continue
    end
    
    planeEq = cadData(i, 10:13);
    vertex1 = cadData(i, 1:3);
    vertex2 = cadData(i, 4:6);
    vertex3 = cadData(i, 7:9);
    
    [intersection, t] = planeIntersectsSegment(planeEq, p1, p2, false);
    if t > 1e-9 && t < 1-1e-9 &&... % intersect, not just "touch"
            pointInTriangle(intersection, vertex1, vertex2, vertex3)
        isObstructed = true;
        return
    end
    
end

isObstructed = false;

end