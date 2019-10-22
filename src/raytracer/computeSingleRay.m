function [exists, dod, doa, multipath, rayLength, dopplerFactor, pathGain]...
    = computeSingleRay(txPos, rxPos, txVel, rxVel, triangIdxList,...
    cadData, materialLibrary, switchMaterial, freq)
%COMPUTESINGLERAY Computes geometry and physics of a ray between txPos and
%rxPos, bouncing over a give list of triangles


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

intersections = methodOfImages(txPos, rxPos, cadData, triangIdxList, 1);

if isempty(intersections)
    % the ray does not exist
    exists = false;
    dod = NaN;
    doa = NaN;
    multipath = NaN;
    rayLength = NaN;
    dopplerFactor = NaN;
    pathGain = NaN;
    return
end

% check if ray exists
exists = verifyRayExists(txPos, intersections, rxPos, cadData, triangIdxList);

% Extract info
dod = intersections(1,:) - txPos;
doa = intersections(end,:) - rxPos;
multipath = getMultipathVector(txPos, intersections, rxPos);
rayLength = getRayLength(txPos, intersections, rxPos);

friisPg = friisPathGain(rayLength,freq);
if switchMaterial
    materialsList = cadData(triangIdxList, 14);
    reflectionLosses = sum(materialLibrary.mu_RL(materialsList)); % TODO: update with Rician distribution
else
    reflectionLosses = 0;
end
pathGain = friisPg - reflectionLosses;

dopplerFactor = getDopplerFactor(txPos, rxPos, txVel, rxVel, cadData, triangIdxList);

end


%% Utils
function exists = verifyRayExists(txPos, intersections, rxPos, cadData, triangIdxList)
points = [txPos; intersections; rxPos];
nRays = size(points,1) - 1;

for i = 1:nRays
    p1 = points(i,:);
    p2 = points(i+1,:);
    
    if i == 1
        % from tx to first bounce
        bounceTriangIdx = triangIdxList(1);
    elseif i == nRays
        % from last bounce to rx
        bounceTriangIdx = triangIdxList(end);
    else
        % from bounce to bounce starting from i=2
        bounceTriangIdx = triangIdxList([i-1, i]);
    end
    
    isObstructed = isRayObstructed(p1,p2,cadData,bounceTriangIdx);
    if isObstructed
        % Early stopping if ray does not exist
        exists = false;
        return
    end
end

exists = true;

end


function multipath = getMultipathVector(txPos, intersections, rxPos)
points = [txPos; intersections; rxPos];
points = flipud(points); % rx first, tx last
points = points.'; % following MATLAB's natural matrix indexing
multipath = points(:).';
end


function l = getRayLength(txPos, intersections, rxPos)
vectors = [intersections; rxPos] - [txPos; intersections];
l = sum(vecnorm(vectors,2,2));
end