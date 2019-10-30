function [output, multipath, currentMaxPathGain] = mymultipath(rxPos, txPos,...
    rxVel, txVel, triangIdxList, cadData, visibilityMatrix,...
    materialLibrary, switchQd, switchMaterial, freq,...
    minAbsolutePathGainThreshold, minRelativePathGainThreshold, currentMaxPathGain)
%MYMULTIPATH Computes ray, if it exists, between rxPos and txPos, bouncing
%over the given lists of triangles. Also computes QD MPCs if requested.


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

triangListLen = size(triangIdxList,1);
reflOrder = size(triangIdxList,2);

% init outputs
if switchQd
    % concatenate D-rays and QD MPCs
    output = nan(0, 21);
else
    % Only D-rays
    output = nan(triangListLen, 21);
end
multipath = nan(triangListLen, 1 + 3 * (reflOrder+2));

multipath(:,1) = reflOrder;

for i = 1:triangListLen
    currentTriangIdxs = triangIdxList(i,:);
    
    [pathExists, dod, doa, rowMultipath, rayLen, dopplerFactor, pathGain,...
        currentMaxPathGain] = computeSingleRay(txPos, rxPos, txVel, rxVel,...
        currentTriangIdxs, cadData, visibilityMatrix, materialLibrary,...
        switchQd, switchMaterial, freq, minAbsolutePathGainThreshold,...
        minRelativePathGainThreshold, currentMaxPathGain);
    
    if pathExists
        deterministicRayOutput = fillOutput(reflOrder, dod, doa, rayLen,...
            pathGain, dopplerFactor, freq);
        
        if switchQd
            % TODO: update QD model
            arrayOfMaterials = materialLibrary.Material(triangIdxList);
            qdOutput = [];%QDGenerator(reflOrder, deterministicRayOutput, arrayOfMaterials, iterateNumberOfRowsArraysOfPlanes?,...
                %materialLibrary, rayLen, freq, indexOutput?,...
                %dod, doa, txVel, velocityTemp?, indexMultipath?, indexReference?);
            output = [output; deterministicRayOutput; qdOutput];
            
        else
            output(i,:) = deterministicRayOutput;
            
        end
        
        multipath(i,2:end) = rowMultipath;
    end
end

% Remove invalid output rows
invalidRows = all(isnan(output), 2);
output(invalidRows,:) = [];
multipath(invalidRows,:) = [];

end