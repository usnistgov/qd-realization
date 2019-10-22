function [output, multipath] = mymultipath(rxPos, txPos, rxVel, txVel, triangIdxList,...
    cadData, materialLibrary, switchQd, switchMaterial, freq)

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
    
    [pathExists, dod, doa, rowMultipath, rayLen, dopplerFactor, pathGain]...
        = computeSingleRay(txPos, rxPos, txVel, rxVel, currentTriangIdxs,...
        cadData, materialLibrary, switchMaterial, freq);
    
    if pathExists
        deterministicRayOutput = fillOutput(reflOrder, dod, doa, rayLen,...
            pathGain - 10, dopplerFactor, freq);
        
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