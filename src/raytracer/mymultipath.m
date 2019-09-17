function [output, multipath] = mymultipath(rxPos, txPos, rxVel, txVel, triangIdxList,...
    cadData, materialLibrary, switchQd, switchMaterial, freq)

triangListLen = size(triangIdxList,1);
reflOrder = size(triangIdxList,2);

output = nan(triangListLen, 21);
multipath = nan(triangListLen, 1 + 3 * (reflOrder+2));

multipath(:,1) = reflOrder;

for i = 1:triangListLen
    currentTriangIdxs = triangIdxList(i,:);
    
    [pathExists, dod, doa, rowMultipath, rayLen, dopplerFactor, pathGain]...
    = computeSingleRay(txPos, rxPos, txVel, rxVel, currentTriangIdxs,...
    cadData, materialLibrary, switchQd, switchMaterial, freq);
    
    if pathExists
        output(i,:) = fillOutput(reflOrder, dod, doa, rayLen,...
            pathGain - 10, dopplerFactor, freq);
        multipath(i,2:end) = rowMultipath;
    end
end

% Remove invalid output rows
invalidRows = all(isnan(output), 2);
output(invalidRows,:) = [];
multipath(invalidRows,:) = [];

end