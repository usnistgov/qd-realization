function output = mymultipath(rxPos, txPos, rxVel, txVel, triangIdxList,...
    cadData, MaterialLibrary, switchQd, switchMaterial, freq)

triangListLen = size(triangIdxList,1);
reflOrder = size(triangIdxList,2);

output = nan(triangListLen, 21);

multipath = nan(triangListLen, 1 + 3 * (reflOrder+2));
multipath(:,1) = reflOrder;
multipath(:,2:4) = repmat(rxPos,triangListLen,1); % TODO: maybe not necessary

for i = 1:triangListLen
    currentTriangIdxs = triangIdxList(i,:);
    % Extract variables for singleMultipathGenerator
    arrayOfPlanes = getArrayOfPlanes(currentTriangIdxs,cadData);
    arrayOfPoints = getArrayOfPoints(currentTriangIdxs,cadData);
    
    reflected = rxPos; % To start recursion. Rx is reflected, not Tx
    
    [isValidPath, ~, dod, doa, multipath, rayLen, dopplerFactor,...
    extraPathLoss, ~, ~, ~, ~, ~] = singleMultipathGenerator(...
        1, reflOrder, 1, arrayOfPlanes, arrayOfPoints,...
        reflected, rxPos, txPos, cadData, multipath, i, txVel, rxVel,...
        [], [], [], [], [], [], []);
    
    computeSingleRay(txPos, rxPos, txVel, rxVel, triangIdxList(i,:),...
    cadData, MaterialLibrary, switchQd, switchMaterial, freq);
    
    if isValidPath == 1
        output(i,:) = fillOutput(reflOrder, dod, doa, rayLen,...
            -extraPathLoss, dopplerFactor, freq);
    end
end

% Remove invalid output rows
invalidRows = all(isnan(output), 2);
output(invalidRows,:) = [];

end