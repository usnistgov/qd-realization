function [exists, dod, doa, multipath, rayLength, dopplerFactor, pathGain]...
    = computeSingleRay(txPos, rxPos, txVel, rxVel, triangIdxList,...
    cadData, materialLibrary, switchQd, switchMaterial, freq)

intersections = methodOfImages(txPos, rxPos, cadData, triangIdxList, 1);

if isempty(intersections)
    % do something
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
    reflectionLosses = sum(materialLibrary.mu_RL(materialsList));
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


function dopplerFactor = getDopplerFactor(txPos, rxPos, txVel, rxVel,...
    cadData, triangIdxList)

txVelRefl = txVel;
txPosRefl = txPos;
for i = 1:length(triangIdxList)
    plane = cadData(triangIdxList(i), 10:13);
    txVelRefl = reflectedVelocity(txVelRefl, plane);
    txPosRefl = reflectedImagePointPlane(txPosRefl, plane);
end

relativeVel = txVelRefl - rxVel;
relativePos = txPosRefl - rxPos;
radialRelativeVel = dot(relativeVel, relativePos) / norm(relativePos);

dopplerFactor = -radialRelativeVel / 3e8;

end