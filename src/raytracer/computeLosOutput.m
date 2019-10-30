function [output, multipath] = computeLosOutput(rxPos, txPos, rxVel, txVel,...
    cadData, freq)

isObstructed = isRayObstructed(rxPos, txPos, cadData, [], []);
if isObstructed
   output = [];
   multipath = [];
   
else
    reflOrder = 0;
    dod = rxPos - txPos;
    doa = -dod;
    rayLen = norm(dod);
    pathGain = friisPathGain(rayLen, freq);
    dopplerFactor = getDopplerFactor(txPos, rxPos, txVel, rxVel, [], []);
    
    output = fillOutput(reflOrder, dod, doa, rayLen, pathGain, dopplerFactor, freq);
    multipath = [rxPos, txPos];
end

end