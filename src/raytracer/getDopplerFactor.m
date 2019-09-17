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