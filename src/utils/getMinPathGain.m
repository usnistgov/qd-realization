function minPathGain = getMinPathGain(txPos, rxPos, cadData, triangIdxList, freq, minDistMatrix)

triangList = cadData(triangIdxList, 1:9);

minDistTx2FirstTriang = simdTriPoint2(triangList(1, :), txPos);
minDistLastTriang2Rx = simdTriPoint2(triangList(end, :), rxPos);

rowIdxs = triangIdxList(1:end-1);
colIdxs = triangIdxList(2:end);
linearIdx = sub2ind(size(minDistMatrix), rowIdxs, colIdxs);
minDistBetweenTriangs = minDistMatrix(linearIdx);

totMinDist = minDistTx2FirstTriang +...
    sum(minDistBetweenTriangs) +...
    minDistLastTriang2Rx;
minPathGain = friisPathGain(totMinDist, freq);

end