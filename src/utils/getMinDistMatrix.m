function minDistMatrix = getMinDistMatrix(CADop)

nTriang = size(CADop, 1);
minDistMatrix = nan(nTriang, nTriang);

for i = 1:nTriang
    for j = i+1:nTriang
        minDist = simdTriTri2(CADop(i, 1:9), CADop(j, 1:9));
        minDistMatrix(i, j) = minDist;
        minDistMatrix(j, i) = minDist;
    end
end

end