function isObstructed = isRayObstructed(p1, p2, cadData, triangIdxs)

for i = 1:size(cadData,1)
    if any(i == triangIdxs)
        % testing a triangle on which the bounce happened: not an
        % obstruction
        continue
    end
    
    triang = cadData(i,:);
    planeEq = triang(10:13);
    vertex1 = triang(1:3);
    vertex2 = triang(4:6);
    vertex3 = triang(7:9);
    
    intersection = planeIntersectSegment(planeEq, p1, p2);
    if ~isempty(intersection) &&...
            PointInTriangle(intersection, vertex1, vertex2, vertex3)
        isObstructed = true;
        return
    end
    
end

isObstructed = false;

end