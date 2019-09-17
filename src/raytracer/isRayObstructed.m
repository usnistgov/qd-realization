function isObstructed = isRayObstructed(p1, p2, cadData, triangIdxs)

for i = 1:size(cadData,1)
    if any(i == triangIdxs)
        % testing a triangle on which the bounce happened: not an
        % obstruction
        continue
    end
    
    planeEq = cadData(i,10:13);
    vertex1 = cadData(i,1:3);
    vertex2 = cadData(i,4:6);
    vertex3 = cadData(i,7:9);
    
    [intersection,t] = planeIntersectSegment(planeEq, p1, p2);
    if t > 1e-9 && t < 1-1e-9 &&... % intersect, not just "touch"
            pointInTriangle(intersection, vertex1, vertex2, vertex3)
        isObstructed = true;
        return
    end
    
end

isObstructed = false;

end