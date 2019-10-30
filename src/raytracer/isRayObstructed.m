function isObstructed = isRayObstructed(p1, p2,...
    cadData, visibilityMatrix, triangIdxs)

if isempty(triangIdxs)
    % No reflections from triangles, thus check all possible obstructions
    possibleObstructors = 1:size(cadData, 1);
else
    % Check obstructions only from triangles visible by all those involved
    % in the path segment
    possibleObstructors = find(all(visibilityMatrix(triangIdxs, :) > 0, 1));
end

for i = possibleObstructors
    
    planeEq = cadData(i,10:13);
    vertex1 = cadData(i,1:3);
    vertex2 = cadData(i,4:6);
    vertex3 = cadData(i,7:9);
    
    [intersection,t] = planeIntersectsSegment(planeEq, p1, p2, false);
    if t > 1e-9 && t < 1-1e-9 &&... % intersect, not just "touch"
            pointInTriangle(intersection, vertex1, vertex2, vertex3)
        isObstructed = true;
        return
    end
    
end

isObstructed = false;

end