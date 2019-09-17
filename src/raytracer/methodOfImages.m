function intersections = methodOfImages(txPos,rxPos,cadData,triangIdxs,recursionDepth)
triang = cadData(triangIdxs(recursionDepth),:);
plane = triang(10:13);
txPosReflected = reflectedImagePointPlane(txPos, plane);

reflOrder = length(triangIdxs);

if reflOrder == recursionDepth
    prevIntersection = rxPos;
    intersections = nan(reflOrder, 3);
    
else % need more reflections
    intersections = methodOfImages(txPosReflected,rxPos,cadData,triangIdxs,recursionDepth+1);
    if isempty(intersections)
        return
    else
        prevIntersection = intersections(recursionDepth+1,:);
    end
end

intersection = planeIntersectSegment(plane,txPosReflected,prevIntersection);
if isempty(intersection)
    intersections = [];
    return
end
% else: valid intersection with triangle's plane
isInsideTriangle = pointInTriangle(intersection, triang(1:3), triang(4:6), triang(7:9));
if isInsideTriangle
    intersections(recursionDepth, :) = intersection;
else
    intersections = [];
    return
end

end