function distance = signedDistanceOfPointFromPlane(point, plane)
distance = ((point * plane(1:3).') + plane(4)) / norm(plane(1:3));
end