function isInFront = isTriangleInFrontOfTriangle(t1, t2, strict)

points = [t1(1:3);...
    t1(4:6);...
    t1(7:9)];
plane = t2(10:13);

dist = signedDistanceOfPointFromPlane(points, plane);
if strict
    isInFront = any(dist > 0);
else
    isInFront = any(dist >= 0);
end

end