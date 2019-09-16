function [intersection, t] = planeIntersectSegment(plane,p1,p2)

v = p2 - p1;
n = plane(1:3);

t = -(dot(p1,n) + plane(4)) / dot(n,v);

if t < 0 || t > 1
    intersection = [];
else
    intersection = p1 + t*v;
end

end