function [intersection, t] = planeIntersectSegment(plane,p1,p2)

v = p2 - p1;
n = plane(1:3);

num = p1 * n.' + plane(4);
den = n * v.';
t = -num / den;

if t < 0 || t > 1
    intersection = [];
else
    intersection = p1 + t*v;
end

end