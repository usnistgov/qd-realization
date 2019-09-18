function [intersection, t] = planeIntersectsSegment(plane,pStart,pEnd,checkDirection)

v = pEnd - pStart;
n = plane(1:3);

num = pStart * n.' + plane(4);
den = n * v.';

if checkDirection && den < 1e-9
    intersection = [];
    t = [];
else
    t = -num / den;
    
    if t < 0 || t > 1
        intersection = [];
    else
        intersection = pStart + t*v;
    end
end

end