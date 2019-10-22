function [intersection, t] = planeIntersectsSegment(plane,pStart,pEnd,checkDirection)
%PLANEINTERSECTSSEGMENT Obtain intersection between plane described by its
%equation and a segmet described by it end points. Optionally check if
%pStart->pEnd follows the same direction as the normal vector given by the
%plane equation.


% Copyright (c) 2019, University of Padova, Department of Information
% Engineering, SIGNET lab.
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%    http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

% segment vector and plane normal
v = pEnd - pStart;
n = plane(1:3);

% numerator and denominator of intersection equation
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