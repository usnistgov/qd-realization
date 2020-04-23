function R = qtrnRotmat(q, f)
%QTRNROTMAT converts quaternion in rotation matrix
%   R = QTRNROTMAT(Q) returns the rotation matrix R for point
%   rotation equivalent to the quaternion Q
%
%   R = QTRNROTMAT(Q, 'frame') returns the rotation matrix R for frame
%   rotation equivalent to the quaternion Q
%
%
%   Copyright 2020 NIST/CLT (steve.blandino@nist.gov)

%#codegen

if ~exist('f', 'var')
    f = 'point';
end
N = size(q,1);
R = zeros(3,3,N);

q = qtrnNormalize(q);
for i = 1: N
    u0 = q(i,1);
    u =  q(1,2:end);    
    R(:,:, i) = (u0^2 - u*u')*eye(3) + 2*(u'*u) + 2*u0* qtrnSmat(q(1,:));
end

if strcmpi(f, 'frame')
    R  = permute(R, [2 1 3]);
end
end