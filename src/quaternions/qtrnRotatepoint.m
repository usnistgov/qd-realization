function P = qtrnRotatepoint(q,v)
%QTRNROTPOINT quaternion point roatation.
%   P = QTRNROTPOINT(Q,V) returns the rotate points U from the points in V
%   using the quaternion Q.
%
%
%   Copyright 2020 NIST/CLT (steve.blandino@nist.gov)

%#codegen


q = qtrnNormalize(q);
N = size(v,1);

v = [zeros(N,1), v];
P = qtrnVectorrotate(q,v);
P = P(:,2:end);


end