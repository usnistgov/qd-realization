function q = qtrnNormalize(q)
%QTRNNORMALIZE returns the normalized quaternion
%   P = QTRNNORMALIZE(Q) returns P the normalized quaternion Q: NORM of P
%   is 1.
%   Q is a 4 dimensional vector or a Nx4 matrix collection N quaternions.
%
%   Copyright 2020 NIST/CTL (steve.blandino@nist.gov)

%#codegen

if isvector(q)
    assert(length(q)==4, 'Expected input to be a vector with 4 elements');
    q = (q(:)./sqrt((q(:)'*q(:)))).';
else
    [~,c] = size(q);
    assert(c==4, 'Expected input to be an array with 4 columns');
    q = q./sqrt(sum(q.*q,2));
end
end