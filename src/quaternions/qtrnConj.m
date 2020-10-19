function Q =  qtrnConj(Q)
%QTRNCONJ Quaternion conjugate operation
%   P = QTRNCONJ(Q) given a quaternion Q returns the quaternion P = Q*
%
%   Copyright 2020 NIST/CTL (steve.blandino@nist.gov)

%#codegen

assert(size(Q, 2) == 4, 'Expected input to be an array with 4 columns');

Q(:, 2:end) = -Q(:, 2:end);
end
