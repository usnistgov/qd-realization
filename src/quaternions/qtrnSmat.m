function S = qtrnSmat(q)
%QTRNSMAT helper function.
%
%   Copyright 2020 NIST/CLT (steve.blandino@nist.gov)

%#codegen


[r,c] = size(q);
assert(c==4, 'Expected input to be an array with 4 columns');
S = [zeros(r,1)   -q(:,4)   q(:,3); ...
    q(:,4)  zeros(r,1)     -q(:,2); ...
    -q(:,3)  q(:, 2)   zeros(r,1)];
S = reshape(reshape(S.', r*3, []).', 3,3,r);
end
