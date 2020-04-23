function P= qtrnVectorrotate(q,v)
%QTRNVECTORROTATE returns the rotated vector.
%
%   P = QTRNVECTORROTATE(Q,V) calculates the rotated vector P obtained by 
%   rotating the vector v by a quaternion q.
%
%
%   Copyright 2020 NIST/CLT (steve.blandino@nist.gov)    
    
P = qtrnMultiply(qtrnMultiply(q,v),qtrnConj(q));

end