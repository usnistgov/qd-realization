function eout = eulFqtrn(q ,rot)
%QPARTS2FEUL - Euler angles from quaternion parts
%   This function is for internal use only. It may be removed in the future.

%   Copyright 2017 The MathWorks, Inc.

%#codegen

%column-ize quaternion parts
qa = q(:,1);
qb = q(:,2);
qc = q(:,3);
qd = q(:,4);
the1 = ones(size(qa), 'like', qa);
the2 = 2*the1;

found = true;
switch upper(rot)
    case 'YZY'
        tmp = qa.^2.*the2 - the1 + qc.^2.*the2;
        tmp(tmp > the1(1)) = the1(1);
        tmp(tmp < -the1(1)) = -the1(1);
        b = acos(tmp);
        a = atan2((qa.*qb.*the2 + qc.*qd.*the2),(qa.*qd.*the2 - qb.*qc.*the2));
        c = -atan2((qa.*qb.*the2 - qc.*qd.*the2),(qa.*qd.*the2 + qb.*qc.*the2));
    case 'YXY'
        tmp = qa.^2.*the2 - the1 + qc.^2.*the2;
        tmp(tmp > the1(1)) = the1(1);
        tmp(tmp < -the1(1)) = -the1(1);
        b = acos(tmp);
        a = -atan2((qa.*qd.*the2 - qb.*qc.*the2),(qa.*qb.*the2 + qc.*qd.*the2));
        c = atan2((qa.*qd.*the2 + qb.*qc.*the2),(qa.*qb.*the2 - qc.*qd.*the2));
    case 'ZYZ'
        tmp = qa.^2.*the2 - the1 + qd.^2.*the2;
        tmp(tmp > the1(1)) = the1(1);
        tmp(tmp < -the1(1)) = -the1(1);
        b = acos(tmp);
        a = -atan2((qa.*qb.*the2 - qc.*qd.*the2),(qa.*qc.*the2 + qb.*qd.*the2));
        c = atan2((qa.*qb.*the2 + qc.*qd.*the2),(qa.*qc.*the2 - qb.*qd.*the2));
    case 'ZXZ'
        tmp = qa.^2.*the2 - the1 + qd.^2.*the2;
        tmp(tmp > the1(1)) = the1(1);
        tmp(tmp < -the1(1)) = -the1(1);
        b = acos(tmp);
        a = atan2((qa.*qc.*the2 + qb.*qd.*the2),(qa.*qb.*the2 - qc.*qd.*the2));
        c = -atan2((qa.*qc.*the2 - qb.*qd.*the2),(qa.*qb.*the2 + qc.*qd.*the2));
    case 'XYX'
        tmp = qa.^2.*the2 - the1 + qb.^2.*the2;
        tmp(tmp > the1(1)) = the1(1);
        tmp(tmp < -the1(1)) = -the1(1);
        b = acos(tmp);
        a = atan2((qa.*qd.*the2 + qb.*qc.*the2),(qa.*qc.*the2 - qb.*qd.*the2));
        c = -atan2((qa.*qd.*the2 - qb.*qc.*the2),(qa.*qc.*the2 + qb.*qd.*the2));
    case 'XZX'
        tmp = qa.^2.*the2 - the1 + qb.^2.*the2;
        tmp(tmp > the1(1)) = the1(1);
        tmp(tmp < -the1(1)) = -the1(1);
        b = acos(tmp);
        a = -atan2((qa.*qc.*the2 - qb.*qd.*the2),(qa.*qd.*the2 + qb.*qc.*the2));
        c = atan2((qa.*qc.*the2 + qb.*qd.*the2),(qa.*qd.*the2 - qb.*qc.*the2));
    case 'XYZ'
        tmp = qa.*qc.*the2 + qb.*qd.*the2;
        tmp(tmp > the1(1)) = the1(1);
        tmp(tmp < -the1(1)) = -the1(1);
        b = asin(tmp);
        a = atan2((qa.*qb.*the2 - qc.*qd.*the2),(qa.^2.*the2 - the1 + qd.^2.*the2));
        c = atan2((qa.*qd.*the2 - qb.*qc.*the2),(qa.^2.*the2 - the1 + qb.^2.*the2));
    case 'YZX'
        tmp = qa.*qd.*the2 + qb.*qc.*the2;
        tmp(tmp > the1(1)) = the1(1);
        tmp(tmp < -the1(1)) = -the1(1);
        b = asin(tmp);
        a = atan2((qa.*qc.*the2 - qb.*qd.*the2),(qa.^2.*the2 - the1 + qb.^2.*the2));
        c = atan2((qa.*qb.*the2 - qc.*qd.*the2),(qa.^2.*the2 - the1 + qc.^2.*the2));
    case 'ZXY'
        tmp = qa.*qb.*the2 + qc.*qd.*the2;
        tmp(tmp > the1(1)) = the1(1);
        tmp(tmp < -the1(1)) = -the1(1);
        b = asin(tmp);
        a = atan2((qa.*qd.*the2 - qb.*qc.*the2),(qa.^2.*the2 - the1 + qc.^2.*the2));
        c = atan2((qa.*qc.*the2 - qb.*qd.*the2),(qa.^2.*the2 - the1 + qd.^2.*the2));
    case 'XZY'
        tmp = qb.*qc.*the2 - qa.*qd.*the2;
        tmp(tmp > the1(1)) = the1(1);
        tmp(tmp < -the1(1)) = -the1(1);
        b = -asin(tmp);
        a = atan2((qa.*qb.*the2 + qc.*qd.*the2),(qa.^2.*the2 - the1 + qc.^2.*the2));
        c = atan2((qa.*qc.*the2 + qb.*qd.*the2),(qa.^2.*the2 - the1 + qb.^2.*the2));
    case 'ZYX'
        tmp = qb.*qd.*the2 - qa.*qc.*the2;
        tmp(tmp > the1(1)) = the1(1);
        tmp(tmp < -the1(1)) = -the1(1);
        b = -asin(tmp);
        a = atan2((qa.*qd.*the2 + qb.*qc.*the2),(qa.^2.*the2 - the1 + qb.^2.*the2));
        c = atan2((qa.*qb.*the2 + qc.*qd.*the2),(qa.^2.*the2 - the1 + qd.^2.*the2));
    case 'YXZ'
        tmp = qc.*qd.*the2 - qa.*qb.*the2;
        tmp(tmp > the1(1)) = the1(1);
        tmp(tmp < -the1(1)) = -the1(1);
        b = -asin(tmp);
        a = atan2((qa.*qc.*the2 + qb.*qd.*the2),(qa.^2.*the2 - the1 + qd.^2.*the2));
        c = atan2((qa.*qd.*the2 + qb.*qc.*the2),(qa.^2.*the2 - the1 + qc.^2.*the2));
    otherwise
        found = false;
        a = zeros(size(qa), 'like', qa);
        b = zeros(size(qa), 'like', qa);
        c = zeros(size(qa), 'like', qa);
end
coder.internal.assert(found, 'shared_rotations:quaternion:NoSeqConv', rot);
eout = [a b c];
end