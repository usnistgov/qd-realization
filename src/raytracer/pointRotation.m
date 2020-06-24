function [P, varargout] = pointRotation(P, Pc, eucl, varargin)
%POINTROTATION rotation in QD software.
%   P = POINTROTATION(P,PC,EUCL) rotates the point P with respect to a
%   reference system centered in PC by the euclidians angles EUCL. Default
%   euclidian sequence is assume 'ZXY'.
%  
%  P = POINTROTATION(P,PC,EUCL, 1) returns the coordinates of the point P
%  when the reference system centered in PC is rotated by the
%  euclidians angles EUCL. Default euclidian sequence is assume 'ZXY'.
% 
%  [P Q] = POINTROTATION(..) returns the quaternion used to obtain the
%  rotated position P.
%  
%   Copyright 2020 NIST/CLT (steve.blandino@nist.gov)

%% Vargin processing
if isempty(varargin)
    rotateFrame = 0;
else
    rotateFrame = varargin{1};
end
assert(ismember(rotateFrame, [0,1]), 'Wrong Input')

%% Convert from euclidian to quaternions
if rotateFrame
    Q = qtrnConj(qtrnFeul(eucl , 'ZXY'));
else
    Q = qtrnFeul(eucl , 'ZXY');
end
% Qnode    = qtrnFeul(orientation , 'ZXY');




%% Rotate position of node
P = qtrnRotatepoint(Q(2,:),P-Pc)+Pc;

%% Rotate to initial orientation 
% Node rotation
P = qtrnRotatepoint(Q(1,:),P-Pc)+Pc;

%% Output
varargout{1} =  eulFqtrn(qtrnMultiply(Q(2,:), Q(1,:)), 'ZXY');

end