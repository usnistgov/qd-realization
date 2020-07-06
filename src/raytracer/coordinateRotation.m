function [P, varargout] = coordinateRotation(P, C, eucl, varargin)
%COORDINATEROTATION rotates a point in the QD software using quaternions.
%
%   P = COORDINATEROTATION(P,PC,EUCL) rotates the point P = (px, py, pz) 
%   with respect to a reference system centered in C = (cx,cy,cz) by the 
%   euclidians angles EUCL. EUCL is the Nx3 matrix where N indicates 
%   consecutive rotations.
%   The point rotation is used to model the device rotation, which brings 
%   the PAAs  in  a  different  position  in  the global  frame.  
%   Default euclidian sequence is assumed 'ZXY'.
%  
%  
%   P = COORDINATEROTATION(P,PC,EUCL, 'frame') returns the coordinates of 
%   the point P when the reference frame centered in PC is rotated by the
%   euclidians angles EUCL. The  frame  rotation  is used to transform AOA 
%   and AOD from global to local coordinates. Default euclidian sequence is 
%   assumed 'ZXY'.
% 
%  [P Q] = POINTROTATION(..) returns the quaternion equivalent to the
%  consecutive rotations in EUCL
%  
%   Copyright 2020 NIST/CLT (steve.blandino@nist.gov)

%% Vargin processing
if isempty(varargin)
    rotateFrame = 'point';
else
    rotateFrame = varargin{1};
end
assert(ismember(rotateFrame, {'point','frame'}), 'Wrong Input')

%% Convert from euclidian to quaternions
switch rotateFrame
    case 'frame'
    Q = qtrnConj(qtrnFeul(eucl , 'ZXY'));
    case 'point'
    Q = qtrnFeul(eucl , 'ZXY');
end

%% Rotate position of node
P = qtrnRotatepoint(Q(2,:),P-C)+C;

%% Rotate to initial orientation 
% Node rotation
P = qtrnRotatepoint(Q(1,:),P-C)+C;

%% Output
varargout{1} =  eulFqtrn(qtrnMultiply(Q(2,:), Q(1,:)), 'ZXY');

end