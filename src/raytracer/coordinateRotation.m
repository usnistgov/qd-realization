function [P, varargout] = coordinateRotation(P, C, euler, varargin)
%COORDINATEROTATION rotates a point in the QD software using quaternions.
%
%   P = COORDINATEROTATION(P,PC,EUL) or 
%   P = COORDINATEROTATION(P,PC,EUL, 'point') rotates the point P = (px, py, pz)
%   with respect to a reference system centered in C = (cx,cy,cz) by the
%   euclidians angles EUL. EUL is the Nx3 matrix where N indicates
%   consecutive rotations. If C is a 1x3 vector the same center is assumed
%   for each of the N rotations. If C is a Nx3 matrix a each rotation uses
%   a different center of rotation.
%   The point rotation is used to model the device rotation, which brings
%   the PAAs  in  a  different  position  in  the global  frame.
%   Default euclidian sequence is assumed 'ZXY'.
%
%   P = COORDINATEROTATION(P,PC,EUL, 'frame') returns the coordinates of
%   the point P when the reference frame centered in PC is rotated by the
%   euclidians angles EUL. The  frame  rotation  is used to transform AOA
%   and AOD from global to local coordinates. Default euclidian sequence is
%   assumed 'ZXY'.
%
%  [P SUCCESSIVE_EUL] = COORDINATEROTATION(..) returns the euler angle
%  SUCCESSIVE_EUL equivalent to the consecutive rotations in EUL
%
%   Copyright 2020 NIST/CTL (steve.blandino@nist.gov)

%% Vargin processing
if isempty(varargin)
    rotateFrame = 'point';
else
    rotateFrame = varargin{1};
end

assert(ismember(rotateFrame, {'point','frame'}), 'Wrong Input')
assert(size(C,1) == size(euler,1) || size(C,1) ==1, 'Provide correct', ...
'centroids to perform rotation')

Nrotations = size(euler,1) ;

if isempty(P)
    return
end

if size(C,1) ==1
    C = repmat(C,[Nrotations,1]);
end

%% Convert from Euler to quaternions
switch rotateFrame
    case 'frame'
        Q = qtrnConj(qtrnFeul(euler , 'ZXY'));
    case 'point'
        Q = qtrnFeul(euler , 'ZXY');
end

%% Apply rotations and compute equivalent quaternion
for j = 1:Nrotations    
    P = qtrnRotatepoint(Q(j,:),P-C(j,:))+C(j,:);    
    if j <= Nrotations - 1
        Q(j,:) = qtrnMultiply(Q(j,:), Q(j+1,:)); % Store at index j the prod 1:j+1
    end
end

%% Output
varargout{1} = eulFqtrn(Q(max(Nrotations-1,1),:), 'ZXY'); % Return the euclidian angle
% correspondent to the quaternion Q at index Nrotations-1
end
