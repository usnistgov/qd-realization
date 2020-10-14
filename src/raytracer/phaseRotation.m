function R = phaseRotation(theta, phi, centrShift, varargin)
%%PHASEROTATION returns the channel phase rotation R with respect the
%%centroid position in a 60GHz channel.
%        R = phaseRotation(theta,phi, centr_shift)
%        
%        **theta, phi are elevation and azimut angles in rad of rays
%        impinging on the array
%        **centr_shift is the shift wrt the centroid position in which R
%        needs to be computed
%
%        R = PHASEROTATION(theta,phi, centr_shift, 'fc', [value])
%        Computes R when the central frequency is set to [value] Hz
%
%Copyright 2019-2020 NIST/CTL (steve.blandino@nist.gov)

%% Varargin processing
var_struct = {'fc'};
for k = 1:2:length(varargin)
    if (~any(strcmp(varargin{k}, var_struct)))
        warning('Cannot specify "%s" as input value - it will be discarded', varargin{k});
    end
    eval([varargin{k},' = varargin{k+1};'])
end
clear('k')
if ~exist('frequency','var'),     fc=60e9;        end

dx = centr_shift(:,1);
dy = centr_shift(:,2);
dz = centr_shift(:,3);

%% 
lambda = 3e8/fc;
k = 2*pi/lambda; 

kx =  k*(sin(theta).*cos(phi))*dx.';
ky =  k*sin(theta).*sin(phi)*dy.';
kz =  k*cos(theta)*dz.';

R = (exp(1i*kz).*exp(1i*ky).*exp(1i*kx)); %6.87A Balanis 4ed
end
