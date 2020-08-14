function [dod, doa, AOD_az, AOD_el, AOA_az, AOA_el] = frameRotation(frmRotMpInfo, orientation)
    
% dod=coordinateRotation(frmRotMpInfo.dod,[0 0 0], orientation.tx, 'frame');
dod = reshape(cell2mat(arrayfun(@(x) coordinateRotation(x.dod, [0 0 0], orientation.tx, 'frame'), frmRotMpInfo, 'UniformOutput', false)), 3, []);
% delay is the total length of multipath
% delay=norm(dod);

% doa= coordinateRotation(frmRotMpInfo.doa, [0 0 0], orientation.rx,'frame');
doa = reshape(cell2mat(arrayfun(@(x) coordinateRotation(x.doa, [0 0 0], orientation.rx, 'frame'), frmRotMpInfo, 'UniformOutput', false)), 3, []);

% Aod azimuth
AOD_az = mod(atan2d(dod(2,:),dod(1,:)), 360).'/180*pi;
% Aod elevation
AOD_el = acosd(dod(3,:)./vecnorm(dod)).'/180*pi;
% Aoa azimuth
AOA_az = mod(atan2d(doa(2,:),doa(1,:)), 360).'/180*pi;
% Aoa elevation
AOA_el = acosd(doa(3,:)./vecnorm(doa)).'/180*pi;
dod = dod.';
doa = doa.';

end