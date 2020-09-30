function [ch, varargout] = ddir2MIMO(ddir, info, frmRotMpInfo, ptr)
%%DDIR2MIMO Converts the double direction impulse response in the MIMO 
% channel matrix assigning phase rotations according with PAA centroids 
% positions and angles of departure/arrival
%
%Copyright 2019-2020 NIST/CLT (steve.blandino@nist.gov)


[~,b] = unique(info{ptr.nt}.centroids);
b = sort(b);
idx = find(info{ptr.nt}.centroids == info{ptr.nt}.centroids(b(ptr.paatx)));
t= idx(ptr.iid_tx); %Pointer to centroid in cell array in info{ptr.nt}

[~,b] = unique(info{ptr.nr}.centroids);
b = sort(b);
idx = find(info{ptr.nr}.centroids == info{ptr.nr}.centroids(b(ptr.paarx)));
r= idx(ptr.iid_rx); %Pointer to centroid in cell array in info{ptr.nr}

%% Overwrite angles

% N_AOD = size(info{ptr.nt}.centroidsShift{t},1);
N_AOD = info{ptr.nt}.nodePAAInfo{ptr.paatx}.rotated_channel(ptr.iid_tx);
% N_AOA = size(info{ptr.nr}.centroidsShift{r},1);
N_AOA = info{ptr.nr}.nodePAAInfo{ptr.paarx}.rotated_channel(ptr.iid_rx);

% orientation.tx = info{ptr.nt}.orientation(b(ptr.paatx),:);
% orientation.rx = info{ptr.nr}.orientation(b(ptr.paarx),:);
% [AOD_az,AOD_el,AOA_az, AOA_el] = frameRotation(frmRotMpInfo, orientation);
% 
% ddir(:,10) = AOD_az*180/pi;
% ddir(:,11) = AOD_el*180/pi;
% ddir(:,12) = AOA_az*180/pi;
% ddir(:,13) = AOA_el*180/pi;
% 
% R_AOD = phaseRotation(AOD_az,AOD_el, info{ptr.nt}.centroidsShift{t} );
% R_AOA = phaseRotation(AOA_az,AOA_el, info{ptr.nr}.centroidsShift{r} );


ch = zeros([size(ddir), N_AOD*N_AOA]);

for id_aod = 1:N_AOD
    for id_aoa = 1:N_AOA
        
        orientation.tx = info{ptr.nt}.orientation{t}(id_aod,:);
        orientation.rx = info{ptr.nr}.orientation{r}(id_aoa,:);
        
        [dod, doa, AOD_az,AOD_el,AOA_az, AOA_el] = frameRotation(frmRotMpInfo, orientation);
        % dod - direction of departsure
        ddir(:, 2:4) = dod;
        % doa - direction of arrival
        ddir(:, 5:7) = doa;
        ddir(:,10) = AOD_az*180/pi;
        ddir(:,11) = AOD_el*180/pi;
        ddir(:,12) = AOA_az*180/pi;
        ddir(:,13) = AOA_el*180/pi;
        R_AOD = phaseRotation(AOD_az,AOD_el, info{ptr.nt}.centroidsShift{t}(id_aod,:) );
        R_AOA = phaseRotation(AOA_az,AOA_el, info{ptr.nr}.centroidsShift{r}(id_aoa,:) );
        
        ch(:,:, (id_aod-1)*(N_AOA)+id_aoa) = ddir;
        ch(:,18, (id_aod-1)*(N_AOA)+id_aoa) = wrapTo2Pi(ch(:,18,id_aod)+angle(R_AOD)+angle(R_AOA));
    end
end
[A,B] = meshgrid(info{ptr.nt}.paaInCluster{t}, info{ptr.nr}.paaInCluster{r});
varargout{1} = [A(:), B(:)];
end