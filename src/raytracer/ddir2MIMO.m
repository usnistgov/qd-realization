function [ch] = ddir2MIMO(ddir, info, ptr)
%%DDIR2MIMO Converts the double direction impulse response in the MIMO 
% channel matrix assigning phase rotations according with PAA centroids 
% positions and angles of departure/arrival
%
%Copyright 2019-2020 NIST/CLT (steve.blandino@nist.gov)


[~,b] = unique(info{ptr.nt}.centroids);
b = sort(b);
idx = find(info{ptr.nt}.centroids == info{ptr.nt}.centroids(b(ptr.paatx)));
t= idx(ptr.rot_tx);

[~,b] = unique(info{ptr.nr}.centroids);
b = sort(b);
idx = find(info{ptr.nr}.centroids == info{ptr.nr}.centroids(b(ptr.paarx)));
r= idx(ptr.rot_rx);

AOD_az = ddir(:,10)/180*pi;
AOD_el = ddir(:,11)/180*pi;
AOA_az = ddir(:,12)/180*pi;
AOA_el = ddir(:,13)/180*pi;

R_AOD = phaseRotation(AOD_az,AOD_el, info{ptr.nt}.centroids_shift{t} );
R_AOA = phaseRotation(AOA_az,AOA_el, info{ptr.nr}.centroids_shift{r} );
N_AOD = size(R_AOD,2);
N_AOA = size(R_AOA,2);

ch = zeros([size(ddir), N_AOD*N_AOA]);

for id_aod = 1:N_AOD
    for id_aoa = 1:N_AOA
        ch(:,:, (id_aod-1)*(N_AOA)+id_aoa) = ddir;
        ch(:,18, (id_aod-1)*(N_AOA)+id_aoa) = wrapTo2Pi(ch(:,18,id_aod)+angle(R_AOD(:,id_aod))+angle(R_AOA(:,id_aoa)));
    end
end

end