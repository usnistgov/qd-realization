function output = qdGenerator(dRayOutput, arrayOfMaterials,...
    MaterialLibrary, freq, vRel)
% output:
% 1. reflection order
% 2:4. DoD
% 5:7. DoA
% 8. delay
% 9. Path Gain
% 10. AoD Azimuth
% 11. AoD Elevation
% 12. AoA Azimuth
% 13. AoA Elevation
%     14:15. PolarizationTx(1,:)
%     16:17. PolarizationTx(2,:)
% 18. phase (reflOrder*pi)
%     19. Cross-pol path gain
% 20. phase (dopplerFactor*freq)
% 21. 0 (?)


warning('TODO: Add randomness to reflection loss of deterministic ray')

% Pre/post cursors output
outputPre = getQdOutput(dRayOutput, arrayOfMaterials, MaterialLibrary,...
    freq, vRel, 'pre');
outputPost = getQdOutput(dRayOutput, arrayOfMaterials, MaterialLibrary,...
    freq, vRel, 'post');

output = [outputPre; dRayOutput; outputPost];

end


%% Utils
function output = getQdOutput(dRayOutput, arrayOfMaterials, MaterialLibrary,...
    freq, vRel, prePostParam)
params = getParams(arrayOfMaterials, MaterialLibrary, prePostParam);

% delays
tau0 = dRayOutput(8); % main cursor's delay

interArrivalTime = rndExp(params.mul, params.nRays, 1); % [s]
taus = tau0 + params.delayMultiplier*cumsum(interArrivalTime);
% TODO: check if ray is arriving before LoS

% path gains
exponent = -abs(taus - tau0)/params.muy;
s = randn(params.nRays, 1)*params.mus; % TODO: check pre/post, check mus dB/lin
pgCursor = 10^(dRayOutput(9) / 10); % path gain of main cursor TODO: lin? db?
pg = pgCursor / params.muk * exp(exponent + s); % TODO: something wrong with units
pg = 10*log10(pg);
% in dB: pg = pgCursorDb - params.mukDb + 10*log10(exp(1)) * (exponent+s);
% PROBLEM: to avoid diffused components to have more power than main
% cursor, s < params.mukDb / (10*log10(exp(1))) - exponent
% TODO: should the cursor's power be maintained while dividing it over the
% diffused components?

% angle spread
aodAzCursor = dRayOutput(10);
aodElCursor = dRayOutput(11);
aoaAzCursor = dRayOutput(12);
aoaElCursor = dRayOutput(13);

[aodAz, aodEl] = getDiffusedAngles(aodAzCursor, aodElCursor,...
    params.mus, params.nRays);
[aoaAz, aoaEl] = getDiffusedAngles(aoaAzCursor, aoaElCursor,...
    params.mus, params.nRays);

% doppler spread
aoaDirections = azel2vec(aoaAz, aoaEl);
vProjected = aoaDirections * vRel.'; % assuming vRel row vector. vRel should be the unwrapped relative velocity between the two nodes wrt to the Rx
dopplerDeltaFreq = vProjected / 3e8 * freq;

% Combine results into output matrix
% Copy D-ray outputs as some columns are repeated (e.g., reflection order)
output = repmat(dRayOutput, params.nRays, 1);
% delay
output(:,8) = taus;
% path gain
output(:,9) = pg;
% AoD azimuth
output(:,10) = aodAz;
% AoD elevation
output(:,11) = aodEl;
% AoA azimuth
output(:,12) = aoaAz;
% AoA elevation
output(:,13) = aoaEl;
% phase
output(:,18) = rand(params.nRays,1) * 2*pi;
% doppler frequency shift
output(:,20) = dopplerDeltaFreq;

end


function params = getParams(arrayOfMaterials, MaterialLibrary, prePostParam)

materialIdx = arrayOfMaterials(end); % QD based on last reflector

switch(prePostParam)
    case 'pre'
        params.muk = 10^(MaterialLibrary.mu_k_Precursor(materialIdx) / 10); % lin
        params.muy = MaterialLibrary.mu_Y_Precursor(materialIdx) * 1e-9; % s
        params.mul = MaterialLibrary.mu_lambda_Precursor(materialIdx) * 1e9; % 1/s
        params.delayMultiplier = -1;
        params.nRays = 3;
        
    case 'post'
        params.muk = 10^(MaterialLibrary.mu_k_Postcursor(materialIdx) / 10); % lin
        params.muy = MaterialLibrary.mu_Y_Postcursor(materialIdx) * 1e-9; % s
        params.mul = MaterialLibrary.mu_lambda_Postcursor(materialIdx) * 1e9; % 1/s
        params.delayMultiplier = 1;
        params.nRays = 16;
        
    otherwise
        error('prePostParam=''%s''. Should be ''pre'' or ''post''', prePostParam)
end

params.mus = MaterialLibrary.mu_sigmaTheta(materialIdx); % TODO: db? lin?

end


function [aodAz, aodEl] = getDiffusedAngles(azCursor, elCursor, angleSpread, nRays)
aodAz = azCursor + rndLaplace(angleSpread, nRays, 1);
aodEl = elCursor + rndLaplace(angleSpread, nRays, 1);
[aodAz, aodEl] = wrapAngles(aodAz, aodEl);

end


function [az, el] = wrapAngles(az, el)
% If elevation is negative, bring it back in [0,180] and rotate azimuth by
% half a turn
negativeElMask = el < 0;
el(negativeElMask) = -el(negativeElMask);
az(negativeElMask) = az(negativeElMask) + 180;

% If elevation is over 180, bring it back in [0,180] and rotate azimuth by
% half a turn
over180ElMask = el > 180;
el(over180ElMask) = 360 - el(over180ElMask);
az(over180ElMask) = az(over180ElMask) + 180;

% Wrap azimuth to [0,360)
az = mod(az, 360);

end