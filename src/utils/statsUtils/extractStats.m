function statsFilePath = extractStats(output,currFolder)
%EXTRACTSTATS Function that extracts the statistics of the variables
%contained in output and saves them in statsFilePath. 
%See the documentation for the list of variables and related statistics.
%
% INPUTS:
% - output: struct in the format returned by readQdFile. Stats related to
% the variables contained in output are extracted. Specifically:
%    - pathGain: path gain of each ray [dB]
%    - deltaPathGain: path gain of each ray [dB]
%    - delay: absolute delay of each ray [s]
%    - delaySpread: delay spread of the rays, weighted on the normalized
%    pathGain of each ray [s]
%    - aodEl: elevation of the angle of departure of each ray [deg]
%    - aodAz: azimuth of the angle of departure of each ray [deg]
%    - aoaEl: elevation of the angle of arrival of each ray [deg]
%    - aoaAz: azimuth of the angle of arrival of each ray [deg]
%    - {angle}AritSpread: for all the angles of arrival/departure, the
%    arithmetic angle spread is computed w.r.t. an arithmetic mean of the
%    angles, weighted on the normalized pathGain of each ray [deg]
%    - {angle}PhaseSpread: for all the angles of arrival/departure, the
%    phase angle spread is computed w.r.t. the angular mean of the angles,
%    weighted on the normalized pathGain of each ray. The angular mean is
%    computed summing the complex vectors having phase equal to the 
%    considered angle [deg]
% OUTPUTS:
% - statsFilePath: path to the file containing the extracted statistics. 
%
% SEE ALSO: READQDFILE.

% Copyright (c) 2019, University of Padova, Department of Information
% Engineering, SIGNET lab.
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%    http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

stats = struct(); % to collect the stats

nSteps = size(output,2);
for i = 1:nSteps % row-by-row inspection to extract the stats
    %%% path gain
    pathGain = output(i).pathGain;
    stats(i).pathGain = pathGain;
    
    totPathGain = sum(pathGain);
    normPathGain = pathGain/totPathGain;
    maxPathGain = max(pathGain);
    
    % delta path gain
    deltaPathGain = pathGain-maxPathGain;
    stats(i).deltaPathGain = deltaPathGain;
    
    %%% delay
    delay = output(i).delay;
    stats(i).delay = delay;
    
    % delay spread
    meanDelay = mean(delay);
    delaySpread = sqrt(sum(normPathGain.*(delay-meanDelay).^2)); % weigthed std w.r.t pathGain
    stats(i).delaySpread = delaySpread;
    
    %%% angles
    % DEPARTURE
    % departure elevation angle
    aodElDeg = output(i).aodEl;
    stats(i).aodEl = aodElDeg;
    aodElRad = deg2rad(aodElDeg);
    
    % departure elevation angle arithmetic spread
    meanAodEl = mean(aodElDeg);
    aodElAritSpread = sqrt(sum(normPathGain.*(aodElDeg-meanAodEl).^2));
    stats(i).aodElAritSpread = aodElAritSpread;
    
    % departure elevation angle phase spread
    phaseMeanAodEl = angle(sum(normPathGain.*exp(1i*aodElRad))); % radians
    aodElPhaseSpread = sqrt(sum(normPathGain.*angle(exp(1i*aodElRad)-exp(1i*phaseMeanAodEl)).^2)); % radians
    stats(i).aodElPhaseSpread = rad2deg(aodElPhaseSpread); % convert to degrees
    
    % departure azimuth angle
    aodAzDeg = output(i).aodAz;
    stats(i).aodAz = aodAzDeg;
    aodAzRad = deg2rad(aodAzDeg);
    
    % departure azimuth angle arithmetic spread
    meanAodAz = mean(aodAzDeg);
    aodAzAritSpread = sqrt(sum(normPathGain.*(aodAzDeg-meanAodAz).^2));
    stats(i).aodAzAritSpread = aodAzAritSpread;
    
    % departure azimuth angle phase spread
    phaseMeanAodAz = angle(sum(normPathGain.*exp(1i*aodAzRad))); % radians
    aodAzPhaseSpread = sqrt(sum(normPathGain.*angle(exp(1i*aodAzRad)-exp(1i*phaseMeanAodAz)).^2)); % radians
    stats(i).aodAzPhaseSpread = rad2deg(aodAzPhaseSpread); % convert to degrees
    
    % ARRIVAL
    % arrival elevation angle
    aoaElDeg = output(i).aoaEl; % degrees
    stats(i).aoaElDeg = aoaElDeg;
    aoaElRad = deg2rad(aoaElDeg);
    
    % arrival elevation angle arithmetic spread
    meanAoaEl = mean(aoaElDeg);
    aoaElAritSpread = sqrt(sum(normPathGain.*(aoaElDeg-meanAoaEl).^2));
    stats(i).aoaElAritSpread = aoaElAritSpread;
    
    % arrival elevation angle phase spread
    phaseMeanAoaEl = angle(sum(normPathGain.*exp(1i*aoaElRad))); % radians
    aoaElPhaseSpread = sqrt(sum(normPathGain.*angle(exp(1i*aoaElRad)-exp(1i*phaseMeanAoaEl)).^2)); % radians
    stats(i).aoaElPhaseSpreadDeg = rad2deg(aoaElPhaseSpread); % convert to degrees
    
    % arrival azimuth angle
    aoaAzDeg = output(i).aoaAz; % degrees
    stats(i).aoaAz = aoaAzDeg;
    aoaAzRad = deg2rad(aoaAzDeg);
    
    % arrival azimuth angle arithmetic spread
    meanAoaAz = mean(aoaAzDeg);
    aoaAzAritSpread = sqrt(sum(normPathGain.*(aoaAzDeg-meanAoaAz).^2));
    stats(i).aoaAzAritSpread = aoaAzAritSpread;
    
    % arrival azimuth angle phase spread
    phaseMeanAoaAz = angle(sum(normPathGain.*exp(1i*aoaAzRad))); % radians
    aoaAzPhaseSpread = sqrt(sum(normPathGain.*angle(exp(1i*aoaAzRad)-exp(1i*phaseMeanAoaAz)).^2)); % radians
    stats(i).aoaAzPhaseSpread = rad2deg(aoaAzPhaseSpread); % convert to degrees
end
statsFolder = fullfile(currFolder,'Stats');
mkdir(statsFolder)
statsFilePath = fullfile(statsFolder,'stats.mat');
save(statsFilePath,'stats')
end