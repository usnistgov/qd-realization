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

clear
close all
clc

addpath('raytracer', 'utils')

% specify node ID to select file to be loaded (default: txId=0, rxId=1)
txId = 0;
rxId = 1;

scenario = 'L-Room';
scenarioPath = fullfile('statsScenarios',scenario);

load = true;
if load
    loadFolder = "10-18-2019 14-06-25";
    currFolder = fullfile(scenarioPath,loadFolder);
    sprintf("Loading files from %s.",currFolder);
else
    % save every execution in a different folder
    currFolder = fullfile(scenarioPath,datestr(now,'mm-dd-yyyy HH-MM-SS'));
    sprintf("Creating new scenario in %s",currFolder);
    
    %% config parameters
    inputPath = fullfile(currFolder,'Input');
    mkdir(inputPath)
    cfgParams.numberOfTimeDivisions = 200; % specify the number of different positions for the stats
    genParameterCfg(cfgParams,inputPath);

    %% create TX-RX grid
    % specify the TX position
    txPos = [0.2,3,2.5]; % the transmitter position is fixed
    visualize=true; % plot scatterer positions
    randomSampling = true; % either random or grid positions for collecting the stats
    createTxRxGrid(currFolder,randomSampling,txPos,visualize);

    %% execute the ray-tracer
    generateRayTracerOutFiles(currFolder)
end
%% extract statistics
outputPath = fullfile(currFolder,'Output','Ns3','QdFiles',sprintf('Tx%dRx%d.txt',txId,rxId));
output = readQdFile(outputPath);

if ~isfield(output,"pathGain") || ~isfield(output,"delay") || ~isfield(output,"aodAz") || ~isfield(output,"aodEl")...
        || ~isfield(output,"aodAz") || ~isfield(output,"aoaEl")
    error('Missing field in Tx%sRx%s.txt',txId,rxId);
end

stats = struct(); % to collect the stats

nSteps = length(output);
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

    %%% departure angles
    % departure elevation angle
    aodEl = output(i).aodEl;
    stats(i).aodEl = aodEl;
    
    % departure elevation angle arithmetic spread
    meanAodEl = mean(aodEl);
    aodElAritSpread = sqrt(sum(normPathGain.*(aodEl-meanAodEl).^2));
    stats(i).aodElAritSpread = aodElAritSpread;
    
    % departure elevation angle phase spread
    phaseMeanAodEl = angle(sum(exp(1i*normPathGain.*aodEl)));
    aodElPhaseSpread = sqrt(sum(normPathGain.*(aodEl-phaseMeanAodEl).^2));
    stats(i).aodElPhaseSpread = aodElPhaseSpread;
    
    % departure azimuth angle
    aodAz = output(i).aodAz;
    stats(i).aodAz = aodAz;
    
    % departure azimuth angle arithmetic spread
    meanAodAz = mean(aodAz);
    aodAzAritSpread = sqrt(sum(normPathGain.*(aodAz-meanAodAz).^2));
    stats(i).aodAzAritSpread = aodAzAritSpread;
    
    % arrival elevation angle
    aoaEl = output(i).aoaEl;
    stats(i).aoaEl = aoaEl;
    
    % arrival elevation angle arithmetic spread
    meanAoaEl = mean(aoaEl);
    aoaElAritSpread = sqrt(sum(normPathGain.*(aoaEl-meanAoaEl).^2));
    stats(i).aoaElAritSpread = aoaElAritSpread;
    
    % arrival azimuth angle
    aoaAz = output(i).aoaAz;
    stats(i).aoaAz = aoaAz;
    
    % arrival azimuth angle arithmetic spread
    meanAoaAz = mean(aoaAz);
    aoaAzAritSpread = sqrt(sum(normPathGain.*(aoaAz-meanAoaAz).^2));
    stats(i).aoaAzAritSpread = aoaAzAritSpread;
end
%% plot statistics
features = ["pathGain"]; % specify the metric to observe for the statistics

for feat = features
    if ~isfield(output,feat)
        error('feature must be a field of Tx%sRx%s.txt',txId,rxId);
    end
    currFeature = stats.(feat);
    
    figure()
    subplot(1,2,1)
    ecdf(currFeature)
    xlabel('$x=$ Path Gain')
    ylabel('$F(x)$')
    % title('Path Gain distribution in a rectangular room. $N_{refl}=2$')

    subplot(1,2,2)
    histogram(currFeature)
    xlabel('Path Gain')
    ylabel('Count')
    % title('Path Gain distribution in a rectangular room. $N_{refl}=2$')
end

function []=generateRayTracerOutFiles(currPath)
    %% Initialization 
    rootFolderPath = pwd;
    fprintf('-------- NIST/CTL QD mmWave Channel Model --------\n');
    fprintf('Current root folder:\n\t%s\n',rootFolderPath);
    [path,folderName] = fileparts(rootFolderPath);
    if strcmp(folderName, 'src')
        fprintf('Start to run.\n');
    else
        error('The root folder should be ''src''');
    end

    %% Input
    endout=regexp(currPath,filesep,'split');
    scenarioNameStr = endout{2};
    sprintf('Use scenario: %s.\n',scenarioNameStr);

    % Check Input Scenario File
    scenarioInputPath = fullfile(currPath,'/Input');

    % Input System Parameters
    paraCfg = parameterCfg(currPath);
    % Input Node related parameters
    [paraCfg, nodeCfg] = nodeProfileCfg(paraCfg);

    % To collect the statistics, we need these two switches
    paraCfg.switchSaveVisualizerFiles = 1;
    paraCfg.mobilitySwitch = 1;
    % Run raytracing function and generate outputs
    outputPath = Raytracer_collectStats(paraCfg, nodeCfg);

    fprintf('Save output data to:\n%s\n',outputPath);
    fprintf('--------- Simulation Complete ----------\n');
end