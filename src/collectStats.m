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

addpath('raytracer', 'utils', 'utils/statsUtils')

% specify node ID to select file to be loaded (default: txId=0, rxId=1)
txId = 0;
rxId = 1;

scenario = 'Camillo_NIST_60-GHz Lecture Room--MPCs_1551899818';
% scenario = 'Indoor1';
scenarioPath = fullfile('realScenarios',scenario);
% scenarioPath = fullfile('statsScenarios',scenario);

%% load output file (measurements or ray tracer output)
loadExternalOutput = true;
if loadExternalOutput
    % specify output path 
    dataFolder = "data";
    fileName = 'output.mat';
    currFolder = fullfile(scenarioPath,dataFolder);
    outputPath = fullfile(currFolder,'Output',fileName);
    
    load(outputPath); % load output

else
    %% TX-RX grid
    createNewRxGrid = true;
    if createNewRxGrid
        % create a new TX-RX grid
        % save every execution in a different folder
        currFolder = fullfile(scenarioPath,datestr(now,'mm-dd-yyyy HH-MM-SS'));
        sprintf("Creating new scenario in %s",currFolder);
        
        %% config parameters
        inputPath = fullfile(currFolder,'Input');
        mkdir(inputPath)
        
        cfgParams.numberOfTimeDivisions = 200; % specify the number of different positions for the stats
        genParameterCfg(cfgParams,inputPath);
        
        % create new TX-RX grid
        txPos = [0.2,3,2.5]; % specify the TX position (fixed)
        visualize=true; % plot scatterer positions
        randomSampling = true; % either random or grid positions for collecting the stats
        createTxRxGrid(currFolder,randomSampling,txPos,visualize);
    
    else % setup inputs based on existing files
        currFolder = fullfile(scenarioPath,'10-25-2019 12-12-49');
        inputPath = fullfile(currFolder,'Input');
        
        cfgParams.numberOfTimeDivisions = 108; % specify the number of different positions for the stats
        genParameterCfg(cfgParams,inputPath);
    end
    %% execute the ray-tracer
    outputPath = generateRayTracerOutFiles(currFolder);
    qdPath = fullfile(outputPath,"Ns3","QdFiles",sprintf("Tx%dRx%d.txt",txId,rxId));
    output = readQdFile(qdPath);
    save(fullfile(outputPath,'output.mat'),'output')
end

% check output file fields
if ~isfield(output,"pathGain") || ~isfield(output,"delay") || ~isfield(output,"aodAz") || ~isfield(output,"aodEl")...
        || ~isfield(output,"aodAz") || ~isfield(output,"aoaEl")
    error('Missing field in Tx%sRx%s.txt',txId,rxId);
end

%% extract statistics
statsFilePath = extractStats(output,currFolder);
sprintf("Statistics saved in %s",statsFilePath);

function outputPath=generateRayTracerOutFiles(currPath)
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