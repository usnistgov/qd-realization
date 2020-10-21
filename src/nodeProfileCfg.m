function [paraCfg, nodeCfg] = nodeProfileCfg(paraCfg)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


% -------------Software Disclaimer---------------
%
% NIST-developed software is provided by NIST as a public service. You may use, copy
% and distribute copies of the software in any medium, provided that you keep intact this
% entire notice. You may improve, modify and create derivative works of the software or
% any portion of the software, and you may copy and distribute such modifications or
% works. Modified works should carry a notice stating that you changed the software
% and should note the date and nature of any such change. Please explicitly
% acknowledge the National Institute of Standards and Technology as the source of the
% software.
%
% NIST-developed software is expressly provided "AS IS." NIST MAKES NO
% WARRANTY OF ANY KIND, EXPRESS, IMPLIED, IN FACT OR ARISING BY
% OPERATION OF LAW, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
% WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE,
% NON-INFRINGEMENT AND DATA ACCURACY. NIST NEITHER REPRESENTS
% NOR WARRANTS THAT THE OPERATION OF THE SOFTWARE WILL BE
% UNINTERRUPTED OR ERROR-FREE, OR THAT ANY DEFECTS WILL BE
% CORRECTED. NIST DOES NOT WARRANT OR MAKE ANY REPRESENTATIONS
% REGARDING THE USE OF THE SOFTWARE OR THE RESULTS THEREOF,
% INCLUDING BUT NOT LIMITED TO THE CORRECTNESS, ACCURACY,
% RELIABILITY, OR USEFULNESS OF THE SOFTWARE.
%
% You are solely responsible for determining the appropriateness of using
% and distributing the software and you assume all risks associated with its use, including
% but not limited to the risks and costs of program errors, compliance with applicable
% laws, damage to or loss of data, programs or equipment, and the unavailability or
% interruption of operation. This software is not intended to be used in any situation
% where a failure could cause risk of injury or damage to property. The software
% developed by NIST employees is not subject to copyright protection within the United
% States.
%
% Modified by: Mattia Lecci <leccimat@dei.unipd.it>, Updated implementation

warning('off', 'MATLAB:MKDIR:DirectoryExists');

scenarioNameStr = paraCfg.inputScenarioName;
% Input Parameters to be Updated
mobilitySwitch = paraCfg.mobilitySwitch;
mobilityType = paraCfg.mobilityType;
numberOfNodes = paraCfg.numberOfNodes;
numberOfTimeDivisions = paraCfg.numberOfTimeDivisions;
switchRandomization = paraCfg.switchRandomization;
numberTracePoints = paraCfg.numberOfTimeDivisions;

% List of paths
inputPath = fullfile(scenarioNameStr, 'Input');

%% Code
nodeRotationTime= zeros(paraCfg.numberOfTimeDivisions,3, paraCfg.numberOfNodes);
nodeInitialPosition = zeros( paraCfg.numberOfNodes, 3);
nodePositionTime = zeros(paraCfg.numberOfTimeDivisions,3, paraCfg.numberOfNodes);

% Try to open previous config: nodes.dat. If it exist convert it into new
% config i.e. NodePositionX.dat
try
    obsoletePosition = csvread(fullfile(inputPath, 'nodes.dat'));
    warning('Configuration obsolete. nodes.dat not used anymore: node information loaded from NodePosition.dat')
    %     obsoleteConfig = true;
    listing = dir(fullfile(scenarioNameStr, 'Input'));
    for nodeId = 1: size(obsoletePosition,1)
        if ~sum(arrayfun(@(x) strcmp(x.name,['NodePosition',num2str(nodeId-1), '.dat']), listing))
            writematrix(obsoletePosition(nodeId,:), fullfile(inputPath, ['NodePosition', num2str(nodeId-1),'.dat']))
        end
    end
catch
    %     obsoleteConfig = false;
end

%% Random generation of node positions
if switchRandomization == 1
    xCoordinateRandomizer = rand * 8 + 1;
    yCoordinateRandomizer = rand * 17 + 1;
    zCoordinateRandomizer = 2.5;
    Tx = [xCoordinateRandomizer, yCoordinateRandomizer, zCoordinateRandomizer];
    
    xCoordinateRandomizer = rand * 8 + 1;
    yCoordinateRandomizer = rand * 17 + 1;
    zCoordinateRandomizer = 1.6;
    Rx = [xCoordinateRandomizer, yCoordinateRandomizer, zCoordinateRandomizer];
    
    if mobilityType ~= 1
        mobilityType = 1;
        warning('Changing mobilityType to %d', mobilityType)
    end
end

%% Extracting data from nodes.dat and nodeVelocities.dat file.
% nodeVelocities contains node velocities
if switchRandomization == 0
    if mobilitySwitch && mobilityType ==1
        try
            nodeVelocities = csvread(fullfile(inputPath, 'nodeVelocities.dat'));
        catch
            switchRandomization = 1;
            warning(['Unable to read nodeVelocities.dat. ',...
                'Changing switchRandomization to %d'], switchRandomization)
        end
    else
        nodeVelocities = zeros(numberOfNodes, 3);
        try
            csvread(fullfile(inputPath, 'nodeVelocities.dat'));
            warning('nodeVelocities.dat not used')
        catch
        end
    end
    
    %% Load NodePositionX.dat and NodeRotationX.dat
    listing = dir(fullfile(scenarioNameStr, 'Input'));
    countListing = 0;
    nodePositionTimeRaw = cell(numberOfNodes,1);
    nodeRotationTimeRaw = cell(numberOfNodes,1);
    
    for iterateNumberOfNodes = 1:numberOfNodes        
            nodePositionFile = sprintf('NodePosition%d.dat',iterateNumberOfNodes-1);
            nodeRotationFile = sprintf('NodeRotation%d.dat',iterateNumberOfNodes-1);
            isNodePosition = any(arrayfun(@(x) strcmp(x.name,nodePositionFile), listing));
            isNodeRotation = any(arrayfun(@(x) strcmp(x.name,nodeRotationFile), listing));

            if ~isNodePosition 
                error([ nodePositionFile, ' not defined.']);
            end
            nodePositionTimeRaw{iterateNumberOfNodes} = readmatrix(fullfile(inputPath,nodePositionFile));

            if ~isNodeRotation
                nlines = size(nodePositionTimeRaw{iterateNumberOfNodes},1);
                nodeRotationTimeRaw{iterateNumberOfNodes} = repmat([0 0  0], nlines,1);
                writematrix(nodeRotationTimeRaw{iterateNumberOfNodes}, fullfile(inputPath,nodeRotationFile));
                warning([nodeRotationFile, ' not defined. Rotation set to [0,0,0] for all time instances.'])
            else
                nodeRotationTimeRaw{iterateNumberOfNodes} = readmatrix(fullfile(inputPath,nodeRotationFile));
            end
    end
    
    %%  Load NodePositionX.dat and NodeRotationX.dat
    for iterateNumberOfNodes = 1:numberOfNodes
        
        % NodePosition processing
        nodePositionTimeTmp = nodePositionTimeRaw{iterateNumberOfNodes};
        if mobilityType == 1 && mobilitySwitch
            nodeInitialPosition(iterateNumberOfNodes,:) = nodePositionTimeTmp;
            countListing = countListing + 1;
            
        else
            
            timeSamplesFile = size(nodePositionTimeTmp,1);
            if  timeSamplesFile< paraCfg.numberOfTimeDivisions &&  ...
                    timeSamplesFile > 1
                nodePositionTime = nodePositionTime(1:timeSamplesFile, :,:);
                paraCfg.numberOfTimeDivisions = timeSamplesFile;
                numberOfTimeDivisions = paraCfg.numberOfTimeDivisions;
                warning('Time divisition too long.')
            end
            
            try
                numberTracePoints = min(timeSamplesFile,paraCfg.numberOfTimeDivisions);
                nodePositionTime(1:numberTracePoints, :, iterateNumberOfNodes) = nodePositionTimeTmp(1:numberTracePoints,:);
                nodePositionTime(numberTracePoints+1:end, :, iterateNumberOfNodes) = repmat(nodePositionTimeTmp, [paraCfg.numberOfTimeDivisions-numberTracePoints,1,1]);
                nodeInitialPosition(iterateNumberOfNodes,:) = squeeze(nodePositionTime(1,:,iterateNumberOfNodes));
                countListing = countListing + 1;
            catch
                mobilityType = 1;
                warning('Node Position input incorrect. Changing mobilityType to 1');
            end
        end
        
        % NodeRotation processing
        nodeRotationTimeTemp = nodeRotationTimeRaw{iterateNumberOfNodes};
        if mobilityType == 1 && mobilitySwitch
            nodeOrientation(iterateNumberOfNodes,:) = squeeze(nodeRotationTime(1,:,iterateNumberOfNodes));
            countListing = countListing + 1;
        else
            
            if  size(nodeRotationTimeTemp,1)< size(nodeRotationTime,1) &&  ...
                    size(nodeRotationTimeTemp,1) > 1
                nodeRotationTime = nodeRotationTime(1:size(nodeRotationTimeTemp,1), :,:);
                paraCfg.numberOfTimeDivisions = size(nodeRotationTimeTemp,1) ;
                numberOfTimeDivisions = paraCfg.numberOfTimeDivisions;
                warning('Time divisition too long.')
            end
            
            try
                numberTracePoints =  min(size(nodeRotationTimeTemp,1),paraCfg.numberOfTimeDivisions);
                nodeRotationTime(1:numberTracePoints, :, iterateNumberOfNodes) = nodeRotationTimeTemp(1:numberTracePoints,:);
                nodeRotationTime(numberTracePoints+1:end, :, iterateNumberOfNodes) = repmat(nodeRotationTimeTemp, [paraCfg.numberOfTimeDivisions-numberTracePoints,1,1]);
                nodeOrientation(iterateNumberOfNodes,:) = squeeze(nodeRotationTime(1,:,iterateNumberOfNodes));
                countListing = countListing + 1;
            catch
                error('Node Rotation config incorrect');
            end
        end
    end
    %%
    
    if mobilityType == 2 && countListing < numberOfNodes
        warning(['Node Position input incorrect. Linear mobility',...
            'model is chosen']);
        mobilityType = 1;
    elseif mobilityType == 2 && countListing == numberOfNodes
        % Cannot compute last velocity, so stop one iteration earlier
        if numberOfTimeDivisions ~= size(nodePositionTime, 1) - 1
            numberOfTimeDivisions = size(nodePositionTime, 1) - 1;
            warning('Changing numberOfTimeDivisions to %d', numberOfTimeDivisions)
        end
    end
    
    if size(nodeInitialPosition,1) ~= size(nodeVelocities,1) && mobilitySwitch == 1 && mobilityType==1
        error(['nodes.dat and nodeVelocities.dat do not have same number ',...
            'of rows. Please check the input files in the Input folder.'])
    end
    
    if ~isempty(numberOfNodes) && numberOfNodes ~= size(nodeInitialPosition, 1)
        warning(['"numberOfNodes" parameter does not match the number of ',...
            'nodes given in file. The "numberOfNodes" is adjusted to ',...
            'the number of nodes given in file (%d)'], size(nodeInitialPosition, 1));
    end
    
end

if mobilitySwitch == 1
    nodeVelocitiesTemp = nodeVelocities;
    clear nodeVelocities;
    nodeVelocities = nodeVelocitiesTemp(1:numberOfNodes, :);
else
    clear nodeVelocities;
    nodeVelocities = zeros(numberOfNodes, 3);
    
    if numberOfTimeDivisions ~= 1
        numberOfTimeDivisions = 1;
        warning('Changing numberOfTimeDivisions to %d', numberOfTimeDivisions)
    end
    
    if paraCfg.totalTimeDuration ~= 0
        paraCfg.totalTimeDuration = 0;
        warning('Changing totalTimeDuration to %d', paraCfg.totalTimeDuration)
    end
end

%% PAA init
iterateNumberOfNodes = 1;
nodePolarization       = zeros(iterateNumberOfNodes, 2);
nodePaaInitialPosition = cell(numberOfNodes,1); %PAA vector position w.r.t node center
nodePaaOrientation     = cell(numberOfNodes,1); %PAA vector position w.r.t node center

while iterateNumberOfNodes <= numberOfNodes
    nodePolarization(iterateNumberOfNodes, :) = [1, 0];
    paaFile = fullfile(inputPath, sprintf('node%dpaa.dat', iterateNumberOfNodes-1));

    % If nodeXpaa.dat is defined
    if isfile(paaFile)
        nodePaaInfo =  readmatrix(paaFile);

        if isempty(nodePaaInfo)
            nodePaaInitialPosition{iterateNumberOfNodes}  = nodePaaInfo;

        else
            nodePaaInitialPosition{iterateNumberOfNodes}  = nodePaaInfo(:, 1:3);

        end

        % PAA orientation not defined
        if size(nodePaaInfo, 2) == 3
            nodePaaOrientation{iterateNumberOfNodes}  = zeros(size(nodePaaInfo));
        
        % PAA orientation
        elseif size(nodePaaInfo, 2) == 6
            nodePaaOrientation{iterateNumberOfNodes}  = nodePaaInfo(:,4:6);

        end

    % If nodeXpaa.dat is not defined
    else
        nodePaaInitialPosition{iterateNumberOfNodes} = zeros(1,3);
        nodePaaOrientation{iterateNumberOfNodes}  = zeros(1,3);

    end
    
    if switchRandomization == 1 && iterateNumberOfNodes > 0
        xCoordinateRandomizer = rand * 8 + 1;
        yCoordinateRandomizer = rand * 17 + 1;
        zCoordinateRandomizer = 1.6;
        nodeInitialPosition(iterateNumberOfNodes, :) =...
            [xCoordinateRandomizer, yCoordinateRandomizer, zCoordinateRandomizer];
        
        xCoordinateRandomizer = rand * 0.7;
        yCoordinateRandomizer = sqrt((0.7^2) - (xCoordinateRandomizer^2));
        zCoordinateRandomizer = 0;
        nodeVelocities(iterateNumberOfNodes, :) =...
            [xCoordinateRandomizer, yCoordinateRandomizer, zCoordinateRandomizer];
        nodeOrientation = zeros(size(nodeInitialPosition));
        nodePositionTime(1, :,iterateNumberOfNodes) = nodeInitialPosition(iterateNumberOfNodes, :);

    end
    
    iterateNumberOfNodes = iterateNumberOfNodes + 1;
end

if switchRandomization ~=0
    switchRandomization = 0;
    warning('Changing switchRandomization to %d', switchRandomization)
end

% If mobility type = 1 nodeTimePosition is generated in raytracer at each
% time step. 
if mobilityType == 1 && mobilitySwitch
    nodePositionTime = [];
    nodePositionTime(1,:,:) = nodeInitialPosition.';
end

%% Process PAA position

paaInfo  = cluster_paa(nodePositionTime, nodePaaInitialPosition, nodePaaOrientation);

%% Output
switchRandomization = 0;

% Check Temp Output Folder
rmdirStatus = rmdir(fullfile(scenarioNameStr, 'Output'), 's');

mkdir(fullfile(scenarioNameStr, 'Output'));
mkdir(fullfile(scenarioNameStr, 'Output/Ns3'));
mkdir(fullfile(scenarioNameStr, 'Output/Visualizer'));

warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames')

paraCfg.mobilityType = mobilityType;
paraCfg.numberOfNodes = numberOfNodes;
paraCfg.numberOfTimeDivisions = numberOfTimeDivisions;
paraCfg.switchRandomization = switchRandomization;

nodeCfg.nodeLoc = nodeInitialPosition;
nodeCfg.nodeAntennaOrientation = nodePaaOrientation;
nodeCfg.nodeOrientation = nodeOrientation;

nodeCfg.nodePolarization = nodePolarization;
nodeCfg.nodePosition = reshape(nodePositionTime, [], 3, numberOfNodes);
if isempty(nodeRotationTime)
    nodeCfg.nodeRotation = zeros(size(nodeCfg.nodeLoc));
else
    nodeCfg.nodeRotation = nodeRotationTime;
end
nodeCfg.nodeVelocities = nodeVelocities;
nodeCfg.paaInfo = paaInfo;
paraCfg.isInitialOrientationOn = any(cellfun(@(x) any(reshape(x, [],1)), nodePaaOrientation));
paraCfg.isDeviceRotationOn = any(nodeRotationTime(:));
paraCfg.isPaaCentered = ~any(cellfun(@(x) any(x(:)), nodePaaInitialPosition));
end
