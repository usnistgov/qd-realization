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
nodeRotation= zeros(paraCfg.numberOfTimeDivisions,3, paraCfg.numberOfNodes);
nodeLoc = zeros( paraCfg.numberOfNodes, 3);

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
% nodes.dat file contains nodes locations and nodeVelocities contains their
% velocities
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
    
   %% Load mobility matrix and write NodePositionX.dat and NodeRotationX.dat   
    listing = dir(fullfile(scenarioNameStr, 'Input'));
    nodePosition = zeros(paraCfg.numberOfTimeDivisions,3, paraCfg.numberOfNodes);
    countListing = 0;
    
    for iterateNumberOfNodes = 1:numberOfNodes
        % If mobility matrix is present
        if sum(arrayfun(@(x) strcmp(x.name,['node',num2str(iterateNumberOfNodes-1),'mobility.mat']), listing))
            savePositionFromTrace(fullfile(inputPath, sprintf('node%dmobility.mat', iterateNumberOfNodes-1)),...
                fullfile(inputPath,sprintf('NodePosition%d.dat', iterateNumberOfNodes-1)) );
            saveRotationFromTrace(fullfile(inputPath, sprintf('node%dmobility.mat', iterateNumberOfNodes-1)),...
                fullfile(inputPath,sprintf('NodeRotation%d.dat', iterateNumberOfNodes-1)) );
        % If mobility matrix is not present check for NodePosition and
        % NodeRotation
        else           
            % If only NodePosition is present, initialize NodeRotation to
            % [0,0,0]
            if sum(arrayfun(@(x) strcmp(x.name,['NodePosition',num2str(iterateNumberOfNodes-1),'.dat']), listing)) && ...
                    ~ sum(arrayfun(@(x) strcmp(x.name,['NodeRotation',num2str(iterateNumberOfNodes-1),'.dat']), listing))
                nlines = size(readmatrix(fullfile(inputPath,sprintf('NodePosition%d.dat', iterateNumberOfNodes-1))),1);
                writematrix(repmat([0 0  0], nlines,1), fullfile(inputPath,sprintf('NodeRotation%d.dat', iterateNumberOfNodes-1)));
                warning('NodeRotation%d.dat not present. Rotation set to [0,0,0] for all time instances.', iterateNumberOfNodes-1)
            % If only NodeRotation return error:  NodePosition is needed
            elseif ~ sum(arrayfun(@(x) strcmp(x.name,['NodePosition',num2str(iterateNumberOfNodes-1),'.dat']), listing)) && ...
                    sum(arrayfun(@(x) strcmp(x.name,['NodeRotation',num2str(iterateNumberOfNodes-1),'.dat']), listing))
                error('NodePosition%d.dat is not present.', iterateNumberOfNodes-1);
            % If NodePosition and NodeRotation are present no further
            % actions are needed
            elseif sum(arrayfun(@(x) strcmp(x.name,['NodePosition',num2str(iterateNumberOfNodes-1),'.dat']), listing)) && ...
                    sum(arrayfun(@(x) strcmp(x.name,['NodeRotation',num2str(iterateNumberOfNodes-1),'.dat']), listing))
                
            % If NodePosition and NodeRotation are not present return error
            else
                error('NodePosition%d.dat and NodeRotation%d.dat are not present.', iterateNumberOfNodes-1,  iterateNumberOfNodes-1);
            end
        end
    end

   %%  Load NodePositionX.dat and NodeRotationX.dat   
   for iterateNumberOfNodes = 1:numberOfNodes
       
       % NodePosition processing
       ln = sprintf('NodePosition%d.dat', iterateNumberOfNodes-1);
       nodePositionTemp = load(fullfile(inputPath, ln));
       
       if mobilityType == 1 && mobilitySwitch
           nodeLoc(iterateNumberOfNodes,:) = nodePositionTemp;
           countListing = countListing + 1;

       else
           
           if  size(nodePositionTemp,1)< size(nodePosition,1) &&  ...
                   size(nodePositionTemp,1) > 1
               nodePosition = nodePosition(1:size(nodeRotationTemp,1), :,:);
               paraCfg.numberOfTimeDivisions = size(nodeRotationTemp,1) ;
               numberOfTimeDivisions = paraCfg.numberOfTimeDivisions;
               warning('Time divisition too long.')
           end
           
           try
               numberTracePoints = min(size(nodePositionTemp,1),paraCfg.numberOfTimeDivisions);
               nodePosition(1:numberTracePoints, :, iterateNumberOfNodes) = nodePositionTemp(1:numberTracePoints,:);
               nodePosition(numberTracePoints+1:end, :, iterateNumberOfNodes) = repmat(nodePositionTemp, [paraCfg.numberOfTimeDivisions-numberTracePoints,1,1]);
               nodeLoc(iterateNumberOfNodes,:) = squeeze(nodePosition(1,:,iterateNumberOfNodes));
               countListing = countListing + 1;
           catch
               mobilityType = 1;
               warning('Node Position input incorrect. Changing mobilityType to 1');
           end
       end
       
       % NodeRotation processing
       ln  = sprintf('NodeRotation%d.dat', iterateNumberOfNodes-1);
       nodeRotationTemp = load(fullfile(inputPath, ln));
       if mobilityType == 1 && mobilitySwitch
           nodeOrientation(iterateNumberOfNodes,:) = squeeze(nodeRotation(1,:,iterateNumberOfNodes));
           countListing = countListing + 1;
       else
           
           if  size(nodeRotationTemp,1)< size(nodeRotation,1) &&  ...
                   size(nodeRotationTemp,1) > 1
               nodeRotation = nodeRotation(1:size(nodeRotationTemp,1), :,:);
               paraCfg.numberOfTimeDivisions = size(nodeRotationTemp,1) ;
               numberOfTimeDivisions = paraCfg.numberOfTimeDivisions;
               warning('Time divisition too long.')
           end
           
           try
               numberTracePoints =  min(size(nodeRotationTemp,1),paraCfg.numberOfTimeDivisions);
               nodeRotation(1:numberTracePoints, :, iterateNumberOfNodes) = nodeRotationTemp(1:numberTracePoints,:);
               nodeRotation(numberTracePoints+1:end, :, iterateNumberOfNodes) = repmat(nodeRotationTemp, [paraCfg.numberOfTimeDivisions-numberTracePoints,1,1]);
               nodeOrientation(iterateNumberOfNodes,:) = squeeze(nodeRotation(1,:,iterateNumberOfNodes));
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
        if numberOfTimeDivisions ~= size(nodePosition, 1) - 1
            numberOfTimeDivisions = size(nodePosition, 1) - 1;
            warning('Changing numberOfTimeDivisions to %d', numberOfTimeDivisions)
        end
    end
        
    if size(nodeLoc,1) ~= size(nodeVelocities,1) && mobilitySwitch == 1 && mobilityType==1
        error(['nodes.dat and nodeVelocities.dat do not have same number ',...
            'of rows. Please check the input files in the Input folder.'])
    end
    
    if ~isempty(numberOfNodes) && numberOfNodes ~= size(nodeLoc, 1)
        warning(['"numberOfNodes" parameter does not match the number of ',...
            'nodes given in file. The "numberOfNodes" is adjusted to ',...
            'the number of nodes given in file (%d)'], size(nodeLoc, 1));
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
nodePAA_position       = cell(numberOfNodes,1); %PAA vector position w.r.t node center
nodePAA_Orientation       = cell(numberOfNodes,1); %PAA vector position w.r.t node center

while iterateNumberOfNodes <= numberOfNodes
    nodePolarization(iterateNumberOfNodes, :) = [1, 0];
    % If nodeXpaa.dat is present
    if isfile(strcat(inputPath, [filesep 'node'], num2str(iterateNumberOfNodes-1), 'paa.dat' ))
        nodePAA_info =  load(strcat(inputPath, [filesep 'node'], num2str(iterateNumberOfNodes-1), 'paa.dat' ));
        if isempty(nodePAA_info)
            nodePAA_position{iterateNumberOfNodes}  = nodePAA_info;
        else
            nodePAA_position{iterateNumberOfNodes}  = nodePAA_info(:, 1:3);
        end
        if size(nodePAA_info,2) ==3
            nodePAA_Orientation{iterateNumberOfNodes}  = zeros(size(nodePAA_info));
        elseif size(nodePAA_info,2) ==6
            nodePAA_Orientation{iterateNumberOfNodes}  = nodePAA_info(:,4:6);
        end
    else
        nodePAA_position{iterateNumberOfNodes} = zeros(1,3);
        nodePAA_Orientation{iterateNumberOfNodes}  = zeros(1,3);
    end
    
    if switchRandomization == 1 && iterateNumberOfNodes > 0
        xCoordinateRandomizer = rand * 8 + 1;
        yCoordinateRandomizer = rand * 17 + 1;
        zCoordinateRandomizer = 1.6;
        nodeLoc(iterateNumberOfNodes, :) =...
            [xCoordinateRandomizer, yCoordinateRandomizer, zCoordinateRandomizer];
        
        xCoordinateRandomizer = rand * 0.7;
        yCoordinateRandomizer = sqrt((0.7^2) - (xCoordinateRandomizer^2));
        zCoordinateRandomizer = 0;
        nodeVelocities(iterateNumberOfNodes, :) =...
            [xCoordinateRandomizer, yCoordinateRandomizer, zCoordinateRandomizer];
        nodeOrientation = zeros(size(nodeLoc));
        
    end
    
    iterateNumberOfNodes = iterateNumberOfNodes + 1;
end

if switchRandomization ~=0
    switchRandomization = 0;
    warning('Changing switchRandomization to %d', switchRandomization)
end

%Check if nodePosition has been generated
if ~exist('nodePosition','var') || (mobilityType == 1 && mobilitySwitch)
    nodePosition = nodeLoc.';
end

%% Process PAA position

[PAA_info]  = cluster_paa(nodePosition, nodePAA_position, nodePAA_Orientation);

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

nodeCfg.nodeLoc = cell2mat(...
    arrayfun(@(x) x{:}.centroid_position, PAA_info, 'UniformOutput', false));

nodeCfg.nodeAntennaOrientation = nodePAA_Orientation;
nodeCfg.nodeOrientation = nodeOrientation;

nodeCfg.nodePolarization = nodePolarization;
nodeCfg.nodePosition = reshape(nodePosition, [], 3,numberOfNodes);
if isempty(nodeRotation)
    nodeCfg.nodeRotation = zeros(size(nodeCfg.nodeLoc));
else
    nodeCfg.nodeRotation = nodeRotation;
end
nodeCfg.nodeVelocities = nodeVelocities;
nodeCfg.PAA_info  = PAA_info;
paraCfg.isInitialOrientationOn =  any(cellfun(@(x) any(reshape(x, [],1)), nodePAA_Orientation));
paraCfg.isDeviceRotationOn = any(nodeRotation(:));
paraCfg.isPAAcentered =  ~any(cellfun(@(x) any(x(:)), nodePAA_position));
end
