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

scenarioNameStr = paraCfg.inputScenarioName;
% Input Parameters to be Updated
% scenarioNameStr = paraCfg.scenarioNameStr;
mobilitySwitch = paraCfg.mobilitySwitch;
mobilityType = paraCfg.mobilityType;
numberOfNodes = paraCfg.numberOfNodes;
numberOfTimeDivisions = paraCfg.numberOfTimeDivisions;
switchRandomization = paraCfg.switchRandomization;

% List of paths
inputPath = fullfile(scenarioNameStr, 'Input');
nodesPositionPath = fullfile(scenarioNameStr, 'Output/Ns3/NodesPosition');
paaPositionPath = strcat(scenarioNameStr,'/Output/Ns3/PAAPosition');
paaPositionPathVisual = strcat(scenarioNameStr,'/Output/Visualizer/PAAPosition');

%% Code
nodePosition = [];
nodeEuclidian =[];

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
    try
        nodeLoc = csvread(fullfile(inputPath, 'nodes.dat'));
        nodeVelocities = csvread(fullfile(inputPath, 'nodeVelocities.dat'));
    catch
        switchRandomization = 1;
        warning(['Unable to read nodes.dat or nodeVelocities.dat. ',...
            'Changing switchRandomization to %d'], switchRandomization)
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
    numberOfNodes = size(nodeLoc, 1);
    
    if mobilitySwitch == 1 && mobilityType == 1
        nodeVelocitiesTemp = nodeVelocities;
        clear nodeVelocities;
        nodeVelocities = nodeVelocitiesTemp(1:numberOfNodes, :);
    else
        clear nodeVelocities;
        nodeVelocities = zeros(numberOfNodes, 3);
    end
    
    if mobilityType == 2
        listing = dir(fullfile(scenarioNameStr, 'Input'));
        nodePosition = zeros(paraCfg.numberOfTimeDivisions,3, paraCfg.numberOfNodes);
        nodeEuclidian= zeros(paraCfg.numberOfTimeDivisions,3, paraCfg.numberOfNodes);
        countListing = 0;
        %         for iterateSizeListing = 1:size(listing, 1)
        %             ln = listing(iterateSizeListing).name;
        
        for iterateNumberOfNodes = 1:numberOfNodes
            if sum(arrayfun(@(x) strcmp(x.name,['node',num2str(iterateNumberOfNodes-1),'mobility.mat']), listing))
                savePositionFromTrace(fullfile(inputPath, sprintf('node%dmobility.mat', iterateNumberOfNodes-1)),...
                    fullfile(inputPath,sprintf('Node%dPosition.dat', iterateNumberOfNodes-1)));
                saveEuclidianFromTrace(fullfile(inputPath, sprintf('node%dmobility.mat', iterateNumberOfNodes-1)),...
                    fullfile(inputPath,sprintf('Node%dEuclidian.dat', iterateNumberOfNodes-1)));
            else
                writematrix(nodeLoc(iterateNumberOfNodes,:), fullfile(inputPath,sprintf('Node%dPosition.dat', iterateNumberOfNodes-1)));
                writematrix([0 0 0], fullfile(inputPath,sprintf('Node%dEuclidian.dat', iterateNumberOfNodes-1)));
            end
            % %                 if strcmp(ln, sprintf('node%dmobility.mat', iterateNumberOfNodes))
            %                     savePositionFromTrace(fullfile(inputPath, sprintf('node%dmobility.mat', iterateNumberOfNodes-1)),...
            %                         fullfile(inputPath,sprintf('Node%dPosition.dat', iterateNumberOfNodes)));
            %                     saveEuclidianFromTrace(fullfile(inputPath, sprintf('node%dmobility.mat', iterateNumberOfNodes-1)),...
            %                         fullfile(inputPath,sprintf('Node%dEuclidian.dat', iterateNumberOfNodes)));
            %                 end
        end
        %         end
        listing = dir(fullfile(scenarioNameStr, 'Input'));

        for iterateSizeListing = 1:size(listing, 1)
            ln = listing(iterateSizeListing).name;
            
            for iterateNumberOfNodes = 1:numberOfNodes
                if strcmp(ln, sprintf('Node%dPosition.dat', iterateNumberOfNodes-1))
                    nodePositionTemp = load(fullfile(inputPath, ln));
                    
                    try
                        numberTracePoints = min(size(nodePositionTemp,1),paraCfg.numberOfTimeDivisions);
                        nodePosition(1:numberTracePoints, :, iterateNumberOfNodes) = nodePositionTemp(1:numberTracePoints,:);
                        nodePosition(numberTracePoints+1:end, :, iterateNumberOfNodes) = repmat(nodePositionTemp, [paraCfg.numberOfTimeDivisions-numberTracePoints,1,1]);
                        countListing = countListing + 1;
                    catch
                        mobilityType = 1;
                        warning('Node Position input incorrect. Changing mobilityType to 1');
                    end
                    
                end
                
                if strcmp(ln, sprintf('Node%dEuclidian.dat', iterateNumberOfNodes-1))
                    nodeEuclidianTemp = load(fullfile(inputPath, ln));
                    
                    try
                        numberTracePoints =  min(size(nodeEuclidianTemp,1),paraCfg.numberOfTimeDivisions);
                        nodeEuclidian(1:numberTracePoints, :, iterateNumberOfNodes) = nodeEuclidianTemp(1:numberTracePoints,:);
                        nodeEuclidian(numberTracePoints+1:end, :, iterateNumberOfNodes) = repmat(nodeEuclidianTemp, [paraCfg.numberOfTimeDivisions-numberTracePoints,1,1]);
                        countListing = countListing + 1;
                    catch
                        mobilityType = 1;
                        warning('Node Position input incorrect. Changing mobilityType to 1');
                    end
                    
                end
                
                
            end
        end
        
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

iterateNumberOfNodes = 1;
nodeAntennaOrientation = zeros(numberOfNodes, 3, 3);
nodePolarization = zeros(iterateNumberOfNodes, 2);
nodePAA_position         = cell(numberOfNodes,1); %PAA vector position w.r.t node center

while iterateNumberOfNodes <= numberOfNodes
    nodeAntennaOrientation(iterateNumberOfNodes, :, :) = [1, 0, 0; 0, 1, 0; 0, 0, 1];
    nodePolarization(iterateNumberOfNodes, :) = [1, 0];
    if isfile(strcat(inputPath, [filesep 'node'], num2str(iterateNumberOfNodes), 'PAAs_position.dat' ))
        nodePAA_position{iterateNumberOfNodes} = load(strcat(inputPath, [filesep 'node'], num2str(iterateNumberOfNodes), 'PAAs_position.dat' ));
    else
        nodePAA_position{iterateNumberOfNodes} = [];
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
    end
    
    iterateNumberOfNodes = iterateNumberOfNodes + 1;
end

if switchRandomization ~=0
    switchRandomization = 0;
    warning('Changing switchRandomization to %d', switchRandomization)
end

%% Process PAA position
if isempty(nodePosition)
    [PAA_info]  = cluster_paa(nodeLoc, nodePAA_position);
else
    [PAA_info]  = cluster_paa(permute(nodePosition, [1 3 2]), nodePAA_position);
end
% nodePosition
switchRandomization = 0;

% Check Temp Output Folder
rmdirStatus = rmdir(fullfile(scenarioNameStr, 'Output'), 's');

mkdir(fullfile(scenarioNameStr, 'Output'));
mkdir(fullfile(scenarioNameStr, 'Output/Ns3'));
mkdir(fullfile(scenarioNameStr, 'Output/Visualizer'));

sizeNode = size(nodeLoc);

if ~isfolder(nodesPositionPath)
    mkdir(nodesPositionPath)
end

if ~isfolder(paaPositionPath)
    mkdir(paaPositionPath)
end

if ~isfolder(paaPositionPathVisual)
    mkdir(paaPositionPathVisual)
end
if paraCfg.jsonOutput == 1
    fNodePosition = fopen(fullfile(nodesPositionPath, 'NodesPosition.json'), 'w');
    for i = 1:numberOfNodes
        s = struct('Node', i-1, 'Position', nodeLoc(i,:));
        json = jsonencode(s);
        fprintf(fNodePosition, '%s\n', json);
    end
    fclose(fNodePosition);
else
    writematrix(nodeLoc,fullfile(nodesPositionPath, 'NodesPosition.csv'));
end

if paraCfg.jsonOutput == 1
    fpp = fopen(fullfile(paaPositionPath,'PAAPosition.json'),'w');
    for i = 1:length(nodePAA_position)
        for ipaa = 1:size(nodePAA_position{i},1)
            s = struct('Node', i-1, 'PAA', ipaa-1, 'posShift',  nodePAA_position{i}(ipaa, :));
            json = jsonencode(s);
            fprintf(fpp, '%s\n', json);
        end
    end
    fclose(fpp);
else
    for i = 1:length(nodePAA_position)
        csvwrite(strcat(paaPositionPath, filesep,...
            'Node', num2str(i) ,'PAAPosition.csv'), nodePAA_position{i});
    end
end

if paraCfg.jsonOutput == 1
    fPaa = fopen(strcat(paaPositionPathVisual, filesep,'PAAPosition.json'), 'w');
    for i = 1:length(nodePAA_position)
        for paaId = 1:PAA_info{i}.nPAA_centroids
            s = struct('Node', i-1, 'PAA',paaId-1, 'Position', squeeze(PAA_info{i}.centroid_position(:,paaId,:)),'Rotation', nodeEuclidian(1:numberTracePoints,:,i));
            json = jsonencode(s);
            fprintf(fPaa, '%s\n', json);
        end
    end
    fclose(fPaa);
else
    for i = 1:length(nodePAA_position)
        writematrix([reshape(squeeze(PAA_info{i}.centroid_position), [], 3), nodeEuclidian(1:numberOfTimeDivisions,:,i)] ,strcat(paaPositionPathVisual, filesep,...
            'Node', num2str(i-1) ,'PAAPosition.csv') );
        %  writematrix([squeeze(PAA_info{i}.centroid_position), ...
        %      permute(repmat(nodeEuclidian(1:numberTracePoints,:,i), [1 1 PAA_info{i}.nPAA_centroids]),[ 1, 3,2])] ,...
        %      strcat(paaPositionPathVisual, filesep,...
        %             'Node', num2str(i-1) ,'PAAPosition.csv') );
        
    end
end


warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames')

paraCfg.mobilityType = mobilityType;
paraCfg.numberOfNodes = numberOfNodes;
paraCfg.numberOfTimeDivisions = numberOfTimeDivisions;
paraCfg.switchRandomization = switchRandomization;

% nodeCfg.nodeLoc = cell2mat(reshape(...
%     arrayfun(@(x) x{:}.centroid_position, PAA_info, 'UniformOutput', false),...
%      [],1));
nodeCfg.nodeLoc = cell2mat(...
    arrayfun(@(x) x{:}.centroid_position, PAA_info, 'UniformOutput', false));

nodeCfg.nodeAntennaOrientation = nodeAntennaOrientation;
nodeCfg.nodePolarization = nodePolarization;
nodeCfg.nodePosition = nodePosition;
if isempty(nodeEuclidian)
    nodeCfg.nodeEuclidian = zeros(size(nodeCfg.nodeLoc));
else
    nodeCfg.nodeEuclidian = nodeEuclidian;
end
nodeCfg.nodeVelocities = nodeVelocities;
nodeCfg.PAA_info  = PAA_info;
% nodeCfg.nodePAAInfo = nodePAAInfo;
% nodeCfg.nodePAA_isSmallScaleIndependent = channelGenerationMethod;
% nodeCfg.nodePAA_node_idxs = node_idxs;
% nodeCfg.nodePAA_node_iid_ss = node_iid_ss;
end
