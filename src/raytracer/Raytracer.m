function [outputPath, varargout] = Raytracer(paraCfgInput, nodeCfgInput)
% Inputs:
% RootFolderPath - it is the current location of the folder where the function is called from
% environmentFileName - it is the CAD file name
% switchRandomization - boolean to either randomly generates nodes and velocity or not
% mobilitySwitch -  is boolean to either have mobility or not
% totalNumberOfReflections - is the highest order of reflections to be computed
% switchQDGenerator - Switch to turn ON or OFF the Qausi dterministic module 1 = ON, 0 = OFF
% nodeLoc - 2d array which contains all node locations
% nodeVelocities - 2d array which contains all node velocities
% nodePolarization - 2d array which contains all node polarization
% nodeAntennaOrientation - 2d array which contains all node antenna orientation
% totalTimeDuration, n1 are for granularity in time domain. t is total period and n is the
% number of divisions of that time period
% mobilityType - This switch lets the user to decide the input to mobility
% 1 = Linear, 2 = input from File
% nodePosition - these are positions of nodes in a 2D array which are
% extracted from a file
% indoorSwitch - This boolean lets user say whether the given CAD file
% is indoor or outdorr. If indoor, then the value is 1 else the value is 0.
% generalizedScenario - This boolean lets user say whether a scenario
% conforms to a regular indoor or outdoor environment or it is a more
% general scenario.
% selectPlanesByDist - This is selection of planes/nodes by distance.
% r = 0 means that there is no limitation.
% referencePoint - Reference point is the center of limiting sphere
%
% Outputs:
% N/A


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
% anddistributing the software and you assume all risks associated with its use, including
% but not limited to the risks and costs of program errors, compliance with applicable
% laws, damage to or loss of data, programs or equipment, and the unavailability or
% interruption of operation. This software is not intended to be used in any situation
% where a failure could cause risk of injury or damage to property. The software
% developed by NIST employees is not subject to copyright protection within the United
% States.
%
% Modified by: Mattia Lecci <leccimat@dei.unipd.it>, Refactored code


%% Input Parameters Management
numberOfNodes = paraCfgInput.numberOfNodes;

% nodeLoc = nodeCfgInput.nodeLoc;
% nodeLoc = reshape(...
%     cell2mat(cellfun(@(x) x.PAA_loc(1,1:3), nodeCfgInput.PAA_info, 'UniformOutput', 0)),...
%     3, []).';
% nodeLoc = reshape(...
%     cell2mat(cellfun(@(x) x.PAA_loc, nodeCfgInput.PAA_info, 'UniformOutput', 0)),...
%     3, []).';
nodeLoc = [];
nodeLocCellArray = cellfun(@(x) reshape(squeeze(x.PAA_loc), size(nodeCfgInput.nodeLoc,1), [], 3), nodeCfgInput.PAA_info, 'UniformOutput', 0).';
for i = 1:length(nodeLocCellArray)
    nodeLoc= cat(2, nodeLoc, nodeLocCellArray{i});
end
nodeAntennaOrientation = nodeCfgInput.nodeAntennaOrientation;
nodePolarization = nodeCfgInput.nodePolarization;
nodePosition = nodeCfgInput.nodePosition;
nodeVelocities = nodeCfgInput.nodeVelocities;
nPAA_centroids = cellfun(@(x) x.nPAA_centroids ,nodeCfgInput.PAA_info );

% Input checking
if paraCfgInput.switchQDGenerator == 1 &&...
        paraCfgInput.carrierFrequency ~= 60e9
    warning(['Please, note that diffuse scattering model is only ',...
        'valid for fc=60 GHz'])
end

% List of paths
inputPath = fullfile(paraCfgInput.inputScenarioName, 'Input');
outputPath = fullfile(paraCfgInput.inputScenarioName, 'Output');

ns3Path = fullfile(outputPath, 'Ns3');
qdFilesPath = fullfile(ns3Path, 'QdFiles');

if paraCfgInput.switchSaveVisualizerFiles == 1
    visualizerPath = fullfile(outputPath, 'Visualizer');
    
    nodePositionsPath = fullfile(visualizerPath, 'NodePositions');
    roomCoordinatesPath = fullfile(visualizerPath, 'RoomCoordinates');
    mpcCoordinatesPath = fullfile(visualizerPath, 'MpcCoordinates');
end

% Subfolders creation
if ~isfolder(qdFilesPath)
    mkdir(qdFilesPath)
end

if paraCfgInput.switchSaveVisualizerFiles == 1
    
    if ~isfolder(nodePositionsPath)
        mkdir(nodePositionsPath)
    end
    if ~isfolder(roomCoordinatesPath)
        mkdir(roomCoordinatesPath)
    end
    if ~isfolder(mpcCoordinatesPath)
        mkdir(mpcCoordinatesPath)
    end
    
end

% Init output files
fids = getQdFilesIds(qdFilesPath, paraCfgInput.numberOfNodes,...
    paraCfgInput.useOptimizedOutputToFile);

%% Init
Tx = reshape(squeeze(nodeLoc(1,1,:)), [],3);
Rx = reshape(squeeze(nodeLoc(1,2,:)), [],3);
vtx = nodeVelocities(1,:);
vrx = nodeVelocities(2,:);
switchPolarization = 0;
switchCp = 0;
polarizationTx = [1, 0];
polarizationRx = [1, 0];

MaterialLibrary = importMaterialLibrary('raytracer/Material_library.txt');

% Extracting CAD file and storing in an XMl file, CADFile.xml
[CADop, switchMaterial] = getCadOutput(paraCfgInput.environmentFileName,...
    inputPath, MaterialLibrary, paraCfgInput.referencePoint,...
    paraCfgInput.selectPlanesByDist, paraCfgInput.indoorSwitch);

if paraCfgInput.switchSaveVisualizerFiles == 1
    % Save output file with room coordinates for visualization
    RoomCoordinates = CADop(:, 1:9);
    csvwrite(fullfile(roomCoordinatesPath, 'RoomCoordinates.csv'),...
        RoomCoordinates);
end


%% Randomization
% if number of nodes is greater than 1 or switch_randomization is set to 1,
% the program generates nodes randomly. If one has more than 2 nodes but
% know the exact locations of nodes, then disable this if statement and
% replace node and node_v with the values of node positions and node
% velocities repsectively

TxInitial = Tx;
RxInitial = Rx;
% t - total time period, n - number of divisions
timeDivisionValue = paraCfgInput.totalTimeDuration / paraCfgInput.numberOfTimeDivisions;

% Finite difference method to simulate mobility. x=x0 + v*dt.
% This method ensures the next position wouldnt collide with any of the
% planes. If that occurs then the velocities are simply reversed (not
% reflected). At every time step the positions of all nodes are updated
for iterateTimeDivision = 1:paraCfgInput.numberOfTimeDivisions
    if mod(iterateTimeDivision,100)==0
    display([num2str(iterateTimeDivision/32070*100),'%'])
    end
    % update mobility
    if paraCfgInput.mobilityType == 1
        if paraCfgInput.numberOfNodes == 2
            [nodeLoc, Tx, Rx, vtx, vrx, nodeVelocities] = LinearMobility...
                (paraCfgInput.numberOfNodes, paraCfgInput.switchRandomization, ...
                iterateTimeDivision-1, nodeLoc, nodeVelocities, vtx,...
                vrx,TxInitial, RxInitial, timeDivisionValue,...
                CADop, Tx, Rx);
        else
            [nodeLoc, Tx, Rx, vtx, vrx, nodeVelocities] = LinearMobility...
                (paraCfgInput.numberOfNodes, paraCfgInput.switchRandomization,...
                iterateTimeDivision-1, nodeLoc, nodeVelocities,...
                [], [], TxInitial, RxInitial, timeDivisionValue, ...
                CADop, Tx, Rx);
        end
        
    elseif paraCfgInput.mobilityType == 2
        [nodeLoc, nodeVelocities] = NodeExtractor...
            (paraCfgInput.numberOfNodes,  paraCfgInput.switchRandomization, ...
            iterateTimeDivision, nodeLoc, nodeVelocities,...
            nodeCfgInput.nodePosition, timeDivisionValue);
    end
    
    % save NodePositionsTrc
    if paraCfgInput.switchSaveVisualizerFiles && ~paraCfgInput.jsonOutput
        filename = sprintf('NodePositionsTrc%d.csv', iterateTimeDivision-1);
        writematrix([squeeze(nodePosition(iterateTimeDivision,:,:)).'...
            squeeze(nodeCfgInput.nodeEuclidian(iterateTimeDivision,:, :)).'],fullfile(nodePositionsPath, filename)...
            );
    end

    % Iterates through all the nodes
    for iterateTx = 1:paraCfgInput.numberOfNodes
        for iterateRx = iterateTx+1:paraCfgInput.numberOfNodes
            for iteratePaaTx = 1:nPAA_centroids(iterateTx)
                for iteratePaaRx = 1:nPAA_centroids(iterateRx)
                    output = [];
                    if (paraCfgInput.numberOfNodes >= 2 || paraCfgInput.switchRandomization == 1)
                        Tx = squeeze(nodeCfgInput.PAA_info{iterateTx}.centroid_position(iterateTimeDivision,iteratePaaTx,:)).';%nodeLoc(iterateTx, :);
                        Rx = squeeze(nodeCfgInput.PAA_info{iterateRx}.centroid_position(iterateTimeDivision,iteratePaaRx,:)).';%nodeLoc(iterateRx, :);
                        QTx.cTx = nodeCfgInput.PAA_info{iterateTx}.node_centroid(iterateTimeDivision,:,:);
                        QTx.euc = nodeCfgInput.nodeEuclidian(iterateTimeDivision,:, iterateTx);
                        QRx.cRx = nodeCfgInput.PAA_info{iterateRx}.node_centroid(iterateTimeDivision,:,:);
                        QRx.euc = nodeCfgInput.nodeEuclidian(iterateTimeDivision,:, iterateRx);
%                         if iterateTimeDivision==1886, keyboard, end
                        [Tx] = pointRotation(Tx, ...
                            QTx.cTx,...
                            QTx.euc);
                        [Rx] = pointRotation(Rx, ...
                            QRx.cRx,...
                            QRx.euc);                               
                        vtx = nodeVelocities(iterateTx, :);
                        vrx = nodeVelocities(iterateRx, :);
                    end
                    
                    % LOS Path generation
                    [switchLOS, output] = LOSOutputGenerator(CADop, Rx, Tx,...
                        output, vtx, vrx, switchPolarization, switchCp,...
                        polarizationTx, paraCfgInput.carrierFrequency, 'qTx', QTx, 'qRx', QRx);

                    if paraCfgInput.switchSaveVisualizerFiles && switchLOS
                        multipath1 = [Tx, Rx];
                        if ~paraCfgInput.jsonOutput
                            filename = sprintf('MpcTx%dPAA%dRx%dPAA%dRefl%dTrc%d.csv',...
                                iterateTx-1, iteratePaaTx-1, iterateRx-1, iteratePaaRx-1, 0, iterateTimeDivision-1);
                            csvwrite(fullfile(mpcCoordinatesPath, filename),...
                                multipath1);
                        else
                            Mpc{iterateTx,iteratePaaTx,iterateRx,iteratePaaRx, 1, iterateTimeDivision+1} =multipath1;
                        end
                    end
                    
                    % Higher order reflections (Non LOS)
                    for iterateOrderOfReflection = 1:paraCfgInput.totalNumberOfReflections
                        numberOfReflections = iterateOrderOfReflection;
                        
                        [ArrayOfPoints, ArrayOfPlanes, numberOfPlanes,...
                            ~, ~, arrayOfMaterials, ~] = treetraversal(CADop,...
                            numberOfReflections, numberOfReflections,...
                            0, 1, 1, 1, Rx, Tx, [], [],...
                            switchMaterial, [], 1, paraCfgInput.generalizedScenario);
                        
                        numberOfPlanes = numberOfPlanes - 1;
                        
                        Nrealizations = nodeCfgInput.PAA_info{iterateTx}.nodePAAInfo{iteratePaaTx}.indep_stoch_channel*...
                            nodeCfgInput.PAA_info{iterateRx}.nodePAAInfo{iteratePaaRx}.indep_stoch_channel;
                        [~, ~, outputTemporary, multipathTemporary,...
                            count, ~] = multipath(...
                            ArrayOfPlanes, ArrayOfPoints, Rx, Tx, ...
                            CADop, numberOfPlanes, ...
                            MaterialLibrary, arrayOfMaterials, ...
                            switchMaterial, vtx, vrx, ...
                            switchPolarization, polarizationTx, [],...
                            polarizationRx, [], switchCp,...
                            paraCfgInput.switchQDGenerator,...
                            paraCfgInput.carrierFrequency,'indStoc', Nrealizations,...
                            'qTx', QTx, 'qRx', QRx);
                        
                        if paraCfgInput.switchSaveVisualizerFiles &&...
                                size(multipathTemporary,1) > 0
                            
                            multipath1 = multipathTemporary(1:count,...
                                2:size(multipathTemporary,2));
                            if  ~paraCfgInput.jsonOutput
                                filename = sprintf('MpcTx%dPAA%dRx%dPAA%dRefl%dTrc%d.csv',...
                                    iterateTx-1, iteratePaaTx-1,...
                                    iterateRx-1, iteratePaaRx-1,...
                                    iterateOrderOfReflection, iterateTimeDivision-1);
                                csvwrite(fullfile(mpcCoordinatesPath, filename),...
                                    multipath1);
                            else
                                Mpc{iterateTx,iteratePaaTx,iterateRx,iteratePaaRx, iterateOrderOfReflection+1, iterateTimeDivision+1} =multipath1;
                            end

                        end
                        if size(output,1)==1
                            output = repmat(output, 1,1,Nrealizations);
                        end
                        if size(output) > 0
                            output = [output;outputTemporary];
                            multipath1 = multipathTemporary;
                        elseif size(outputTemporary) > 0
                            output = outputTemporary;
                            multipath1 = multipathTemporary;
                        end
                        
                    end
                    eval(['outputPAA{iterateTx, iterateRx}.paaTx',num2str(iteratePaaTx),'paaRx', num2str(iteratePaaRx), '= output;'] );
                    eval(['outputPAA{iterateRx, iterateTx}.paaTx',num2str(iteratePaaRx),'paaRx', num2str(iteratePaaTx), '= reverseOutputTxRx(output);'] );
                end
            end
            
        end
    end

outputPAATime(:,:,iterateTimeDivision) = generateChannelPaa(outputPAA, nodeCfgInput.PAA_info); %#ok<NODEF>

for iterateTx = 1:paraCfgInput.numberOfNodes
    for iterateRx = iterateTx+1:paraCfgInput.numberOfNodes
        writeQdFileOutput(outputPAATime{iterateTx, iterateRx,iterateTimeDivision}, paraCfgInput.useOptimizedOutputToFile, fids, iterateTx, iterateRx,...
            qdFilesPath,...
            paraCfgInput.qdFilesFloatPrecision);
        writeQdFileOutput(outputPAATime{iterateRx,iterateTx,iterateTimeDivision}, paraCfgInput.useOptimizedOutputToFile, fids, iterateRx,iterateTx,...
            qdFilesPath,...
            paraCfgInput.qdFilesFloatPrecision);
    end
end
clear outputPAA
end

if paraCfgInput.jsonOutput
 writeQdJsonOutput(outputPAATime,...
            qdFilesPath,...
            paraCfgInput.qdFilesFloatPrecision);

Mpc(:,:,:,:,:,1) = [];
end
if paraCfgInput.switchSaveVisualizerFiles && paraCfgInput.jsonOutput
    % MPC       
    fmpc = fopen(fullfile(mpcCoordinatesPath, 'Mpc.json'), 'w');
    for iterateTx = 1:paraCfgInput.numberOfNodes
        for iterateRx = iterateTx+1:paraCfgInput.numberOfNodes
            for iteratePaaTx = 1:nPAA_centroids(iterateTx)
                for iteratePaaRx = 1:nPAA_centroids(iterateRx)
                    for reflOrd = 1:paraCfgInput.totalNumberOfReflections+1
                        Mpc_t = squeeze(cell2mat(Mpc(iterateTx,iteratePaaTx,...
                            iterateRx,iteratePaaRx,reflOrd,:)));
                            s = struct('TX', iterateTx-1, 'PAA_TX', iteratePaaTx-1,...
                                   'RX', iterateRx-1, 'PAA_RX', iteratePaaRx-1, ...
                                   'Rorder', reflOrd-1, 'MPC', Mpc_t);
                               json = jsonencode(s);
                               fprintf(fmpc, '%s\n', json);  
                    end
                end
            end
        end
    end
    fclose(fmpc);
    
    
    % Node Position
    if paraCfgInput.switchSaveVisualizerFiles && paraCfgInput.jsonOutput
        filename = sprintf('NodePositions.json');
        f = fopen(fullfile(nodePositionsPath, filename), 'w');
        fprintf(f, '%s', jsonencode(nodePosition));
        fclose(f);
    end
end

closeQdFilesIds(fids, paraCfgInput.useOptimizedOutputToFile);
% varargout{1} = outputPAA;
end