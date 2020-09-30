function [outputPath] = Raytracer(paraCfgInput, nodeCfgInput)
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


%% Input Parameters Management and preallocation
%nodeLocCellArray is a cell array. The i-th entry of the array represents
%the position each PAA in node i. Time x PAA x 3 (coordinates)
% nodeLocCellArray = cellfun(@(x) reshape(squeeze(x.PAA_loc), size(nodeCfgInput.nodeLoc,1), [], 3), ...
%     nodeCfgInput.PAA_info, 'UniformOutput', 0).';
% nodeLocCellArray = cellfun(@(x) x.centroidTimePosition(1,:,:),nodeCfgInput.PAA_info,'UniformOutput',false);
% nodeLoc = [];
% for i = 1:length(nodeLocCellArray)
%     nodeLoc= cat(2, nodeLoc, nodeLocCellArray{i});
% end
nodeLoc(1,:,:) = nodeCfgInput.nodeLoc;

nodePosition = nodeCfgInput.nodePosition;
nodeVelocities = nodeCfgInput.nodeVelocities;
nPAA_centroids = cellfun(@(x) x.nPAA_centroids ,nodeCfgInput.PAA_info );
Mpc = cell(paraCfgInput.numberOfNodes,...
    max(nPAA_centroids),...
    paraCfgInput.numberOfNodes,...
    max(nPAA_centroids),...
    paraCfgInput.totalNumberOfReflections+1,...
    paraCfgInput.numberOfTimeDivisions+1 );
frmRotMpInfo = cell(1, paraCfgInput.totalNumberOfReflections+1);
keepBothQDOutput =0; % If 1 when JSON output will keep also previous QD output

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
end

% Subfolders creation
if ~isfolder(qdFilesPath)
    mkdir(qdFilesPath)
end

% Init output files
if ~paraCfgInput.jsonOutput || keepBothQDOutput
    fids = getQdFilesIds(qdFilesPath, paraCfgInput.numberOfNodes,...
        paraCfgInput.useOptimizedOutputToFile);
end

%% Init
Tx = reshape(squeeze(nodeLoc(1,1,:)), [],3);
Rx = reshape(squeeze(nodeLoc(1,2,:)), [],3);
vtx = nodeVelocities(1,:);
vrx = nodeVelocities(2,:);
switchPolarization = 0;
switchCp = 0;
polarizationTx = [1, 0];
polarizationRx = [1, 0];

% Define Material library for various environments
switch(paraCfgInput.environmentFileName)
    case 'DataCenter.xml'
        MaterialLibrary = importMaterialLibrary('raytracer/Material_library_DataCenter.txt');
    otherwise
        MaterialLibrary = importMaterialLibrary('raytracer/Material_library_Default.txt');
        warning('Environment file ''%s'' not recognized. Using default material library.',...
            paraCfgInput.environmentFileName)
end

% Extracting CAD file and storing in an XMl file, CADFile.xml
[CADop, switchMaterial] = getCadOutput(paraCfgInput.environmentFileName,...
    inputPath, MaterialLibrary, paraCfgInput.referencePoint,...
    paraCfgInput.selectPlanesByDist, paraCfgInput.indoorSwitch);

if paraCfgInput.switchSaveVisualizerFiles == 1
    % Save output file with room coordinates for visualization
    RoomCoordinates = CADop(:, 1:9);
    csvwrite(fullfile(visualizerPath, 'RoomCoordinates.csv'),...
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
        disp([fprintf('%2.2f', iterateTimeDivision/paraCfgInput.numberOfTimeDivisions*100),'%'])
    end
    
    % Update Linear mobility
    if paraCfgInput.mobilityType == 1 && paraCfgInput.mobilitySwitch ==1
        if paraCfgInput.numberOfNodes == 2
            [nodeLoc, Tx, Rx, vtx, vrx, nodeVelocities,nodeCfgInput.PAA_info] = LinearMobility...
                (paraCfgInput.numberOfNodes, paraCfgInput.switchRandomization, ...
                iterateTimeDivision-1, nodeLoc, nodeVelocities, vtx,...
                vrx,TxInitial, RxInitial, timeDivisionValue,...
                CADop, Tx, Rx, nodeCfgInput.PAA_info);
        else
            [nodeLoc, Tx, Rx, vtx, vrx, nodeVelocities,nodeCfgInput.PAA_info] = LinearMobility...
                (paraCfgInput.numberOfNodes, paraCfgInput.switchRandomization,...
                iterateTimeDivision-1, nodeLoc, nodeVelocities,...
                [], [], TxInitial, RxInitial, timeDivisionValue, ...
                CADop, Tx, Rx, nodeCfgInput.PAA_info);
        end
        nodePosition(iterateTimeDivision,:,:) = permute(nodeLoc, [1 3 2]);
    elseif paraCfgInput.mobilityType == 2  && paraCfgInput.mobilitySwitch ==1
        [nodeLoc, nodeVelocities] = NodeExtractor...
            (paraCfgInput.numberOfNodes,  paraCfgInput.switchRandomization, ...
            iterateTimeDivision, nodeLoc, nodeVelocities,...
            nodeCfgInput.nodePosition, timeDivisionValue);
    end
         
    % Compute rotation
    for nodeId = 1:paraCfgInput.numberOfNodes
            centerRotation = nodePosition(iterateTimeDivision,:, nodeId);
            nodeRotationEucAngles = nodeCfgInput.nodeRotation(iterateTimeDivision,:, nodeId);
            paaInitialPosition = reshape(squeeze(...
                nodeCfgInput.PAA_info{nodeId}.centroidTimePosition(iterateTimeDivision,:,:)), [], 3);
            [paaRotatedPosition, nodeEquivalentRotationAngle] = coordinateRotation(paaInitialPosition, ...
                centerRotation,...
                nodeRotationEucAngles ...
                );
            nodeCfgInput.nodeRotationTot(iterateTimeDivision,:, nodeId) = nodeEquivalentRotationAngle;
            nodeCfgInput.PAA_info{nodeId}.centroid_position_rot(iterateTimeDivision,:,:) =paaRotatedPosition;
            isPAACentroidValid(RoomCoordinates,paaRotatedPosition);
    end
    
    % save NodePositionsTrc
    if paraCfgInput.switchSaveVisualizerFiles && ~paraCfgInput.jsonOutput
        filename = sprintf('NodePositionsTrc%d.csv', iterateTimeDivision-1);
        writematrix([squeeze(nodePosition(iterateTimeDivision,:,:)).'...
            squeeze(nodeCfgInput.nodeRotationTot(iterateTimeDivision,:, :)).'],fullfile(visualizerPath, filename)...
            );
    end
    
    % Iterates through all the PAA centroids
    for iterateTx = 1:paraCfgInput.numberOfNodes
        for iterateRx = iterateTx+1:paraCfgInput.numberOfNodes
            for iteratePaaTx = 1:nPAA_centroids(iterateTx)
                for iteratePaaRx = 1:nPAA_centroids(iterateRx)
                    output = [];
                    if (paraCfgInput.numberOfNodes >= 2 || paraCfgInput.switchRandomization == 1)
                        
                        Tx = squeeze(nodeCfgInput.PAA_info{iterateTx}.centroid_position_rot(iterateTimeDivision,iteratePaaTx,:)).';
                        Rx = squeeze(nodeCfgInput.PAA_info{iterateRx}.centroid_position_rot(iterateTimeDivision,iteratePaaRx,:)).';
                        
                        % Rotation Tx struct
                        QTx.center(1,:) = nodePosition(iterateTimeDivision,:,iterateTx);%nodeCfgInput.PAA_info{iterateTx}.node_centroid(iterateTimeDivision,:,:);
                        QTx.angle(1,:) = nodeCfgInput.nodeRotation(iterateTimeDivision,:, iterateTx);
                        
                        % Rotation Rx struct
                        QRx.center(1,:) = nodePosition(iterateTimeDivision,:,iterateRx);%nodeCfgInput.PAA_info{iterateRx}.node_centroid(iterateTimeDivision,:,:);            
                        QRx.angle(1,:) = nodeCfgInput.nodeRotation(iterateTimeDivision,:, iterateRx);
                        
                        vtx = nodeVelocities(iterateTx, :);
                        vrx = nodeVelocities(iterateRx, :);
                    end
  
                    % LOS Path generation
                    [switchLOS, output, frmRotMpInfo{1}] = LOSOutputGenerator(CADop, Rx, Tx,...
                        output, vtx, vrx, switchPolarization, switchCp,...
                        polarizationTx, paraCfgInput.carrierFrequency, 'qTx', QTx, 'qRx', QRx);
                    
                    if paraCfgInput.switchSaveVisualizerFiles && switchLOS
                        multipath1 = [Tx, Rx];
                        if ~paraCfgInput.jsonOutput
                            filename = sprintf('MpcTx%dPAA%dRx%dPAA%dRefl%dTrc%d.csv',...
                                iterateTx-1, iteratePaaTx-1, iterateRx-1, iteratePaaRx-1, 0, iterateTimeDivision-1);
                            csvwrite(fullfile(visualizerPath, filename),...
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
                            count, ~, frmRotMpInfo{iterateOrderOfReflection+1}] = multipath(...
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
                                csvwrite(fullfile(visualizerPath, filename),...
                                    multipath1);
                            else
                                Mpc{iterateTx,iteratePaaTx,iterateRx,iteratePaaRx, iterateOrderOfReflection+1, iterateTimeDivision+1} =multipath1;
                            end
                            
                        end
                        if size(output,1)==1
                            output = repmat(output, 1,1,Nrealizations);
                        end
                        if size(output) > 0
                            output = [output;outputTemporary]; %#ok<AGROW>
                        elseif size(outputTemporary) > 0
                            output = outputTemporary;
                        end
                        
                    end
                    
                    % Create outputPAA array of struct. Each entry of the
                    % array is a struct relative to a NodeTx-NodeRx 
                    % combination. Each struct has the entries 
                    % - paaTxXXpaaRxYY: channel between paaTx XX and paaRx
                    % YY.
                    % -frmRotMpInfopaaTxXXpaaRxXX. Information to perform
                    % frame rotation after raytracing                    
                    eval(['outputPAA{iterateTx, iterateRx}.frmRotMpInfopaaTx',num2str(iteratePaaTx-1),'paaRx', num2str(iteratePaaRx-1), '= [frmRotMpInfo{:}];'] );
                    eval(['outputPAA{iterateRx, iterateTx}.frmRotMpInfopaaTx',num2str(iteratePaaRx-1),'paaRx', num2str(iteratePaaTx-1), '= reverseFrmRotMpInfo([frmRotMpInfo{:}]);'] );
                    frmRotMpInfo = {};                    
                    eval(['outputPAA{iterateTx, iterateRx}.paaTx',num2str(iteratePaaTx-1),'paaRx', num2str(iteratePaaRx-1), '= output;'] );
                    eval(['outputPAA{iterateRx, iterateTx}.paaTx',num2str(iteratePaaRx-1),'paaRx', num2str(iteratePaaTx-1), '= reverseOutputTxRx(output);'] );
                    
                end
            end
            
        end
    end
    
    outputPAATime(:,:,iterateTimeDivision) = generateChannelPaa(outputPAA, nodeCfgInput.PAA_info);  %#ok<AGROW>
    
    % Write QD output in CSV files
    if ~paraCfgInput.jsonOutput || keepBothQDOutput
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
    end
    clear outputPAA
end

% Write QD output in JSON files
if paraCfgInput.jsonOutput
    writeQdJsonOutput(outputPAATime,cellfun(@(x) x.nPaa,  nodeCfgInput.PAA_info),...
        qdFilesPath,...
        paraCfgInput.qdFilesFloatPrecision);
    
    Mpc(:,:,:,:,:,1) = [];
end

if paraCfgInput.switchSaveVisualizerFiles && paraCfgInput.jsonOutput
    %% Write MPC.json
    f = fopen(fullfile(visualizerPath, 'Mpc.json'), 'w');
    for iterateTx = 1:paraCfgInput.numberOfNodes
        for iterateRx = iterateTx+1:paraCfgInput.numberOfNodes
            for iteratePaaTx = 1:nPAA_centroids(iterateTx)
                nodeTxCluster  = nodeCfgInput.PAA_info{iterateTx}.paaInCluster{iteratePaaTx};
                for txPaaCluster = 1:length(nodeTxCluster)
                    for iteratePaaRx = 1:nPAA_centroids(iterateRx)
                        nodeRxCluster  = nodeCfgInput.PAA_info{iterateRx}.paaInCluster{iteratePaaRx};
                        for rxPaaCluster = 1:length(nodeRxCluster)
                            for reflOrd = 1:paraCfgInput.totalNumberOfReflections+1
                                Mpc_t = squeeze((Mpc(iterateTx,iteratePaaTx,...
                                    iterateRx,iteratePaaRx,reflOrd,:)));
                                s = struct('TX', iterateTx-1, 'PAA_TX', nodeTxCluster(txPaaCluster)-1,...
                                    'RX', iterateRx-1, 'PAA_RX', nodeRxCluster(rxPaaCluster)-1, ...
                                    'Rorder', reflOrd-1);
                                s.MPC =  Mpc_t;
                                json = jsonencode(s);
                                fprintf(f, '%s\n', json);
                            end
                        end
                    end
                end
            end
        end
    end
    fclose(f);
    
    %% Write Node Position
    filename = sprintf('NodePositions.json');
    f = fopen(fullfile(visualizerPath, filename), 'w');
    for i = 1:paraCfgInput.numberOfNodes
        s = struct('Node' , i-1, ...
            'Position', [nodePosition(:,:,i); [inf inf inf]], ...
            'Rotation', [nodeCfgInput.nodeRotationTot(:,:,i); [inf inf inf]]);
        json = jsonencode(s); % Add a temporary inf vector to make sure
        % more than a single vector will be encoded. Matlab json
        % encoder lose the square brackets when encoding vectors.
        str2remove =',[null,null,null]'; %Temporary string to remove
        rem_ind_start = num2cell(strfind(json, str2remove)); % Find start string to remove
        index2rm = cell2mat(cellfun(@(x) x:x+length(str2remove)-1,rem_ind_start,'UniformOutput',false)); % Create index of char to remove
        json(index2rm) = []; % Remove temporary vector.
        fprintf(f, '%s\n', json);
    end
    fclose(f);
    
    %% Write PAAPosition.json
    f = fopen(strcat(visualizerPath, filesep,'PAAPosition.json'), 'w');
    for i = 1:paraCfgInput.numberOfNodes
        idOrientation = 0;
        for paaId = 1:nPAA_centroids(i)
            nodeTxCluster  = nodeCfgInput.PAA_info{i}.paaInCluster{paaId};
            for paaCentroid = 1:length(nodeTxCluster)
                idOrientation = idOrientation+1;
                s = struct('Node', i-1, 'PAA',nodeTxCluster(paaCentroid)-1, ...
                    'Centroid', nodeCfgInput.PAA_info{i}.centroids(paaId)-1, ...
                    'Orientation', nodeCfgInput.nodeAntennaOrientation{i}(idOrientation,:) ,...
                'Position', [reshape(squeeze(nodeCfgInput.PAA_info{i}.centroid_position_rot(:,paaId,:)), [],3); [inf inf inf]]);
                json = jsonencode(s);% Add a temporary inf vector to make sure
                % more than a single vector will be encoded. Matlab json
                % encoder lose the square brackets when encoding vectors.
                str2remove =',[null,null,null]'; %Temporary string to remove
                rem_ind_start = num2cell(strfind(json, str2remove)); % Find start string to remove
                index2rm = cell2mat(cellfun(@(x) x:x+length(str2remove)-1,rem_ind_start,'UniformOutput',false)); % Create index of char to remove
                json(index2rm) = []; % Remove temporary vector.
                fprintf(f, '%s\n', json);
            end
        end
    end
    fclose(f); 
end

if ~paraCfgInput.jsonOutput || keepBothQDOutput
    closeQdFilesIds(fids, paraCfgInput.useOptimizedOutputToFile);
end

% Write useful output information. Set to 0 to allow succeful test.
writeReportOutput =0 ; 
if writeReportOutput
    f = fopen(strcat(outputPath, filesep,'report.dat'), 'w'); %#ok<UNRCH>
%     elapsedTime = toc;
%     fprintf(f, 'Elapsed Time:\t%f\n', elapsedTime);    
%     isDeviceRotationOn = 'true';
    fprintf(f, 'Device Rotation:\t%d\n', paraCfgInput.isDeviceRotationOn);
%     isInitialOrientationOn = 'true';
    fprintf(f, 'Initial Orientation:\t%d\n', paraCfgInput.isInitialOrientationOn);
%     isPAAcentered = 'true';
    fprintf(f, 'PAA centered:\t%d\n', paraCfgInput.isPAAcentered);
    fclose(f);
end