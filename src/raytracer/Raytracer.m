function [outputPath] = Raytracer(paraCfgInput, nodeCfgInput)
%%RAYTRACER generates the QD channel model.
% Inputs:
% paraCfgInput - Simulation configuration
% nodeCfgInput - Node configuration 


%% -------------Software Disclaimer---------------
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
%              Steve Blandino <steve.blandino@nist.gov>

%% Input Parameters Management and preallocation
nodeLoc(1,:,:) = nodeCfgInput.nodeLoc;
nodePosition = nodeCfgInput.nodePosition;
nPAA_centroids = cellfun(@(x) x.nPAA_centroids ,nodeCfgInput.paaInfo );
Mpc = cell(paraCfgInput.numberOfNodes,...
    max(nPAA_centroids),...
    paraCfgInput.numberOfNodes,...
    max(nPAA_centroids),...
    paraCfgInput.totalNumberOfReflections+1,...
    paraCfgInput.numberOfTimeDivisions+1 );
frmRotMpInfo = cell(1, paraCfgInput.totalNumberOfReflections+1);
keepBothQDOutput =0; % If 1 when using JSON output, it will keep also previous QD output
displayProgress = 1;
ts = paraCfgInput.totalTimeDuration/paraCfgInput.numberOfTimeDivisions;

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
vtx = [0,0,0]; %nodeVelocities(1,:);
vrx = [0,0,0]; %nodeVelocities(2,:);
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

%% Loop over time instances
for iterateTimeDivision = 1:paraCfgInput.numberOfTimeDivisions
    if mod(iterateTimeDivision,100)==0 && displayProgress 
        disp([fprintf('%2.2f', iterateTimeDivision/paraCfgInput.numberOfTimeDivisions*100),'%'])
    end
             
    %% Point rotation: PAAs not centered in the center of the node have a
    % different position in the global frame if the node rotate. Compute
    % the new PAAs position as well as the equivalent angle resulting from
    % successive transformations (initial PAA orientation + rotation of the
    % node over time)
    for nodeId = 1:paraCfgInput.numberOfNodes
            centerRotation = nodePosition(iterateTimeDivision,:, nodeId);
            nodeRotationEucAngles = nodeCfgInput.nodeRotation(iterateTimeDivision,:, nodeId);
            paaInitialPosition = reshape(squeeze(...
                nodeCfgInput.paaInfo{nodeId}.centroidTimePosition(iterateTimeDivision,:,:)), [], 3);
            [paaRotatedPosition, nodeEquivalentRotationAngle] = coordinateRotation(paaInitialPosition, ...
                centerRotation,...
                nodeRotationEucAngles ...
                );
            nodeCfgInput.nodeEquivalentRotationAngle(iterateTimeDivision,:, nodeId) = nodeEquivalentRotationAngle;
            nodeCfgInput.paaInfo{nodeId}.centroid_position_rot(iterateTimeDivision,:,:) =paaRotatedPosition;
    end
    
    % save NodePositionsTrc
    if paraCfgInput.switchSaveVisualizerFiles && ~paraCfgInput.jsonOutput
        filename = sprintf('NodePositionsTrc%d.csv', iterateTimeDivision-1);
        writematrix(squeeze(nodePosition(iterateTimeDivision,:,:)).',fullfile(visualizerPath, filename)...
            );
    end
    
    %% Iterates through all the PAA centroids
    for iterateTx = 1:paraCfgInput.numberOfNodes
        
        for iterateRx = iterateTx+1:paraCfgInput.numberOfNodes
            
            for iteratePaaTx = 1:nPAA_centroids(iterateTx)
                
                for iteratePaaRx = 1:nPAA_centroids(iterateRx)
                    output = [];
                    if (paraCfgInput.numberOfNodes >= 2 || paraCfgInput.switchRandomization == 1)
                        
                        % Update centroids position
                        Tx = squeeze(nodeCfgInput.paaInfo{iterateTx}.centroid_position_rot(iterateTimeDivision,iteratePaaTx,:)).';
                        Rx = squeeze(nodeCfgInput.paaInfo{iterateRx}.centroid_position_rot(iterateTimeDivision,iteratePaaRx,:)).';
                        
                        % Update rotation Tx struct
                        QTx.center(1,:) = nodePosition(iterateTimeDivision,:,iterateTx);
                        QTx.angle(1,:) = nodeCfgInput.nodeRotation(iterateTimeDivision,:, iterateTx);
                        
                        % Update rotation Rx struct
                        QRx.center(1,:) = nodePosition(iterateTimeDivision,:,iterateRx);         
                        QRx.angle(1,:) = nodeCfgInput.nodeRotation(iterateTimeDivision,:, iterateRx);
                        
                        % Update node velocity
                        previousTxPosition =  squeeze(nodeCfgInput.paaInfo{iterateTx}.centroid_position_rot(max(iterateTimeDivision-1,1),iteratePaaTx,:)).';
                        previousRxPosition =  squeeze(nodeCfgInput.paaInfo{iterateRx}.centroid_position_rot(max(iterateTimeDivision-1,1),iteratePaaTx,:)).';

                        vtx = (Tx-previousTxPosition)./ts;
                        vrx = (Rx-previousRxPosition)./ts;
                    end
  
                    % LOS Path generation
                    [switchLOS, output, frmRotMpInfo{1}] = LOSOutputGenerator(CADop, Rx, Tx,...
                        output, vtx, vrx, switchPolarization, switchCp,...
                        polarizationTx, paraCfgInput.carrierFrequency, 'qTx', QTx, 'qRx', QRx);
                    
                    % Store MPC
                    if paraCfgInput.switchSaveVisualizerFiles && switchLOS
                        multipath1 = [Tx, Rx];
                            Mpc{iterateTx,iteratePaaTx,iterateRx,iteratePaaRx, 1, iterateTimeDivision+1} =multipath1;
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
                        
                        Nrealizations = nodeCfgInput.paaInfo{iterateTx}.nodePAAInfo{iteratePaaTx}.indep_stoch_channel*...
                            nodeCfgInput.paaInfo{iterateRx}.nodePAAInfo{iteratePaaRx}.indep_stoch_channel;
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
                        
                        %Store MPC
                        if paraCfgInput.switchSaveVisualizerFiles &&...
                                size(multipathTemporary,1) > 0                            
                            multipath1 = multipathTemporary(1:count,...
                                2:size(multipathTemporary,2));
                                Mpc{iterateTx,iteratePaaTx,iterateRx,iteratePaaRx, iterateOrderOfReflection+1, iterateTimeDivision+1} =multipath1;
                            
                        end
                        
                        %Store QD output
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
    
    %% Generate channel for each PAA given the channel of the centroids
    outputPAATime(:,:,iterateTimeDivision) = generateChannelPaa(outputPAA, nodeCfgInput.paaInfo);  %#ok<AGROW>
    
    %% Write QD output in CSV files
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

%% Write output in JSON files
% QD output
if paraCfgInput.jsonOutput
    writeQdJsonOutput(outputPAATime,cellfun(@(x) x.nPaa,  nodeCfgInput.paaInfo),...
        qdFilesPath);
    
    Mpc(:,:,:,:,:,1) = [];
end

if paraCfgInput.switchSaveVisualizerFiles && paraCfgInput.jsonOutput
    % Write MPC.json
    f = fopen(fullfile(visualizerPath, 'Mpc.json'), 'w');
    for iterateTx = 1:paraCfgInput.numberOfNodes
        for iterateRx = iterateTx+1:paraCfgInput.numberOfNodes
            for iteratePaaTx = 1:nPAA_centroids(iterateTx)
                nodeTxCluster  = nodeCfgInput.paaInfo{iterateTx}.paaInCluster{iteratePaaTx};
                for txPaaCluster = 1:length(nodeTxCluster)
                    for iteratePaaRx = 1:nPAA_centroids(iterateRx)
                        nodeRxCluster  = nodeCfgInput.paaInfo{iterateRx}.paaInCluster{iteratePaaRx};
                        for rxPaaCluster = 1:length(nodeRxCluster)
                            for reflOrd = 1:paraCfgInput.totalNumberOfReflections+1                                
                                s = struct('TX', iterateTx-1, 'PAA_TX', nodeTxCluster(txPaaCluster)-1,...
                                    'RX', iterateRx-1, 'PAA_RX', nodeRxCluster(rxPaaCluster)-1, ...
                                    'Rorder', reflOrd-1);
                                Mpc_t = squeeze((Mpc(iterateTx,iteratePaaTx,...
                                    iterateRx,iteratePaaRx,reflOrd,:)));
                                maxSize = max(cell2mat(cellfun(@size ,Mpc_t ,'UniformOutput' ,false)));
                                    Mpc_t= cellfun(@(x) padNan(x, maxSize),Mpc_t,'UniformOutput',false);                                     
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
    
    % Write Node Position
    filename = sprintf('NodePositions.json');
    f = fopen(fullfile(visualizerPath, filename), 'w');
    for i = 1:paraCfgInput.numberOfNodes
        s = struct('Node' , i-1, ...
            'Position', [nodePosition(:,:,i); [inf inf inf]], ...
            'Rotation', [nodeCfgInput.nodeEquivalentRotationAngle(:,:,i); [inf inf inf]]);
        json = jsonencode(s); % Add a temporary inf vector to make sure
        % more than a single vector will be encoded. Matlab json
        % encoder loses the square brackets when encoding vectors.
        str2remove =',[null,null,null]'; %Temporary string to remove
        rem_ind_start = num2cell(strfind(json, str2remove)); % Find start string to remove
        index2rm = cell2mat(cellfun(@(x) x:x+length(str2remove)-1,rem_ind_start,'UniformOutput',false)); % Create index of char to remove
        json(index2rm) = []; % Remove temporary vector.
        fprintf(f, '%s\n', json);
    end
    fclose(f);
    
    % Write PAAPosition.json
    f = fopen(strcat(visualizerPath, filesep,'PAAPosition.json'), 'w');
    for i = 1:paraCfgInput.numberOfNodes
        idOrientation = 0;
        for paaId = 1:nPAA_centroids(i)
            nodeTxCluster  = nodeCfgInput.paaInfo{i}.paaInCluster{paaId};
            for paaCentroid = 1:length(nodeTxCluster)
                idOrientation = idOrientation+1;
                s = struct('Node', i-1, 'PAA',nodeTxCluster(paaCentroid)-1, ...
                    'Centroid', nodeCfgInput.paaInfo{i}.centroids(paaId)-1, ...
                    'Orientation', nodeCfgInput.nodeAntennaOrientation{i}(idOrientation,:) ,...
                'Position', [reshape(squeeze(nodeCfgInput.paaInfo{i}.centroid_position_rot(:,paaId,:)), [],3); [inf inf inf]]);
                json = jsonencode(s);% Add a temporary inf vector to make sure
                % more than a single vector will be encoded. Matlab json
                % encoder loses the square brackets when encoding vectors.
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

%% Write useful output information. 
writeReportOutput = 0 ; %Set to 0 to allow succeful test.
if writeReportOutput
    f = fopen(strcat(outputPath, filesep,'report.dat'), 'w'); %#ok<UNRCH>
%     elapsedTime = toc;
%     fprintf(f, 'Elapsed Time:\t%f\n', elapsedTime);    
    fprintf(f, 'Device Rotation:\t%d\n', paraCfgInput.isDeviceRotationOn);
%     isInitialOrientationOn = 'true';
    fprintf(f, 'Initial Orientation:\t%d\n', paraCfgInput.isInitialOrientationOn);
%     isPaaCentered = 'true';
    fprintf(f, 'PAA centered:\t%d\n', paraCfgInput.isPaaCentered);
    fclose(f);
end
end

function x = padNan(x, maxSize)
if ~all(size(x) == maxSize)
    sizeX = size(x);
    x(sizeX(1)+1:maxSize(1),1:maxSize(2)) = nan;
end
end