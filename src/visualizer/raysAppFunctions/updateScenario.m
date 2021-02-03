function updateScenario(app, mainPath)
%UPDATESCENARIO Prepare visualization of new scenario. Plots new
%environment, read and load Output/ files.
%
%SEE ALSO: RAYSAPPPLOTROOM, SETUPTXSLIDER, SETUPRXSLIDER, UPDATETIMESTEP


% Copyright (c) 2020, University of Padova, Department of Information
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

if nargin < 2
    mainPath = uigetdir(app.srcPath);
end

if mainPath == 0
    return
end

scenarioName = getScenarioNameFromPath(mainPath, true);
if strcmp(scenarioName, '')
    uialert(app.Visualizer,...
        'Selected path is not a valid scenario to be visualized',...
        'Invalid path')
    return
else
    app.scenarioName = scenarioName;
end

if strcmp(app.UIAxes.Title.Interpreter, 'latex')
    app.UIAxes.Title.String = strrep(scenarioName, '_', '\_');
else
    app.UIAxes.Title.String = scenarioName;
end

app.outputPath = fullfile(mainPath,'Output');
app.visualizerPath = fullfile(app.outputPath,'Visualizer');
app.ns3Path = fullfile(app.outputPath,'Ns3');

raysAppPlotRoom(app)

setupNodes(app);
setupPaas(app);

setupTimestepInfo(app);

end


%% Utils
function setupNodes(app)
initialPos = readNodeJsonFile(sprintf('%s/NodePositions.json',...
    app.visualizerPath));

app.numNodes = size([initialPos.Node],2);

if app.numNodes < 2
    error('There should be at least 2 nodes in the scenario')
end

% else
setupTxRxSliders(app);
setupTxSlider(app,1);
setupRxSlider(app,2);

end


function setupTxRxSliders(app)
app.TxDropdown.Items = array2cellstr(0:app.numNodes-1);
app.RxDropdown.Items = array2cellstr(0:app.numNodes-1);
end


function setupTimestepInfo(app)
app.timestepInfo = struct();

extractQdFilesInfo(app);
extractMpcCoordinatesInfo(app);
extractNodePositionsInfo(app);
extractPaaInfo(app);

totalTimesteps = length(app.timestepInfo);
app.currentTimestep = 1; % added listener to this variable 

if totalTimesteps < 2
    app.TimestepSlider.Limits = [0,1];
    app.TimestepSpinner.Limits = [0,1];
    
    app.TimestepSlider.Enable = 'off';
    app.TimestepSpinner.Enable = 'off';
    app.PlayButton.Enable = 'off';
    app.PauseButton.Enable = 'off';
else
    app.TimestepSlider.Limits = [1, totalTimesteps];
    app.TimestepSpinner.Limits = [1, totalTimesteps];
    
    app.TimestepSlider.Enable = 'on';
    app.TimestepSpinner.Enable = 'on';
    app.PlayButton.Enable = 'on';
    app.PauseButton.Enable = 'on';
end
end


function extractQdFilesInfo(app)
qdFiles = dir(sprintf('%s/QdFiles',app.ns3Path));
qdfileName = qdFiles(3).name;
[~,~,extension]=fileparts(qdfileName);
switch extension
    case '.json'
        qdFile = readQdJsonFile(sprintf('%s/QdFiles/qdOutput.json',...
            app.ns3Path));
        for i = 1:size(qdFile,2)
            for timestep = 1:size(qdFile(1).Delay,1)
                out = struct('TX', cell(1,1), ...
                    'RX', cell(1,1), ...
                    'PAA_TX', cell(1,1), ...
                    'PAA_RX', cell(1,1), ...
                    'Delay', cell(1,1), ...
                    'Gain', cell(1,1), ...
                    'Phase', cell(1,1), ...
                    'AODEL', cell(1,1), ...
                    'AODAZ', cell(1,1), ...
                    'AOAEL', cell(1,1), ...
                    'AOAAZ', cell(1,1) ...
                    );
                out.TX = qdFile(i).TX + 1;
                out.RX = qdFile(i).RX + 1;
                out.PAA_TX = qdFile(i).PAA_TX + 1;
                out.PAA_RX = qdFile(i).PAA_RX + 1;
                if size(qdFile(1).Delay,1)> 1
                    out.Delay =  cell2mat(qdFile(i).Delay(timestep, :));
                    out.Gain =   cell2mat(qdFile(i).Gain(timestep, :));
                    out.Phase =  cell2mat(qdFile(i).Phase(timestep, :));
                    out.AODEL =  cell2mat(qdFile(i).AODEL(timestep, :));
                    out.AODAZ =  cell2mat(qdFile(i).AODAZ(timestep, :));
                    out.AOAEL =  cell2mat(qdFile(i).AOAEL(timestep, :));
                    out.AOAAZ =  cell2mat(qdFile(i).AOAAZ(timestep, :));
                else
                    out.Delay =  (qdFile(i).Delay(timestep, :));
                    out.Gain =   (qdFile(i).Gain(timestep, :));
                    out.Phase =  (qdFile(i).Phase(timestep, :));
                    out.AODEL =  (qdFile(i).AODEL(timestep, :));
                    out.AODAZ =  (qdFile(i).AODAZ(timestep, :));
                    out.AOAEL =  (qdFile(i).AOAEL(timestep, :));
                    out.AOAAZ =  (qdFile(i).AOAAZ(timestep, :));
                end
                
                app.timestepInfo(timestep).paaInfo(out.PAA_TX,out.PAA_RX).qdInfo(out.TX,out.RX) = out;
            end
        end
    case '.txt'
        for i = 1:length(qdFiles)
            token = regexp(qdFiles(i).name,'Tx(\d+)Rx(\d+).txt','tokens');
            if isempty(token)
                continue
            end
            
            % else
            Tx = str2double(token{1}{1}) + 1;
            Rx = str2double(token{1}{2}) + 1;
            
            qd = readQdFile(sprintf('%s/%s',...
                qdFiles(i).folder,qdFiles(i).name));
            
            for t = 1:length(qd)
                app.timestepInfo(t).qdInfo(Tx,Rx) = qd(t);
            end
        end
end

end


function extractMpcCoordinatesInfo(app)
mpcFile = readMpcJsonFile(sprintf('%s/Mpc.json',...
    app.visualizerPath));
TxList = [mpcFile.TX];
RxList = [mpcFile.RX];
PaaTxList = [mpcFile.PAA_TX];
PaaRxList = [mpcFile.PAA_RX];
Rorder = [mpcFile.Rorder];
for i = 1:size(mpcFile,2)

    Tx = TxList(i) + 1;
    Rx = RxList(i) + 1;
    PaaTx = PaaTxList(i) + 1;
    PaaRx = PaaRxList(i) + 1;    
    Refl = Rorder(i) + 1;
    
    for timestep = 1:size(mpcFile(i).MPC,1)
        if ~iscell(mpcFile(i).MPC)
            switch Refl 
                case 1
                    numColumns = 6;
                case 2
                    numColumns = 9;
                case 3 
                    numColumns = 12;
                case 4
                    numColumns = 15;
            end
            app.timestepInfo(timestep).paaInfo(PaaTx,PaaRx).mpcs{Tx,Rx,Refl} ...
                        = reshape(mpcFile(i).MPC(timestep,:,:),[],numColumns);
         else
             app.timestepInfo(timestep).PaaInfo(PaaTx,PaaRx).mpcs{Tx,Rx,Refl} = [];
        end

    end
end

end


function extractNodePositionsInfo(app)
posFile = readNodeJsonFile(sprintf('%s/NodePositions.json',app.visualizerPath));

for timestep = 1:size(posFile(1).Position,1) 
    pos = zeros(size(posFile,2),3);
    rot = zeros(size(posFile,2),3);
    for i = 1:size(posFile,2)
        pos(i,:) = posFile(i).Position(timestep,:,:);
        rot(i,:) = posFile(i).Rotation(timestep,:,:);

    end
    app.timestepInfo(timestep).pos = pos;
    app.timestepInfo(timestep).rot = rot;

end

end

function setupPaas(app)
% get number of PAAs per node e.g., [2;2]- 2 nodes each has 2 PAAs
paaPosFile = readPaaJsonFile(sprintf('%s/PAAPosition.json',app.visualizerPath));

getPaaInfo = tabulate([paaPosFile.Node]);

app.numPaas = getPaaInfo(:,2);

if app.numPaas < 1
    error('There should be at least 1 paa per node in the scenario')
end

end

function extractPaaInfo(app)
paaPosFile = readPaaJsonFile(sprintf('%s/PAAPosition.json',app.visualizerPath));
getPaaInfo = tabulate([paaPosFile.Node]);
for timestep = 1:size(paaPosFile(1).Position,1) 
    index  = 1;
    paaPos = {};
    paaOri = {};
    for iNode = 1:size(getPaaInfo,1)       
        for iPaa = 1:getPaaInfo(iNode,2)           
            paaPos{iPaa,iNode} = paaPosFile(index).Position(timestep,:,:);
            paaOri{iPaa,iNode} = reshape(paaPosFile(index).Orientation,1,[]);

            index = index+1;
        end        
    end
    app.timestepInfo(timestep).paaPos = paaPos;
    app.timestepInfo(timestep).paaOri = paaOri; 

end

end