function [] = createTxRxGrid(mainPath,randomSampling,txPos,visualize)
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

endout=regexp(mainPath,filesep,'split');
cadPath = fullfile(endout{1:2},'cachedCadOutput.mat');
CAD = load(cadPath); % load the CAD file
CADcoord = CAD.CADop; % extract only the coordinates

% extract the (x,y,z) coordinates of the edges of the triangles
triangles = CADcoord(:,1:9);

lRoom = contains(mainPath,'L-Room');
if randomSampling
    paraCfg = parameterCfg(mainPath);
    nSteps = paraCfg.numberOfTimeDivisions;
    [txPosSteps,rxPosSteps] = createGrid(triangles,randomSampling,nSteps,txPos,lRoom);
else
    [txPosSteps,rxPosSteps] = createGrid(triangles,randomSampling,[],txPos,lRoom);
end

% generate also fictious velocity and location files (just for compatibility reasons)
txVel = [0,0,0];
rxVel = [0,0,0];
velocities = [txVel;rxVel];

txInitPos = txPos;
rxInitPos = rxPosSteps(1,:);
initPos = [txInitPos;rxInitPos];

%% save files
% save every grid in a different folder
inputPath = fullfile(mainPath, 'Input');

save([inputPath,'/NodePosition1','.dat'],'txPosSteps','-ascii') % tx position

save([inputPath,'/NodePosition2','.dat'],'rxPosSteps','-ascii') % rx positions

csvwrite(fullfile(inputPath,'nodeVelocities.dat'),velocities) % (fictious) velocities

csvwrite(fullfile(inputPath,'nodes.dat'),initPos) % (fictious) initial positions

%% plots
if visualize
    [faces,edges] = roomCoords2edges(triangles);
    fv = struct();
    fv.faces = faces;
    fv.vertices = edges;
    mask = inpolyhedron(fv,rxPosSteps);
    
    rxPosStepsIn = rxPosSteps(mask,:);
    rxPosStepsOut = rxPosSteps(~mask,:);
    
    [Tri,X,Y,Z] = roomCoords2triangles(triangles);
    figure()
    % scatter3(rxPosSteps(:,1),rxPosSteps(:,2),rxPosSteps(:,3),[],'k')
    % hold on
    scatter3(rxPosStepsIn(:,1),rxPosStepsIn(:,2),rxPosStepsIn(:,3),10,'k','filled')
    hold on
    scatter3(rxPosStepsOut(:,1),rxPosStepsOut(:,2),rxPosStepsOut(:,3),10,'r','filled')
    scatter3(txPosSteps(:,1),txPosSteps(:,2),txPosSteps(:,3),25,'b','filled')
    trisurf(Tri,X,Y,Z,'FaceAlpha',0.2)
end
end
%% utils
function [txPosSteps,rxPosSteps]=createGrid(triangles,randomSampl,nSteps,txPos,lRoom)
x = triangles(:,1:3:end);
y = triangles(:,2:3:end);
z = triangles(:,3:3:end);

%% RX position
minOffset = [0.5, 0.5, 0.5]; % distance of the first sample from the wall
maxOffset = [0.5, 0.5, 0.5]; % distance of the last sample from the wall
resolution = [1, 1, 1]; % step between adjacent samples ([x,y,z])

if lRoom  % split the L shape in two rectangular parts
    % part 1
    maxX(1) = 10;
    maxY(1) = 6;
    maxZ(1) = 3;
    maxs(1,1:3) = [maxX(1), maxY(1), maxZ(1)];
    
    minX(1) = 0;
    minY(1) = 0;
    minZ(1) = 0;
    mins(1,1:3) = [minX(1), minY(1), minZ(1)];
    
    % part 2
    maxX(2) = 10;
    maxY(2) = 19;
    maxZ(2) = 3;
    maxs(2,1:3) = [maxX(2), maxY(2), maxZ(2)];
    
    minX(2) = 6;
    minY(2) = 6-2*minOffset(1); % remove the offset between the two parts
    minZ(2) = 0;
    mins(2,1:3) = [minX(2), minY(2), minZ(2)];
else
    [maxX, maxXInd] = max(max(x));
    [maxY, maxYInd] = max(max(y));
    [maxZ, maxZInd] = max(max(z));
    maxs = [maxX, maxY, maxZ];
    
    [minX, minXInd] = min(min(x));
    [minY, minYInd] = min(min(y));
    [minZ, minZInd] = min(min(z));
    mins = [minX, minY, minZ];
end

% sampling
grid = [];
for j=1:lRoom+1
    if ~randomSampl
        
        x = mins(j,1) + minOffset(1):resolution(1):maxs(j,1) - maxOffset(1);
        y = mins(j,2) + minOffset(2):resolution(2):maxs(j,2) - maxOffset(2);
        z = mins(j,3) + minOffset(3):resolution(3):maxs(j,3) - maxOffset(3);
        [xRange, yRange, zRange] = meshgrid(x, y, z); % create the grid
        
        xRange = xRange(:);
        yRange = yRange(:);
        zRange = zRange(:);
        
        grid = [grid;[xRange,yRange,zRange]];
    else
        % TODO: don't split samples in half
        gridTmp = mins(j,:) + minOffset + (maxs(j,:)-maxOffset-(mins(j,:)+minOffset)).*rand(nSteps/(lRoom+1),3);
        grid = [grid;gridTmp];
    end
end
if ~randomSampl
    nSteps = size(grid,1);
end
%% each TX-RX position is tested in a different time step

% position grids
txPosSteps = repmat(txPos,nSteps,1); % the transmitter position is fixed
rxPosSteps = grid;

end

function [Tri,edges] = roomCoords2edges(roomCoords)
nTriang = size(roomCoords,1);

X=roomCoords(:,1:3:end)';
xEdges = X(:);

Y=roomCoords(:,2:3:end)';
yEdges = Y(:);

Z=roomCoords(:,3:3:end)';
zEdges = Z(:);

edges = unique([xEdges,yEdges,zEdges],'rows');
nEdges = size(edges,1);

Tri = nan(nEdges,3);
for i=1:nTriang
    edge1 = find(all((roomCoords(i,1:3)==edges),2));
    edge2 = find(all((roomCoords(i,4:6)==edges),2));
    edge3 = find(all((roomCoords(i,7:9)==edges),2));
    Tri(i,:) = [edge1,edge2,edge3];
end

end

function [Tri,X,Y,Z] = roomCoords2triangles(roomCoords)
nTriangles = size(roomCoords,1);

X=roomCoords(:,1:3:end)';
Y=roomCoords(:,2:3:end)';
Z=roomCoords(:,3:3:end)';

Tri = reshape(1:nTriangles*3, 3, [])';

end
