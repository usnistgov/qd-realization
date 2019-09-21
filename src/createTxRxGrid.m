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
% Modified by: Mattia Lecci <leccimat@dei.unipd.it>, Automatic path check

clear
close all
clc

addpath('raytracer', 'utils')

scenario = 'ScenarioTest';
CADPath = fullfile(scenario, 'Input', 'cachedCadOutput');

CAD = load(CADPath); % load the CAD file
CADcoord = CAD.CADop; % extract only the coordinates

% extract the (x,y,z) coordinates of the edges of the triangles
triangles = CADcoord(:,1:9);
x = CADcoord(:,[1,4,7]);
y = CADcoord(:,[2,5,8]);
z = CADcoord(:,[3,6,9]);

[maxX, maxXInd] = max(x);
[maxY, maxYInd] = max(y);
[maxZ, maxZInd] = max(z);
maxs = [maxX, maxY, maxZ];

[minX, minXInd] = min(x);
[minY, minYInd] = min(y);
[minZ, minZInd] = min(z);
mins = [minX, minY, minZ];

%% RX position
% TODO: only works for parallelepipeds
minOffset = [0.5, 0.5, 0.5]; % distance of the first sample from the wall
maxOffset = [0.5, 0.5, 0.5]; % distance of the last sample from the wall
resolution = [1, 1, 1]; % step between adjacent samples ([x,y,z])

xRange = minX + minOffset(1):resolution(1):maxX - maxOffset(1); 
yRange = minY + minOffset(2):resolution(2):maxY - maxOffset(2);
zRange = 1.5;%minZ + minOffset(3):resolution(3):maxZ - maxOffset(3);

[xGrid, yGrid, zGrid] = meshgrid(xRange, yRange, zRange); % create the grid

coordGrid = cat(length(size(zGrid))+1,xGrid,yGrid,zGrid);
%% TX position
txPos = [8,8,2.5]; % the transmitter position is fixed

%% each TX-RX position is tested in a different time step
nSteps = numel(xGrid);

txPosSteps = repmat(txPos,nSteps,1);
rxPosSteps = [xGrid(:),yGrid(:),zGrid(:)];

savePath = fullfile(scenario, 'Input'); % tx position
save([savePath,'/NodePosition1','.dat'],'txPosSteps','-ascii')

savePath = fullfile(scenario, 'Input'); % rx position
save([savePath,'/NodePosition2','.dat'],'rxPosSteps','-ascii')

