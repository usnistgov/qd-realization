function plotFrame(app)
%PLOTFRAME Plot the frame of the current timestamp


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

% Updated by: Neeraj Varshney <neeraj.varshney@nist.gov> for JSON format, 
% node rotation and paa orientation

plotNodes(app)
plotRays(app)
plotQd(app,'aod')
plotQd(app,'aoa')
end


function plotNodes(app)
delete(app.nodesPlotHandle)
delete(app.paasTxPlotHandle)
delete(app.paasRxPlotHandle)

t = app.currentTimestep;
tx = app.txIndex;
rx = app.rxIndex;
paaPostx = cell2mat(app.timestepInfo(t).paaPos(:,tx));
paaPosrx = cell2mat(app.timestepInfo(t).paaPos(:,rx));
paaOritx = cell2mat(app.timestepInfo(t).paaOri(:,tx));
paaOrirx = cell2mat(app.timestepInfo(t).paaOri(:,rx));
nodeRottx = app.timestepInfo(t).rot(tx,:);
nodeRotrx = app.timestepInfo(t).rot(rx,:);
nodePostx = app.timestepInfo(t).pos(tx,:);
nodePosrx = app.timestepInfo(t).pos(rx,:);

app.paasTxPlotHandle = plotPaa(app.UIAxes,paaPostx,paaOritx,nodePostx,...
    nodeRottx);
app.paasRxPlotHandle = plotPaa(app.UIAxes,paaPosrx,paaOrirx,nodePosrx,...
    nodeRotrx);

nodesPos = app.timestepInfo(t).pos([tx,rx],:);
app.nodesPlotHandle = scatter3(app.UIAxes,...
    nodesPos(:,1), nodesPos(:,2), nodesPos(:,3),50,'s',...
    'r', 'filled');

end

function plotRays(app)
delete(app.raysPlotHandle)

t = app.currentTimestep;
tx = app.txIndex;
rx = app.rxIndex;
paatx = app.numPaas(tx);
paarx = app.numPaas(rx);
for ipaatx = 1:paatx
    for ipaarx = 1:paarx
        
        timestepInfo = app.timestepInfo(t);
        if isempty(timestepInfo.paaInfo(ipaatx,ipaarx).mpcs)
            % no rays
            return
        end
        % if and else loop to deal with discarded combinations in Mpc.json
        if size(timestepInfo.paaInfo(ipaatx,ipaarx).mpcs,1) >= tx
            mpcs = timestepInfo.paaInfo(ipaatx,ipaarx).mpcs(tx,rx,:);
        else
            % use reverse rx/tx 
            mpcs = timestepInfo.paaInfo(ipaarx,ipaatx).mpcs(rx,tx,:);
        end
        if all(cellfun(@isempty, mpcs))
            % use reverse rx/tx
            mpcs = timestepInfo.paaInfo(ipaarx,ipaatx).mpcs(rx,tx,:);
        end

        for i = 1:length(mpcs)
            relfOrder = i - 1;

            coords = mpcs{i};
            [color, width] = getRayAspect(relfOrder);

            app.raysPlotHandle = [app.raysPlotHandle;...
                plot3(app.UIAxes,...
                coords(:,1:3:end)',coords(:,2:3:end)',coords(:,3:3:end)',...
                'Color',color,...
                'LineWidth',width)];
        end
    end
end
end


function plotQd(app,direction)

Tx = app.txIndex;
Rx = app.rxIndex;
t = app.currentTimestep;
paatx = app.numPaas(Tx);
paarx = app.numPaas(Rx);
timestampInfo = app.timestepInfo(t);
for ipaatx = 1:paatx
    for ipaarx = 1:paarx
        if ~isempty(timestampInfo.paaInfo(ipaatx,ipaarx).qdInfo)
            qd = timestampInfo.paaInfo(ipaatx,ipaarx).qdInfo(Tx,Rx);
            raysAppPlotQdStruct(app, qd, direction)
        else
            switch(direction)
                case 'aoa'
                    delete(app.aoaPlotHandle)
                case 'aod'
                    delete(app.aodPlotHandle)
                otherwise
                    error('direction should be either ''aoa'' or ''aod''')
            end
        end
    end
end

end

function paasNodePlotHandle = plotPaa(UIAxes,paalocation,paaorientation,...
    nodeposition,noderotation)

paasNodePlotHandle = zeros(1,size(paalocation,1));
for ipaa = 1:size(paalocation,1)
    left = paalocation(ipaa,2) - 0.2;
    right = paalocation(ipaa,2) + 0.2;
    bottom = paalocation(ipaa,3) - 0.1;
    top = paalocation(ipaa,3) + 0.1;
    y = [left left right right];
    z = [bottom top top bottom];
    x = zeros(size(y)) + paalocation(ipaa,1);
    paasNodePlotHandle(ipaa) = fill3(UIAxes,x, y, z, 'k');
    
    % For PAA Orientation
    if paaorientation(ipaa,1)~=0 % z-axis
        rotate(paasNodePlotHandle(ipaa),[0,0,1],...
            paaorientation(ipaa,1)*180/pi,paalocation(ipaa,:));
        
    end
    if paaorientation(ipaa,2)~=0 % x-axis
        rotate(paasNodePlotHandle(ipaa),[1,0,0],...
            paaorientation(ipaa,2)*180/pi,paalocation(ipaa,:));

    end
    if paaorientation(ipaa,3)~=0 % y-axis
        rotate(paasNodePlotHandle(ipaa),[0,1,0],...
            paaorientation(ipaa,3)*180/pi,paalocation(ipaa,:));

    end
    % For Node Rotation
    if noderotation(1) ~= 0
       rotate(paasNodePlotHandle(ipaa),[0,0,1],...
           noderotation(1)*180/pi,paalocation(ipaa,:));
       
    end
    if noderotation(2) ~= 0
        rotate(paasNodePlotHandle(ipaa),[1,0,0],...
            noderotation(2)*180/pi,paalocation(ipaa,:));
        
    end
    if noderotation(3) ~= 0
        rotate(paasNodePlotHandle(ipaa),[0,1,0],...
            noderotation(3)*180/pi,paalocation(ipaa,:));
        
    end
    

end
end
