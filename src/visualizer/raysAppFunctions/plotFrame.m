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
app.paasTxPlotHandle = plotPaa(app.UIAxes,paaPostx);
app.paasRxPlotHandle = plotPaa(app.UIAxes,paaPosrx);

pos = app.timestepInfo(t).pos([tx,rx],:);
app.nodesPlotHandle = scatter3(app.UIAxes,...
    pos(:,1), pos(:,2), pos(:,3),50,'s',...
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

function paasNodePlotHandle = plotPaa(UIAxes,paalocation)
paasNodePlotHandle = zeros(1,size(paalocation,1));
for ipaa = 1:size(paalocation,1)
    left = paalocation(ipaa,1) - 0.25;
    right = paalocation(ipaa,1) + 0.25;
    bottom = paalocation(ipaa,2) - 0.25;
    top = paalocation(ipaa,2) + 0.25;
    x = [left left right right];
    y = [bottom top top bottom];
    z = zeros(size(x)) + paalocation(ipaa,3);
    paasNodePlotHandle(ipaa) = fill3(UIAxes,x, y, z, 'k');
end

end
