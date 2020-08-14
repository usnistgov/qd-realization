function writeQdJsonOutput(output, paaNodes, qdFilesPath, precision)
%WRITEQDFILEOUTPUT Writes timestamp information to QdFile
%
% INPUTS:
% - output: output matrix formatted as in MULTIPATH and LOSOUTPUTGENERATOR
% - useOptimizedOutputToFile: see PARAMETERCFG
% - fids: see GETQDFILESIDS
% - iTx: index of the TX
% - iRx: index of the RX
% - qdFilesPath: path to Output/Ns3/QdFiles
% - precision: floating point output precision in number of digits
%
% SEE ALSO: GETQDFILESIDS, CLOSEQDFILESIDS, MULTIPATH, LOSOUTPUTGENERATOR, PARAMETERCFG


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

% if ~useOptimizedOutputToFile
%     filename = sprintf('Tx%dRx%d.txt', iTx - 1, iRx - 1);
filepath = fullfile(qdFilesPath, 'qdOutput.json');
fid = fopen(filepath, 'w');
% else
% %     fid = fids(iTx, iRx);
% end
NODES = size(output,1);
ITER  = size(output,3);
floatFormat = sprintf('%%.%dg',precision);
nodeList = 1:NODES;
Noutput = 21;

for tx = nodeList
    for rx = nodeList(nodeList~=tx)
            for txPaa = 1:paaNodes(tx)
                for rxPaa = 1:paaNodes(rx)
                    mimoCh = squeeze(output(tx,rx,:));
                    mimoCh = cellfun(@(x) appendNan(x,Noutput,paaNodes(tx)*paaNodes(rx)), mimoCh, 'UniformOutput', false);
                    sisoCh =cell2mat(cellfun(@(x) x(:,:,(txPaa-1)*paaNodes(rx)+rxPaa), mimoCh,'UniformOutput', false));
                    rowDist = cellfun(@(x) size(x,1), mimoCh);
                    s = struct('TX', tx-1, 'RX', rx-1,...
                        'PAA_TX', txPaa-1, 'PAA_RX', rxPaa-1);
                    s.Delay = mat2cell(single(sisoCh(:,8)), rowDist);
                    s.Gain  = mat2cell(single(sisoCh(:,9)), rowDist);
                    s.Phase = mat2cell(single(sisoCh(:,18)), rowDist);
                    s.AODEL = mat2cell(single(sisoCh(:,11)), rowDist);
                    s.AODAZ = mat2cell(single(sisoCh(:,10)), rowDist);
                    s.AOAEL = mat2cell(single(sisoCh(:,13)), rowDist);
                    s.AOAAZ = mat2cell(single(sisoCh(:,12)), rowDist);
                    json = jsonencode(s);
                    str2remove =',null'; %Temporary string to remove
                    rem_ind_start = num2cell(strfind(json, str2remove)); % Find start string to remove
                    index2rm = cell2mat(cellfun(@(x) x:x+length(str2remove)-1,rem_ind_start,'UniformOutput',false)); % Create index of char to remove
                    json(index2rm) = [];
                    str2remove ='null';
                    rem_ind_start = num2cell(strfind(json, str2remove)); % Find start string to remove
                    index2rm = cell2mat(cellfun(@(x) x:x+length(str2remove)-1,rem_ind_start,'UniformOutput',false)); % Create index of char to remove
                    json(index2rm) = [];
                    fprintf(fid, '%s\n', json);
                end
            end
    end
end
fclose(fid);


end

function x = appendNan(x,n,m)
if isempty(x)
    x = nan(2,n,m);
elseif size(x,3)<m
    x(:, :, size(x,3)+1:m) = nan;
else
    x(end+1,:,:) = nan;
end
end