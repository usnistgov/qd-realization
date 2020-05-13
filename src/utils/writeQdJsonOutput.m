function writeQdJsonOutput(output, nPAA_centroids, qdFilesPath, precision)
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


for tx = 1:NODES
    for rx = tx+1:NODES
        for txPaa = 1:nPAA_centroids(tx)
            for rxPaa = 1:nPAA_centroids(rx)
                mimoCh = squeeze(output(tx,rx,:));
                sisoCh =cell2mat(cellfun(@(x) x(:,:,(txPaa-1)*rxPaa+rxPaa), mimoCh,'UniformOutput', false));
                rowDist = cellfun(@(x) size(x,1), mimoCh);
                s = struct('TX', tx-1, 'RX', rx-1,...
                    'PAA_TX', txPaa-1, 'PAA_RX', rxPaa-1, ...
                    'Delay',mat2cell(single(sisoCh(:,8)), rowDist),...
                    'Gain', mat2cell(single(sisoCh(:,9)), rowDist),...
                    'Phase',mat2cell(single(sisoCh(:,18)), rowDist),...
                    'AODEL',mat2cell(single(sisoCh(:,11)), rowDist),...
                    'AODAZ',mat2cell(single(sisoCh(:,10)), rowDist),...
                    'AOAEL',mat2cell(single(sisoCh(:,13)), rowDist),...
                    'AOAAZ',mat2cell(single(sisoCh(:,12)), rowDist)...
                    );
                json = jsonencode(s);
                fprintf(fid, '%s\n', json);
            end
        end
    end
end
fclose(fid);


end