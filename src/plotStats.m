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

folderList = ["Indoor1",fullfile("Camillo_NIST_60-GHz Lecture Room--MPCs_1551899818","data","TX1_measurements")];

% select plot metrics
features = ["pathGain","delay","aoaAz"]; % specify the metric to observe for the statistics
if strcmp(features,"all") % if "all" all the stats will be plotted
    features = fieldnames(stats);
end

% whether to save the figures or not
savefigs = 1;
if savefigs
    endout=regexp(statsFilePath,filesep,'split');
    statsFigFolder = endout{1:end-1}; % path where to save the plots
end

for folder = folderList
    
    statsFilePath = fullfile(folder,"Stats","stats.mat");
    load(statsFilePath)
    
    for featID = 1:length(features)
        feat = features{featID};
        
        if ~isfield(stats,feat)
            error('feature must be a field of Tx%sRx%s.txt',txId,rxId);
        end
        currFeature = extractfield(stats,feat);
        
        figure(featID)
        subplot(1,2,1)
        ecdf(currFeature)
        xlabel(sprintf('$x=$ %s',feat))
        ylabel('$F(x)$')
        % title('Path Gain distribution in a rectangular room. $N_{refl}=2$')
        hold on
        
        subplot(1,2,2)
        histogram(currFeature,'Normalization','pdf')
        xlabel(feat)
        ylabel('Count')
        % title('Path Gain distribution in a rectangular room. $N_{refl}=2$')
        hold on
        
        fileName = fullfile(statsFigFolder,sprintf('%s.fig',feat));
        savefig(fileName)
    end
end