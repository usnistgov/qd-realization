clear
close all
clc

%% Input
scenario = 'L-Room4nodes-fewSamples';
% thresholdType = 'absolute';
% thresholds = [-120, -110, -100, -90, -80];
thresholdType = 'relative';
thresholds = [-50, -40, -30, -25, -20, -10];
qdFilename = 'Tx0Rx1.txt';

%%
qdFilepath = fullfile(scenario, 'Output/Ns3/QdFiles', qdFilename);
simTime = nan(size(thresholds));

for i = 1:length(thresholds)
    % run simulation
    forcedCfgParams = getForcedCfgParams(thresholdType, thresholds(i));
    
    t0 = tic;
    launchRaytracer(scenario, forcedCfgParams);
    simTime(i) = toc(t0);
    
    % prepare plots
    paraCfg = parameterCfg(scenario);
    paraCfg = applyForcedCfgParams(paraCfg, forcedCfgParams);
    name = sprintf('abs %.0f dB, rel %.0f dB',...
        paraCfg.minAbsolutePathGainThreshold,...
        paraCfg.minRelativePathGainThreshold);
    
    data = readQdFile(qdFilepath);
    
    figure(1)
    [y,x] = ecdf([data.pathGain]);
    plot(x,y,'DisplayName', name); hold on

    figure(2)
    pgs = cellfun(@(x) x-max(x), {data.pathGain}, 'UniformOutput', false);
    [y,x] = ecdf([pgs{:}]);
    plot(x,y,'DisplayName', name); hold on

end

figure(1)
legend('show', 'Location', 'northwest')
xlabel('PG [dB]')
ylabel('CDF')
hold off

figure(2)
legend('show', 'Location', 'northwest')
xlabel('$\Delta$PG [dB]')
ylabel('CDF')
hold off

% speedup plot
if thresholdType == "absolute"
    figure(1)
elseif thresholdType == "relative"
    figure(2)
else
    error()
end

yyaxis right
plot(thresholds, simTime(1)./simTime, '-o', 'DisplayName', 'Total speedup')
ylabel('Speedup')


%% UTILS
function forcedCfgParams = getForcedCfgParams(thresholdType, threshold)

forcedCfgParams = struct();

if thresholdType == "absolute"
    forcedCfgParams.minAbsolutePathGainThreshold = threshold;
elseif thresholdType == "relative"
    forcedCfgParams.minRelativePathGainThreshold = threshold;
else
    error('Threshold type ''%s'' not recorgnized', thresholdType)
end

end