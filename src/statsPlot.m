scenario = 'L-Room4nodes-fewSamples';
qdFiles = fullfile(scenario,'Output/Ns3/QdFiles/Tx0Rx1.txt');

data = readQdFile(qdFiles);

figure(1)
[y,x] = ecdf([data.pathGain]);
plot(x,y,'DisplayName','baseline')

legend('show','Location','northwest')
xlabel('PG [dB]')
ylabel('CDF')
hold on

figure(2)
pgs = cellfun(@(x) x-max(x), {data.pathGain}, 'UniformOutput', false);
[y,x] = ecdf([pgs{:}]);
plot(x,y,'DisplayName','baseline')

legend('show','Location','northwest')
xlabel('$\Delta$PG [dB]')
ylabel('CDF')
hold on

%%
scenario = 'L-Room4nodes-fewSamples-rel10';
qdFiles = fullfile(scenario,'Output/Ns3/QdFiles/Tx0Rx1.txt');

data = readQdFile(qdFiles);
paraCfg = parameterCfg(scenario);
name = sprintf('abs %.0f dB, rel %.0f dB',...
    paraCfg.minAbsolutePathGainThreshold,...
    paraCfg.minRelativePathGainThreshold);

figure(3)
[y,x] = ecdf([data.pathGain]);
plot(x,y,'DisplayName',name)

figure(4)
pgs = cellfun(@(x) x-max(x), {data.pathGain}, 'UniformOutput', false);
[y,x] = ecdf([pgs{:}]);
plot(x,y,'DisplayName',name)

%%
PG = [-120, -110, -100, -90, -80];
absThr = [24.887, 22.744;...
    19.922, 18.045;...
    14.524, 13.044;...
    13.019, 11.703;...
    12.378, 11.420];

deltaPG = [-50, -40, -30, -25, -20, -10];
relThr = [23.976, 21.945;...
    24.063, 22.086;...
    22.943, 21.057;...
    21.370, 19.561;...
    18.752, 17.099;...
    16.293, 14.817];

figure(1)
yyaxis right
% plot(PG, absThr(:,1), '-o', 'DisplayName', 'Tot sim. time')
% plot(PG, absThr(:,2), '--o', 'DisplayName', 'Multipath time')
% ylabel('Time [s]')
plot(PG, absThr(1,1)./absThr(:,1), '-o', 'DisplayName', 'Total speedup')
ylabel('Speedup')

figure(4)
yyaxis right
% plot(deltaPG, relThr(:,1), '-o', 'DisplayName', 'Tot sim. time')
% plot(deltaPG, relThr(:,2), '--o', 'DisplayName', 'Multipath time')
% ylabel('Time [s]')
plot(deltaPG, relThr(1,1)./relThr(:,1), '-o', 'DisplayName', 'Total speedup')
ylabel('Speedup')