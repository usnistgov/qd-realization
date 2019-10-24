function [] = genParameterCfg(cfgParams,inputPath)

endout=regexp(inputPath,filesep,'split');
scenarioNameStr = endout{2};

% xml room file
% = "Box.xml" (Default file name)
if ~isfield(cfgParams,'environmentFileName')
    xmlPath = fullfile(endout{1},scenarioNameStr);
    dirList = dir(xmlPath);
    fileList = {dirList.name};
    xmlFiles = endsWith(fileList,'.xml');
    if sum(xmlFiles)>1
        error('Too many .xml files in %s.',inputPath)
    elseif sum(xmlFiles)==0
        error('There is no .xml file in %s.',inputPath)
    else
        cfgParams.environmentFileName=string(fileList{xmlFiles}); % convert to string for saving via table
    end
end

% Generalized Scenario
% = 1 (Default)
if ~isfield(cfgParams,'generalizedScenario')
    cfgParams.generalizedScenario=1;
end

% Switch Indoor
% = 1;
if ~isfield(cfgParams,'indoorSwitch')
    cfgParams.indoorSwitch=1;
end

% This is switch to turn on or off mobility.
% 1 = mobility ON, 0 = mobility OFF (Default)
% TODO: can be made smart checking if file exists
if ~isfield(cfgParams,'mobilitySwitch')
    defaultMobilitySwitch = getDefaultMobilitySwitch(scenarioNameStr);
    cfgParams.mobilitySwitch= defaultMobilitySwitch;
end

% This switch lets the user to decide the input to mobility
% 1 = Linear (Default), 2 = input from File
% TODO: make it smart: is a valid mobility file is present, set 2
if ~isfield(cfgParams,'mobilityType')
    defaultMobilityType = getDefaultMobilityType(scenarioNameStr, cfgParams.mobilitySwitch);
    cfgParams.mobilityType= defaultMobilityType;
end

% n is the total number of time divisions. If n  = 100 and t  = 10, then we
% have 100 time divisions for 10 seconds. Each time division is 0.1 secs in
% length
% = 10 (Default)
if ~isfield(cfgParams,'numberOfTimeDivisions')
    error('you have to specify the number of time divisions');
end

% Reference point is the center of limiting sphere. 
% = [3,3,2] (Default)
% TODO: default value is arbitrary and not generic at all
if ~isfield(cfgParams,'referencePoint')
    cfgParams.referencePoint = "[3,3,2]"; % TODO: save it as numeric
end

% This is selection of planes/nodes by distance. r = 0 means that there is
% no limitation (Default). 
if ~isfield(cfgParams,'selectPlanesByDist')
    cfgParams.selectPlanesByDist = 0;
end

% Switch to turn ON or OFF the Qausi dterministic module
% 1 = ON, 0 = OFF (Default)
if ~isfield(cfgParams,'switchQDGenerator')
    cfgParams.switchQDGenerator = 0;
end

% This is switch to turn ON or OFF randomization.
% 1 = random (Default), 0 = Tx,Rx are determined by Tx,Rx paramters
if ~isfield(cfgParams,'switchRandomization')
    cfgParams.switchRandomization = 0;
end

% This parameter denotes the number of nodes
% = 2  (Default)
if ~isfield(cfgParams,'numberOfNodes')
    switch(cfgParams.switchRandomization)
        case 0
            defaultNumberOfNodes = ""; % specified by nodes.dat
            % TODO: save as defaultNumberOfNodes=[]
        case 1
            defaultNumberOfNodes = 2;
        otherwise
            error('Cannot handle switchRandomization=%f',cfgParams.switchRandomization)
    end
    cfgParams.numberOfNodes = defaultNumberOfNodes;
end

% Switch to enable or disable the visuals
% = 0 (Default)
if ~isfield(cfgParams,'switchVisuals')
    cfgParams.switchVisuals = 0;
end

% Order of reflection.
% 1 = multipath until first order, 2 = multipath until second order (Default)
if ~isfield(cfgParams,'totalNumberOfReflections')
    cfgParams.totalNumberOfReflections = 2;
end

% t is the time period in seconds. The time period for which the simulation
% has to run when mobility is ON
% = 1 (Default)
if ~isfield(cfgParams,'totalTimeDuration')
    cfgParams.totalTimeDuration = 1;
end

cfgParamsTab = sortrows(rows2vars(struct2table(cfgParams)));
writetable(cfgParamsTab,fullfile(inputPath,'paraCfgCurrent.txt'),'FileType','text',...
    'Delimiter','tab')  
end
function defaultMobilityType = getDefaultMobilityType(scenarioNameStr, mobilitySwitch)
if ~mobilitySwitch
    defaultMobilityType = 2;
    return
end

% Mobility switch activated
if isNodePositionPresent(scenarioNameStr)
    defaultMobilityType = 1;
else
    defaultMobilityType = 0;
end

end
function defaultMobilitySwitch = getDefaultMobilitySwitch(scenarioNameStr)
if isNodePositionPresent(scenarioNameStr)
    defaultMobilitySwitch = 1;
else
    defaultMobilitySwitch = 0;
end

end
function b = isNodePositionPresent(path)
files = dir(sprintf('%s/Input',path));

b = any(startsWith({files.name},'NodePosition'));
end