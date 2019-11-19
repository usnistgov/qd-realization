function outputPath = launchRaytracer(varargin)
%LAUNCHRAYTRACER Launch raytracer with given parameters. Acts as a
%functionalized main script
%
%LAUNCHRAYTRACER(): runs default scenario ('ScenarioTest')
%LAUNCHRAYTRACER(scenarioName): runs the given scenario. scenarioName can
%contain a path relative to this function's parent folder, e.g.,
%examples/Indoor1
%LAUNCHRAYTRACER(scenarioName, forcedParaCfg): additionally, forcedParaCfg
%overwrites the configuration parameters obtained by the configuration
%input file.

% Input handling
p = inputParser;

addOptional(p, 'scenarioName', '', @(x) isStringScalar(x) || ischar(x));
addOptional(p, 'forcedParaCfg', struct(), @isstruct);

parse(p, varargin{:});

scenarioName = p.Results.scenarioName;
forcedParaCfg = p.Results.forcedParaCfg;

% Init
functionPath = fileparts(mfilename('fullpath'));
addpath(fullfile(functionPath, 'raytracer'),...
    fullfile(functionPath, 'utils'))

% Input
if ~isempty(scenarioName)
    fprintf('Use customized scenario: %s.\n',scenarioName);
else
    scenarioName = 'ScenarioTest';
    fprintf('Use default scenario: ScenarioTest.\n');
end
scenarioPathStr = fullfile(functionPath,scenarioName);

% Check input scenario file
if ~isfolder(scenarioPathStr)
    scenarioInputPath = fullfile(functionPath, scenarioName, 'Input');
    mkdir(scenarioInputPath);
    
    copyfile(fullfile(functionPath, 'Input'), scenarioInputPath);
    
    fprintf(['%s folder does not exist, creating a new folder with',...
        ' default scenario from root Input folder.\n'],scenarioName);
    
else
    fprintf('%s folder already exists and using this scenario to process.\n',...
        scenarioName);
    
end

% Input system and node-related parameters
paraCfg = parameterCfg(scenarioName);
[paraCfg, nodeCfg] = nodeProfileCfg(paraCfg);

% Apply forced configuration parameters
paraCfg = applyForcedCfgParams(paraCfg, forcedParaCfg);

% Run raytracing function and generate outputs
outputPath = Raytracer(paraCfg, nodeCfg);

end