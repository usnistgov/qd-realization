function outputPath = launchRaytracer(varargin)
%LAUNCHRAYTRACER Launch raytracer with given parameters. Acts as a
%functionalized main script
%
%LAUNCHRAYTRACER(): runs default scenario ('ScenarioTest')
%LAUNCHRAYTRACER(scenarioName): runs the given scenario. scenarioName can
%contain a path relative to this function's parent folder, e.g.,
%examples/Indoor1. Default: ''.
%LAUNCHRAYTRACER(scenarioName, forcedParaCfg): additionally, forcedParaCfg
%overwrites the configuration parameters obtained by the configuration
%input file. Default: struct().
%LAUNCHRAYTRACER(__, 'verbose', v): 'verbose' parameter with level v. If
%v==0, no information is written on the command window. Default: v = 1.

% Input handling
p = inputParser;

addOptional(p, 'scenarioName', '', @(x) isStringScalar(x) || ischar(x));
addOptional(p, 'forcedParaCfg', struct(), @isstruct);
addParameter(p, 'verbose', 1, @(x) validateattributes(x, ...
    {'numeric'}, {'scalar', 'nonempty', 'integer', 'nonnegative'},...
    mfilename, 'verbose'));

parse(p, varargin{:});

scenarioName = p.Results.scenarioName;
forcedParaCfg = p.Results.forcedParaCfg;
verbose = p.Results.verbose;

% Init
functionPath = fileparts(mfilename('fullpath'));
addpath(fullfile(functionPath, 'raytracer'),...
    fullfile(functionPath, 'utils'))

% Input
if ~isempty(scenarioName)
    if verbose > 0
        fprintf('Use customized scenario: %s.\n',scenarioName);
    end
else
    scenarioName = 'ScenarioTest';
    
    if verbose > 0
        fprintf('Use default scenario: ScenarioTest.\n');
    end
end
scenarioPathStr = fullfile(functionPath,scenarioName);

% Check input scenario file
if ~isfolder(scenarioPathStr)
    scenarioInputPath = fullfile(functionPath, scenarioName, 'Input');
    mkdir(scenarioInputPath);
    
    copyfile(fullfile(functionPath, 'Input'), scenarioInputPath);
    
    if verbose > 0
        fprintf(['%s folder does not exist, creating a new folder with',...
            ' default scenario from root Input folder.\n'],scenarioName);
    end
    
else
    if verbose > 0
        fprintf('%s folder already exists and using this scenario to process.\n',...
            scenarioName);
    end
    
end

% Input system and node-related parameters
paraCfg = parameterCfg(scenarioName);
[paraCfg, nodeCfg] = nodeProfileCfg(paraCfg);

% Apply forced configuration parameters
paraCfg = applyForcedCfgParams(paraCfg, forcedParaCfg);

% Run raytracing function and generate outputs
outputPath = Raytracer(paraCfg, nodeCfg);

end