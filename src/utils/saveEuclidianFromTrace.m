function saveEuclidianFromTrace(nameTraceFile, nameFile, varargin)
if isempty(varargin)
    captureId = 's1';
else
    captureId = ['s', num2str(varargin{1})];
end
trace = load(nameTraceFile);
traceField = fieldnames(trace);
traceSelect = ~cellfun(@isempty,(regexp(traceField,['.', captureId])));
if sum(traceSelect)
    trace = eval(['trace.',traceField{traceSelect}]);
    angles = trace(:,[4 6 5])/180*pi;
    writematrix(angles, nameFile)
else
    warning('Trace not found')
end
end
