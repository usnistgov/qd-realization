function  [] =  savePositionFromTrace(nameTraceFile,nameFile, varargin)
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
    position = trace(:,[1 3 2])*2.54/100;
    writematrix(position, nameFile)
else
    warning('Trace not found')
end
end
