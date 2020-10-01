function out = readQdJsonFile(path)
%READQDJSONFILE extracts the QD information from the output JSON file
%
% INPUTS:
% - path: file path. This could be either and absolute path, or a relative
% path, starting at least from the Output/ folder.
%
% OUTPUTS:
% - out: array of JSON struct
% Copyright 2019-2020 NIST/CTL (steve.blandino@nist.gov)

fid = fopen(path);

out = struct('TX', cell(1,1), ...
    'RX', cell(1,1), ...
    'PAA_TX', cell(1,1), ...
    'PAA_RX', cell(1,1), ...
    'Delay', cell(1,1), ...
    'Gain', cell(1,1), ...
    'Phase', cell(1,1), ...
    'AODEL', cell(1,1), ...
    'AODAZ', cell(1,1), ...
    'AOAEL', cell(1,1), ...
    'AOAAZ', cell(1,1) ...
    );

line = 1;
isEof = false;
while ~isEof
    tline = fgetl(fid);
    if ischar( tline )
        out(line) = jsondecode(tline);
        line = line + 1;
    else
        isEof = true;
    end
end
fclose(fid);

end