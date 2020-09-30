function out = readNodeJsonFile(path)

fid = fopen(path);

out = struct('Node', [], ...
    'Position',[], ...
    'Rotation', [] ...
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