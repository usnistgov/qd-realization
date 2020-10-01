function b = isPaaPositionsJson(path)
%ISPAAPOSITIONSJSON Function that checks whether the given path matches with
% the expected QD file output file position (as given by the
% documentation). This allows to safely read the file later.
%
% INPUTS:
% - path: file path. This could be either and absolute path, or a relative
% path, starting at least from the Output/ folder.
%
% SEE ALSO: READQDFILE


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
% Copyright 2019-2020 NIST/CTL (steve.blandino@nist.gov)

[remainingPath, fileName, extension] = fileparts(path);
b = strcmp([fileName, extension], 'PAAPosition.json');

[remainingPath, qdFilesFolder] = fileparts(remainingPath);
b = b && strcmp(qdFilesFolder, 'Visualizer');


[~, outputFolder] = fileparts(remainingPath);
b = b && strcmp(outputFolder, 'Output');

end