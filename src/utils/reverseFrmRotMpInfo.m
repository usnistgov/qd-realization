function frmRotMpInfoOut = reverseFrmRotMpInfo(frmRotMpInfo)
%REVERSEFRMROTMPINFO return symmetric link info, eg replace DOA with DOD
%
% INPUTS:
% - frmRotMpInfo
%
% OUTPUTS:
% - frmRotMpInfo
%
% Copyright 2019-2020 NIST/CTL (steve.blandino@nist.gov)

L = length(frmRotMpInfo);

if isempty(frmRotMpInfo)
    frmRotMpInfoOut = frmRotMpInfo;
else

    frmRotMpInfoOut = struct('dod', cell(1,L), ...
        'doa',  cell(1,L), ...
        'TxVel', cell(1,L), ...
        'RxVel', cell(1,L) ...
        );
    [frmRotMpInfoOut.dod] = deal(frmRotMpInfo.doa);
    [frmRotMpInfoOut.doa] = deal(frmRotMpInfo.dod);
    [frmRotMpInfoOut.TxVel] = deal(frmRotMpInfo.RxVel);
    [frmRotMpInfoOut.RxVel] = deal(frmRotMpInfo.TxVel);


end
end