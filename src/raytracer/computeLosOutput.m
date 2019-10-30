function [output, multipath] = computeLosOutput(rxPos, txPos, rxVel, txVel,...
    cadData, freq, minAbsolutePathGainThreshold)
%COMPUTELOSOUTPUT Computes LoS ray


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

dod = rxPos - txPos;
rayLen = norm(dod);
pathGain = friisPathGain(rayLen, freq);

if pathGain < minAbsolutePathGainThreshold ||...
        isRayObstructed(rxPos, txPos, cadData, [], [])
    
    output = [];
    multipath = [];
    
else
    reflOrder = 0;
    doa = -dod;
    dopplerFactor = getDopplerFactor(txPos, rxPos, txVel, rxVel, [], []);
    
    output = fillOutput(reflOrder, dod, doa, rayLen, pathGain, dopplerFactor, freq);
    multipath = [rxPos, txPos];
end

end