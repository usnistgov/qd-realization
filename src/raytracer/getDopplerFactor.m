function dopplerFactor = getDopplerFactor(txPos, rxPos, txVel, rxVel,...
    cadData, triangIdxList)
%GETDOPPLERFACTOR Get doppler factor from given positions, velocities and
%triangle bounces.


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

txVelRefl = txVel;
txPosRefl = txPos;
for i = 1:length(triangIdxList)
    plane = cadData(triangIdxList(i), 10:13);
    txVelRefl = reflectedVelocity(txVelRefl, plane);
    txPosRefl = reflectedImagePointPlane(txPosRefl, plane);
end

relativeVel = txVelRefl - rxVel;
relativePos = txPosRefl - rxPos;
radialRelativeVel = dot(relativeVel, relativePos) / norm(relativePos);

dopplerFactor = -radialRelativeVel / 3e8;

end