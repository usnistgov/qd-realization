function output = fillOutput(reflOrder, dod, doa, rayLen, pathGain,...
    dopplerFactor, freq)
%FILLOUTPUT Systematically creates a consistent output vector for a given
%ray. Columns are as follows:
% 1. Reflection order
% 2-4. Direction of Departure (AoD, vector pointing from TX to RX, if LoS
% ray, or to first reflection point, if NLoS ray)
% 5-7. Direction of Arrival (AoA, vector pointing from RX to TX, if LoS
% ray, or to last reflection point, if NLoS ray)
% 8. Time delay [s] (ray length at speed of light)
% 9. Path Gain [dB]
% 10-11. AoD azimuth/elevation [deg]
% 12-13. AoA azimuth/elevation [deg]
% 14-17. (TBD)
% 18. Phase shift (caused by reflections, i.e. reflOrder*pi)
% 19. (TBD)
% 20. Doppler frequency
% 21. 0 (for backward compatibility)


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

output = nan(1,21);

% Reflection Order
output(1) = reflOrder;
% Direction of Departure
output(2:4) = dod;
% Direction of Arrival
output(5:7) = doa;
% Time delay
output(8) = rayLen / 3e8;
% Path gain
output(9) = pathGain;
% AoD azimuth
output(10) = mod(atan2d(dod(2),dod(1)), 360);
% AoD elevation
output(11) = acosd(dod(3) / norm(dod));
% AoA azimuth
output(12) = mod(atan2d(doa(2),doa(1)), 360);
% AoA elevation
output(13) = acosd(doa(3) / norm(doa));
% output(14)
% output(15)
% output(16)
% output(17)
% Phase shift caused by reflections
output(18) = reflOrder * pi;
% output(19)
% Doppler frequency
output(20) = dopplerFactor * freq;

output(21) = 0;

end