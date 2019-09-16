function output = fillOutput(reflOrder, dod, doa, rayLen, pathGain,...
    dopplerFactor, freq)

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