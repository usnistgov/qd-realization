function [az, el] = vectorToAngle(angleVector)

az = mod(atan2d(angleVector(:,2),angleVector(:,1)), 360);

el = acosd(angleVector(:,3) ./ vecnorm(angleVector,2,2));
end