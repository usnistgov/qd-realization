function dir = azel2vec(az, el)
dir = [sind(el).*cosd(az), sind(el).*sind(az), cosd(el)];
end