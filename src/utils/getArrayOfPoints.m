function arrayOfPoints = getArrayOfPoints(triangIdxs,cadData)
triangles = cadData(triangIdxs,:);
points = triangles(:, 1:9);
arrayOfPoints = reshape(points.', [], 1).';
end
