function arrayOfPlanes = getArrayOfPlanes(triangIdxs,cadData)
reflOrder = size(triangIdxs,2);
triangles = cadData(triangIdxs,:);
planes = triangles(:, 10:13);
arrayOfPlanes = reshape(planes.', [], 1).';
arrayOfPlanes = [reflOrder, arrayOfPlanes];
end
