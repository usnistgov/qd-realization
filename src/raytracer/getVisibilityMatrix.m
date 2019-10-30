function visibilityMatrix = getVisibilityMatrix(cadOutput)
nTriangles = size(cadOutput,1);

rows = [];
cols = [];
values = [];

for i = 1:nTriangles
    t1 = cadOutput(i,:);
    for j = i+1:nTriangles
        t2 = cadOutput(j,:);
        
        areVisible = isTriangleInFrontOfTriangle(t1, t2, true) &&...
            isTriangleInFrontOfTriangle(t2, t1, true);
        
        if areVisible
            rows = [rows, i, j];
            cols = [cols j, i];
            values = [values, true, true];
        end
        
        
    end
end

visibilityMatrix = sparse(rows,cols,values);

end