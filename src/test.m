N = 1000;

a = rand(N,14);

b = table(a(:,1:3), a(:,4:6), a(:,7:9), a(:,10:13), a(:,14),...
    'VariableNames',{'p1','p2','p3','plane','material'});
c = table2struct(b);


for i = 1:1e6
    idx = randi(N);
    
    p1a = a(idx, 1:3);
    p2a = a(idx, 4:6);
    p3a = a(idx, 7:9);
    plane_a = a(idx, 10:13);
    material_a = a(idx, 14);
    
%     p1b = b.p1(idx, :);
%     p2b = b.p2(idx, :);
%     p3b = b.p3(idx, :);
%     plane_b = b.plane(idx, :);
%     material_b = b.material(idx);
    
    p1c = c(idx).p1;
    p2c = c(idx).p2;
    p3c = c(idx).p3;
    plane_c = c(idx).plane;
    material_c = c(idx).material;
    
end
% 3233, 1663 -> 0.5144
% 1594, 3139 -> 0.5078