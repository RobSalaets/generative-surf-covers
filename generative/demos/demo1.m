clear
close all

names_cell = dir('generative/meshes/all/*.ply'); names_cell = {names_cell.name};
labels_cell = dir('generative/meshes/all/*.txt'); labels_cell = {labels_cell.name};
for i = 1:3
    fname = names_cell{i};
    flabel = labels_cell{i};
    flID = fopen(strcat('meshes/all/',flabel),'r');
    L = fscanf(flID,'%i\n');
    fclose(flID);
    [V,F] = read_ply(strcat('meshes/all/',fname));
    V(:,3) = -V(:,3);
    m = mean(V, 1);
    V = V-m;
    s = max(abs(V), [], 'all');
    V = (V)./s;
    
    figure
    trisurf(F, V(:,1), V(:,2), V(:,3), L)
end
