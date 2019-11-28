clear
close all

load('minAGD_points.mat')
load('names.mat')
n = length(names_cell);
minAGD_points = zeros(n,1);
for i = 43:size(names_cell, 2)
    i
    fname = names_cell{i};
    [V,F] = read_ply(strcat('meshes/all/',fname));
    V(:,3) = -V(:,3);
    [V, F] = delaunayize(V,F);
    m = mean(V, 1);
    V = V-m;
    adj = triangulation2adjacency_change(F,V');
    dist = graphallshortestpaths(adj,'directed',false);
    faceAreas = computeSurfAreas(V,F);
    nV = length(V);
    oneRingAreas = zeros(nV,1);
    for ii = 1:nV
        ff = any(F'==ii); % indices of faces in the ii'th vertex 1-ring
        oneRingAreas(ii) = (1/3)*sum(faceAreas(ff));
    end
    AGD_raw = dist*oneRingAreas;
    [~, mIdx] = min(AGD_raw);
    minAGD_points(i) = mIdx;
end
save('minAGD_points.mat','minAGD_points')