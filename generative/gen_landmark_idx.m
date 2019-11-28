clear
close all


names_cell = dir('generative/meshes/all/*.ply'); names_cell = {names_cell.name};
labels_cell = dir('generative/meshes/all/*.txt'); labels_cell = {labels_cell.name};
landmarks = zeros(size(names_cell, 2), 5);
for i = 1:size(names_cell, 2)
    fname = names_cell{i};
    flabel = labels_cell{i};
    flID = fopen(strcat('meshes/all/',flabel),'r');
    L = fscanf(flID,'%i\n');
    fclose(flID);
    [V,F] = read_ply(strcat('meshes/all/',fname));
    V(:,3) = -V(:,3);
    [V, F] = delaunayize(V,F);
    m = mean(V, 1);
    V = V-m;
    s = max(abs(V), [], 'all');
    V = (V)./s;
    
    land_mark_ids = zeros(5,1);
    for lm=2:6
        lmidx = find(L == lm);
        Flm = F(lmidx,:);
        map = unique(Flm(:));
        Vlm = V(unique(Flm(:)),:);
        for ii= 1:length(Flm)
            for jj = 1:3
                Flm(ii,jj) = find(map == Flm(ii,jj));
            end
        end
        
        adj = triangulation2adjacency_change(Flm,Vlm');
        dist = graphallshortestpaths(adj,'directed',false);
        faceAreas = computeSurfAreas(Vlm,Flm);
        nV = length(Vlm);
        oneRingAreas = zeros(nV,1);
        for ii = 1:nV
            ff = any(Flm'==ii); % indices of faces in the ii'th vertex 1-ring
            oneRingAreas(ii) = (1/3)*sum(faceAreas(ff));
        end
        AGD_raw = dist*oneRingAreas;
        [~, mIdx] = min(AGD_raw);
        land_mark_ids(lm-1) = map(mIdx);
    end
    landmarks(i,:) = land_mark_ids;
%     f = figure;
%     trisurf(F,V(:,1), V(:,2), V(:,3), L)
%     colorbar
%     uiwait(f);
end
save('landmarks.mat','landmarks');
save('names.mat', 'names_cell');
