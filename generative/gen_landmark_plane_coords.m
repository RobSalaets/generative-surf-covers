clear
close all
load('names.mat');
load('rotmat.mat');
load('minAGD_points.mat');
load('landmarks.mat');
n = length(names_cell);

tuple1 = {{[6, 7, 8, 9, 10],[1],[2],[3],[4],[5]},{[1, 7, 3, 4, 9],[2],[5],[6],[8],[10]...
    },{[1, 8, 4, 3, 7],[2],[5],[6],[9],[10]},{[2, 5, 4, 7, 6],[1],[3],[8],[9],[10]},{[2, 10, 9, 4, 5],...
    [1],[3],[6],[7],[8]}};


lmlocs = zeros(n, 5,5);
lm2D = zeros(n, 5,5,2);
for ii = 1:n
    ii
    fname = names_cell{ii};
    
    [V,F] = read_ply(strcat('meshes/all/',fname));
    [V,F] = preprocess(V,F,R,ii);
    lms = landmarks(ii,:);
    minAGD = minAGD_points(ii);
    [cutMesh, ~, ~] = flatten_sphere(V ,F, lms, minAGD, tuple1);
    
    for jj = 1:5
        cell = cutMesh.inds_mesh_divided_to_inds_plane(lms(jj));
        locs = cell{1};
        [~,perm] = sort(cutMesh.vertex_scale(locs), 'ascend');
        locs = locs(perm);
        locs = locs(1:5);
        lmlocs(ii, jj,:)=locs;
        pl = cutMesh.V(locs,:);
        pl(pl<0) = pl(pl<0)+1;
        pl(pl>1) = pl(pl>1)-1;
        lm2D(ii,jj,:,:) = pl;
    end
end
save('generative/landmark_plane_coords.mat', 'lm2D','lmlocs')

function [Vr, Fr] = preprocess(V, F, R, ii)
    Vr = V;
    Vr(:,3) = -Vr(:,3);
    [Vr, Fr] = delaunayize(Vr,F);
    m = mean(Vr, 1);
    Vr = Vr-m;
    Vr =(R(:,:,ii)*Vr')';
end