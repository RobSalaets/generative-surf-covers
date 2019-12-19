clear
close all
load('evalues.mat')
load('evectors.mat')
load('meanShape.mat')
load('facesShapeModel.mat')
load('human_lms.mat')
evectors = reshape(evectors, [4300 6449 3]);

tuple1 = {{[6, 7, 8, 9, 10],[1],[2],[3],[4],[5]},{[1, 7, 3, 4, 9],[2],[5],[6],[8],[10]...
    },{[1, 8, 4, 3, 7],[2],[5],[6],[9],[10]},{[2, 5, 4, 7, 6],[1],[3],[8],[9],[10]},{[2, 10, 9, 4, 5],...
    [1],[3],[6],[7],[8]}};
scalef =0.0045*5e2;
minAGD = 1110;

for ii=1:500
    ii
    weights= ((rand(4300,1)*2 -1).*sqrt(evalues)') * scalef;
    variation= squeeze(sum(weights.*evectors,1));
    newpoints = points + variation;
    lms = newpoints(inds, :);
    visualize(faces, newpoints, lms);
    V = newpoints / 2000;
    V = V - mean(V, 1);
    F = faces;
    [cutMesh, ~, ~] = flatten_sphere(V ,F, inds, minAGD, tuple1);
    params.sz = 128;
    pushed_function = push_functions_to_flattening_AE(cutMesh, V, params);
    save(strcat('generative/meshes/humans/h', num2str(ii)), 'pushed_function')
    
end

function visualize(faces, newpoints, lms)
    figure
    trisurf(triangulation(faces, newpoints));
    axis vis3d
    axis equal
    hold on
    scatter3(lms(:,1), lms(:,2),lms(:,3))
    scatter3(newpoints(5396,1), newpoints(5396,2), newpoints(5396,3))
end