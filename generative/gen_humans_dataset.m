clear
close all
load('evalues.mat')
load('evectors.mat')
load('meanShape.mat')
load('facesShapeModel.mat')
load('human_lms.mat')
evectors = reshape(evectors, [4300 6449 3]);

% k = 6, r = 2, d = 3 ook proberen

tuple1 = {{[6, 7, 8, 9, 10],[1],[2],[3],[4],[5]},{[1, 7, 3, 4, 9],[2],[5],[6],[8],[10]...
    },{[1, 8, 4, 3, 7],[2],[5],[6],[9],[10]},{[2, 5, 4, 7, 6],[1],[3],[8],[9],[10]},{[2, 10, 9, 4, 5],...
    [1],[3],[6],[7],[8]}};
scalef =0.0045*5e2;
minAGD = 1110;

normF = 2000;
V = points / normF;
V = V - mean(V, 1);
inds = [inds(1) inds(3) inds(5) inds(4) inds(2)];
gluer = Gluer(V,faces,inds,tuple1,minAGD);
flattenerO = Torus_Flattener(gluer.V_torus,gluer.T_torus);
cutMeshO = CutMesh(gluer, flattenerO);
for ii=1:2000
    ii
    weights= ((rand(4300,1)*2 -1).*sqrt(evalues)') * scalef;
    variation= squeeze(sum(weights.*evectors,1)) / normF;
    newpoints = V + variation;
%     lms = newpoints(inds, :);
%     visualize(faces, newpoints, lms);
    
    flattener = Torus_Flattener(gluer.V_torus,gluer.T_torus);
    cutMesh = CutMesh(gluer, flattener);

    
    params.sz = 128;
    pushed_function = push_functions_to_flattening_AE(cutMesh, newpoints, params);
    save(strcat('generative/meshes/humans/h', num2str(ii)), 'pushed_function')
end