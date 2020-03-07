clear
% close all

load('meanShape.mat')
load('facesShapeModel.mat')
load('human_lms.mat')

% k = 6, r = 2, d = 3 ook proberen

tuple1 = {{[6, 7, 8, 9, 10],[1],[2],[3],[4],[5]},{[1, 7, 3, 4, 9],[2],[5],[6],[8],[10]...
    },{[1, 8, 4, 3, 7],[2],[5],[6],[9],[10]},{[2, 5, 4, 7, 6],[1],[3],[8],[9],[10]},{[2, 10, 9, 4, 5],...
    [1],[3],[6],[7],[8]}}; %d=10
tuple2 = {{[3,4,5],[1],[2]},{[2,3,5],[1],[4]},{[1,5,2],[3],[4]},{[1,2,5],[3],[4]},{[2,4,3],[1],[5]}};%d=5
tuple3 = {{[1,2]},{[1,2]},{[1,2]},{[1,2]}}; %d=2
tuple4 = {{[1,2,3],[4]},{[2,3,4],[1]},{[1,2,4],[3]},{[1,2,3],[4]}};
inds = [inds(1) inds(3) inds(5) inds(4) inds(2)];
% inds = inds(2:end);
tuple = tuple1;

scalef =0.0045*5e2;
minAGD = 1110;

normF = 2000;
V = points / normF;
V = V - mean(V, 1);
gluer = Gluer(V,faces,inds,tuple,minAGD);
flattenerO = Torus_Flattener(gluer.V_torus,gluer.T_torus);
cm = CutMesh(gluer, flattenerO, 1);

%%

all_scales = zeros(size(cm.V_divided,1),10);
for ii = 1:size(cm.V_divided,1)
    ixs = cm.inds_mesh_divided_to_inds_plane{ii};
    scales = cm.vertex_scale(ixs);
    [maxS(ii), mloc] = max(scales);
    max_scale_idx(ii) = ixs(mloc);
    scales = scales./sum(scales);
    all_scales(ii,1:length(scales)) = scales';
    
end
meanMaxS = mean(maxS)
[minMaxS, mms_iloc] = min(maxS)
figure
histogram(maxS, 100)

%%
sm_idx = maxS > 0.5;
ids = max_scale_idx(sm_idx);
figure
scatter(mod(cm.V(:,1),1), mod(cm.V(:,2),1), 20, cm.vertex_scale, '.')
hold on
scatter(mod(cm.V(ids,1),1), mod(cm.V(ids,2),1), 20, 'red', '.')

visualize_with_lms(faces,V,[],  find(sm_idx(1:6449)))
