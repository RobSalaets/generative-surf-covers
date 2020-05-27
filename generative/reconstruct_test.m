clear
close all
load('evalues.mat')
load('evectors.mat')
load('meanShape.mat')
load('facesShapeModel.mat')
load('human_lms.mat')
% load('reg_lms.mat')
% [V,faces] = read_ply('tr_reg_020.ply');

tuple1 = {{[6, 7, 8, 9, 10],[1],[2],[3],[4],[5]},{[1, 7, 3, 4, 9],[2],[5],[6],[8],[10]...
    },{[1, 8, 4, 3, 7],[2],[5],[6],[9],[10]},{[2, 5, 4, 7, 6],[1],[3],[8],[9],[10]},{[2, 10, 9, 4, 5],...
    [1],[3],[6],[7],[8]}};

% scalef =0.0045*5e2;
minAGD = 1110;
normF = 1000;
V = points / normF;
V = V - mean(V, 1);
inds = [inds(1) inds(3) inds(5) inds(4) inds(2)];
% evectors = reshape(evectors, [4300 6449 3]);
% weights= ((rand(4300,1)*2 -1).*sqrt(evalues)') * scalef*1.3;
% variation= squeeze(sum(weights.*evectors,1)) / normF;
% V = V + variation;

ivert = inds(1);

% BASE MODEL

[cm, gluer, flattenerO, v1, v2] = get_consistent_mapping(tuple1, V, faces, inds, minAGD, ivert);

params = struct();
params.sz = 512;
[pushed_function_og] = push_functions_to_flattening_edit(cm, V, params);
visualize_with_lms(faces, V, [], inds)
% Get lms
lms_plane = {};
for ii=1:length(inds)
    icmV = cm.inds_mesh_divided_to_inds_plane_unsorted{inds(ii)};
    lms_plane{ii} = {cm.V(icmV,:), inds(ii)};
end
%% Uncorrelated mesh
params.n = 5;
params.nFarthest = 5;
params.doplot = 0;

inds = [6619,5907,3217,2447,413];
[V2,faces2] = read_ply('tr_reg_010.ply');
V2 = V2 - mean(V2, 1);
V2 = [V2(:,1) -V2(:,3) V2(:,2)];
[cones, AGD] = getPointTripletsByFarthestPointSampling(V2, faces2, params);
[cones, min_agd_point] = sort_cones_in_plane(V2,faces2, cones, AGD);

ivert = cones(end);
visualize_with_lms(faces2, V2, [], cones)
%%
[cm2, gluer2, flattener2] = get_consistent_mapping(tuple1, V2, faces2, cones, min_agd_point, ivert, v1, v2);
[pushed_function] = push_functions_to_flattening_edit(cm2, V2, params);

figure
imagesc(0.5+pushed_function_og)
title('og')
figure
imagesc(0.5+pushed_function)

%%
% V = points / normF;
% V = V - mean(V, 1);
% [cm3, ~, ~] = get_consistent_mapping(tuple1, V, faces, inds, minAGD, ivert, v1, v2);
%%
figure
quickscatter2d(cm.V, 1, '.')
hold on
quickscatter2d(cm2.V, 1, '.')
% quickscatter2d(cm3.V, 1, '.')
%%
figure(1)
for kk=1:32
load(sprintf('C:\\Users\\Rob\\Desktop\\Thesis\\Geometry\\progressive_growing_of_gans\\generated\\512-fmb2048-zm-3993\\Gs%i.mat', kk-1))
% load(sprintf('bh (1%03i)', kk))
% pushed_function = pushed_function_og;
% load(sprintf('zmh%i', kk))
load('hstats10k.mat')

recon = reconstruct_function(cm2, pushed_function+pf_avg);

visualize_with_lms(faces, recon{3}, [], {}, 0,30,0)
title('Recon Avg')
% visualize_with_lms(faces, recon{2}, [], {}, 0,30,0)
% title('Recon Max BC')

end