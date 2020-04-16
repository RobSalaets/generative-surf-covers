clear
close all
load('evalues.mat')
load('evectors.mat')
load('meanShape.mat')
load('facesShapeModel.mat')
load('human_lms.mat')
% load('reg_lms.mat')
% [V,faces] = read_ply('tr_reg_020.ply');

%%%%%
% LOOK AT lift_image_AE to ''improve'' reconstruction
%%%%%


tuple1 = {{[6, 7, 8, 9, 10],[1],[2],[3],[4],[5]},{[1, 7, 3, 4, 9],[2],[5],[6],[8],[10]...
    },{[1, 8, 4, 3, 7],[2],[5],[6],[9],[10]},{[2, 5, 4, 7, 6],[1],[3],[8],[9],[10]},{[2, 10, 9, 4, 5],...
    [1],[3],[6],[7],[8]}};
% inds = [6619,5907,3217,2447,413];
scalef =0.0045*5e2;
minAGD = 1110;
normF = 1000;
V = points / normF;
V = V - mean(V, 1);
inds = [inds(1) inds(3) inds(5) inds(4) inds(2)];
evectors = reshape(evectors, [4300 6449 3]);
% weights= ((rand(4300,1)*2 -1).*sqrt(evalues)') * scalef*1.3;
% save('tmpweights', 'weights')
load('tmpweights.mat')
variation= squeeze(sum(weights.*evectors,1)) / normF;
V = V + variation;

ivert = inds(1);
[cm, gluer, flattenerO, v1, v2] = get_consistent_mapping(tuple1, V, faces, inds, minAGD, ivert);

% gluer = Gluer(V,faces,inds,tuple1,minAGD);
% flattenerO = Torus_Flattener(gluer.V_torus,gluer.T_torus);
% cm = CutMesh(gluer, flattenerO, 1);

params = struct();
params.sz = 512;
[pushed_function_og] = push_functions_to_flattening_edit(cm, V, params);
visualize_with_lms(faces, V, [], {}, 0,30,0)
%% Get lms

lms_plane = {};
for ii=1:length(inds)
    icmV = cm.inds_mesh_divided_to_inds_plane_unsorted{inds(ii)};
    lms_plane{ii} = {cm.V(icmV,:), inds(ii)}
end
%%
V = points / normF;
V = V - mean(V, 1);
[cm2, gluer2, flattener2] = get_consistent_mapping(tuple1, V, faces, inds, minAGD, ivert, v1, v2, lms_plane);
figure
quickscatter2d(cm.V, 1, '.')
hold on
quickscatter2d(cm2.V, 1, '.')
%%
V = points / normF;
V = V - mean(V, 1);
[cm3, ~, ~] = get_consistent_mapping(tuple1, V, faces, inds, minAGD, ivert, v1, v2);


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
drawnow
pause(0.1)
frame = getframe(1);
im = frame2im(frame);
[imind,qsdf] = rgb2ind(im,256);
if kk == 1
  imwrite(imind,qsdf,'testgif2.gif','gif', 'Loopcount',inf);
else
  imwrite(imind,qsdf,'testgif2.gif','gif','WriteMode','append');
end
% visualize_with_lms(faces, recon{2}, [], {}, 0,30,0)
% title('Recon Max wBC')

% recon2 = reconstruct_function(cm3, pushed_function_og);
% visualize_with_lms(faces, recon2{3}, [], {}, 0,30,0)
% title('Recon Avg woBC')
% visualize_with_lms(faces, recon2{2}, [], {}, 0,30,0)
% title('Recon Max woBC')


% figure
% imagesc(pushed_function_og+0.5)
% title('og + var')
% figure
% imagesc(pushed_function + 0.5)
% title('bh')
% figure
% quickscatter2d(cm.V, 1, '.')
% hold on
% quickscatter2d(cm2.V, 1, '.')
% hold on
% quickscatter2d(cm3.V, 1, '.')
end