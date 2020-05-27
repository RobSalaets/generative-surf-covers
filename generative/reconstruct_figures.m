clear
close all
load('evalues.mat')
load('evectors.mat')
load('meanShape.mat')
load('facesShapeModel.mat')
load('human_lms.mat')

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

%DO AS GEN HUMAN DATASET
gluer = Gluer(V,faces,inds,tuple1,minAGD);
flattener = Torus_Flattener(gluer.V_torus,gluer.T_torus);
cm = CutMesh(gluer, flattener, 1);
params = struct();
%%
params.sz = 256;
[pushed_function_og] = push_functions_to_flattening_edit(cm, V, params);
% flat_tr = triangulation(cm.T, cm.V);
% triplot(flat_tr)
%%
load(sprintf('C:\\Users\\Rob\\Desktop\\Thesis\\Geometry\\progressive_growing_of_gans\\generated\\512-fmb2048-zm-3993\\Gs%i.mat', 1))
% load(sprintf('bh (1%03i)', 500))
% pushed_function = pushed_function_og;
% load(sprintf('zmh%i', 200))
load('hstats10k.mat')
load('hstats8000.mat')

recon = reconstruct_function(cm,pushed_function+ pf_avg);
% recon = reconstruct_function(cm, pushed_function+mn);
% recon = reconstruct_function(cm, pushed_function_og);


visualize_with_lms(faces, recon{3}, [], {},0,30,45)
% title('Reconstructie - Gewogen gemiddelde')
xlim([-0.3 0.3])
ylim([-0.3 0.3])
zlim([0.3 0.9])
% axis off
visualize_with_lms(faces, recon{2}, [], {},0,30,45)
% title('Reconstructie - Beste hoekpunt')
xlim([-0.3 0.3])
ylim([-0.3 0.3])
zlim([0.3 0.9])
% axis off

%%
a = pi/4;
R = [cos(a) -sin(a) 0; sin(a) cos(a) 0; 0 0 1];
f=figure('units','normalized','outerposition',[0 0 1 1]);
axes('Units', 'normalized', 'Position', [0 0 1 1])
f.WindowState = 'maximized';
V = points / normF;
V = V - mean(V, 1);
visualize_with_lms(faces, V*R, [], {},0,30,45)
hold on
ii = 1;
for sz = [64 128 256 512]
    params.sz = sz;
    [pushed_function_res] = push_functions_to_flattening_edit(cm, V, params);
    
    recon = reconstruct_function(cm, pushed_function_res);
    V2 = recon{3}*R + [0.6 0 0] * ii;
    h(ii)=trisurf(triangulation(faces, V2));
%     shading interp
    h(ii).FaceLighting = 'flat';
    h(ii).AmbientStrength = 0.7;
    h(ii).DiffuseStrength = 0.8;
    h(ii).SpecularStrength = 0.1;
    h(ii).SpecularExponent = 25;
    h(ii).FaceColor = [0.1 0.3 0.25];
    h(ii).EdgeColor = 'none';
    ii = ii + 1;
end
xlim([-0.5 ii*0.6-0.1])
axis off