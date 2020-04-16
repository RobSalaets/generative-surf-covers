clear
close all

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
tuple3 = {{[1,2]},{[1,2]},{[1,2]},{[1,2]}};
% inds = [6619,5907,3217,2447,413];
scalef =0.0045*5e2;
minAGD = 1110;
normF = 2000;
V = points / normF;
V = V - mean(V, 1);
inds = [inds(1) inds(3) inds(5) inds(4) inds(2)];
gluer = Gluer(V,faces,inds,tuple1,minAGD);
flattenerO = Torus_Flattener(gluer.V_torus,gluer.T_torus);
cm = CutMesh(gluer, flattenerO, 1);
%%
%Equidistant topology
sz = 128;
[Xgrid, Ygrid] = meshgrid(linspace(0, 1-1/sz, sz), linspace(0,1-1/sz, sz));
V_grid(:,1) = Xgrid(:);
V_grid(:,2) = Ygrid(:);
tsz = length(V_grid);
T_grid = [(1:tsz-sz)' (2:1+tsz-sz)' (2+sz:1+tsz)';
          (1:tsz-sz)' (2+sz:1+tsz)' (1+sz:tsz)'];
halfsz = tsz-sz;
T_grid(sz:sz:tsz-sz,:) = T_grid(sz:sz:tsz-sz,:)-[0 sz sz];
T_grid(halfsz+sz:sz:halfsz+tsz-sz,:) = T_grid(halfsz+sz:sz:halfsz+tsz-sz,:)-[0 sz 0];
T_grid = [T_grid; 
          (sz*sz-sz+1:sz*sz-1)' (sz*sz-sz+2:sz*sz)' (2:sz)';
          (sz*sz-sz+1:sz*sz-1)' (2:sz)' (1:sz-1)'];
trisurf(triangulation(T_grid, [V_grid zeros(length(V_grid), 1)]))
%%

for kk=1:5
load(sprintf('C:\\Users\\Rob\\Desktop\\Thesis\\Geometry\\progressive_growing_of_gans\\generated\\fm128-1800ticks-zm-fp16\\G%i.mat', kk))
load('hstats8000.mat')
pushed_function = double(pushed_function) + mn;
params.sz = sz;

u = V_grid(:,1);
v = V_grid(:,2);
% Not necassary is size are compatible
% interpV(:,1) = interp2(X,Y,pushed_function(:,:,1),u,v);
% interpV(:,2) = interp2(X,Y,pushed_function(:,:,2),u,v);
% interpV(:,3) = interp2(X,Y,pushed_function(:,:,3),u,v);
interpV = reshape(pushed_function, [sz*sz 3]);

%Get scale grid
si_scale = scatteredInterpolant(mod(cm.V,1),cm.vertex_scale);
scale_grid = reshape(si_scale(u, v), [sz sz]);
figure
imagesc(scale_grid)
%Max
% reconVmax = interpV(max_scale_idx,:);
% reconVavg = reconVmax;
% %AVG
% for ii = 1:size(cm.V_divided,1)
%     ixs = cm.inds_mesh_divided_to_inds_plane{ii};
%     reconVavg(ii,:) = all_scales(ii, 1:length(ixs)) * interpV(ixs,:); 
% end

% visualize_with_lms(faces, reconVavg)
visualize_with_lms(T_grid, interpV)
end

function s = perlin2D (m)
  s = zeros([m,m]);     % Prepare output image (size: m x m)
  w = 16;
  i = 4;
  while w > 3
    i = i + 1;
    d = interp2(randn([m,m]), i-1, 'spline');
    s = s + i * d(1:m, 1:m);
    w = w - ceil(w/2 - 1);
  end
  s = (s - min(min(s(:,:)))) ./ (max(max(s(:,:))) - min(min(s(:,:))));
end