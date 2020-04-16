clear


load('meanShape.mat')
load('facesShapeModel.mat')
load('human_lms.mat')
close all
% load('reg_lms.mat')
% [V,faces] = read_ply('tr_reg_020.ply');




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

params = struct();
params.sz = 64;
[pushed_function_og] = push_functions_to_flattening_edit(cm, V, params);
%%
close all
%Equidistant topology
sz = size(pushed_function_og, 1);
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
% visualize_with_lms(T_grid, reshape(pushed_function_og, [sz*sz 3]))
%%
for kk=1:10
load(sprintf('C:\\Users\\Rob\\Desktop\\Thesis\\Geometry\\progressive_growing_of_gans\\generated\\fm128-zm-2355\\G%i.mat', kk-1))
load('hstats8000.mat')
pushed_function = double(pushed_function)*1.2 + mn;
% load('h11.mat')
% load('512h1.mat')

params.sz = size(pushed_function, 1);

all_scales = zeros(size(cm.V_divided,1),10);
for ii = 1:size(cm.V_divided,1)
    ixs = cm.inds_mesh_divided_to_inds_plane{ii};
    scales = cm.vertex_scale(ixs);
    scales = scales./sum(scales);
    all_scales(ii,1:length(scales)) = scales';
    [maxS(ii), mloc] = max(scales);
    max_scale_idx(ii) = ixs(mloc);
end


[X,Y] = meshgrid(linspace(0,1-1/params.sz,params.sz),linspace(0,1-1/params.sz, params.sz));
u = mod(cm.V(:,1),1.0);
v = mod(cm.V(:,2),1.0);
interpV(:,1) = interp2(X,Y,pushed_function(:,:,1),u,v, 'spline');
interpV(:,2) = interp2(X,Y,pushed_function(:,:,2),u,v, 'spline');
interpV(:,3) = interp2(X,Y,pushed_function(:,:,3),u,v, 'spline');
% figure
% pcshow(interpV)
% title('Point cloud - Flat torus')
%Max
reconVmax = interpV(max_scale_idx,:);
reconVavg = reconVmax;
reconVavg2 = reconVavg;
%AVG
for ii = 1:size(cm.V_divided,1)
    ixs = cm.inds_mesh_divided_to_inds_plane{ii};
    reconVavg(ii,:) = all_scales(ii, 1:length(ixs)) * interpV(ixs,:);
    reconVavg2(ii,:) = all_scales(ii, 1:3) * interpV(ixs(1:3),:) / sum(all_scales(ii, 1:3));
end

% visualize_with_lms(faces, reconVmax, [], {}, 0,30,0)
% title('Max')
% visualize_with_lms(faces, reconVavg, [], {}, 0,30,0)
% title('Avg')
visualize_with_lms(faces, reconVavg2, [], {}, 0,30,0)
title('Avg 3 comp')
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

function [dataFunctions] = push_functions_to_flattening_edit(cutMesh, functions, params)

    params.null = [];
    sz = getoptions(params,'sz', 512);
    numFunctions = size(functions,2);

    for ii = 1:length(cutMesh.divided_edges)
        functions = [functions ; ...
            (functions(cutMesh.divided_edges{ii}(:,1),:) + functions(cutMesh.divided_edges{ii}(:,2),:)) / 2];
    end

    % process first function and set good k (how many tiles to put on each side)
    has_nan = 1;
    k = 1;
%     while any(has_nan(:))
        disp([datestr(datetime('now')) ' tiling with k=', num2str(k)]);
        [V_merged, ~, f, vals] = tile(cutMesh, functions, k);
        X = linspace(0, 1-1/sz, sz);
        Y = linspace(0, 1-1/sz, sz);    
        SI = scatteredInterpolant(V_merged,f(vals));
        [mX, mY] = meshgrid(X,Y);
        out2 = SI(mX,mY);
%         [out,tn,al2,al3] = mytri2grid(V_merged', T_merged', f(vals), X, Y);
%         dataFunctions(:,:,1) = out;
        dataFunctions(:,:,1) = out2;
        has_nan = isnan(dataFunctions);
        assert(~any(has_nan(:)));
%         k = k + 1;
%     end

    for ii=2:numFunctions
        f=functions(:,ii);
%         [out,tn,al2,al3] = mytri2grid(V_merged', T_merged', f(vals), tn, al2, al3);
        
        SI = scatteredInterpolant(V_merged,f(vals));
        out2 = SI(mX,mY);
%         dataFunctions(:,:,ii) = out;
        dataFunctions(:,:,ii) = out2;
    end
end


function [V_merged, T_merged, f, vals] = tile(cutMesh, functions, k)
    tot_tiles = (k*2+1)^2;

    f = functions(:,1);
    [V_merged, T_merged] = make_tiling(cutMesh, k);
    vals = repmat(cutMesh.inds_plane_to_divided_inds_mesh, tot_tiles, 1);
%     seg = seg(cutMesh.dividedTs2Ts);
%     seg_merged = repmat(seg, (size(cutMesh.T, 1) / size(seg, 1)) * tot_tiles, 1);
end
