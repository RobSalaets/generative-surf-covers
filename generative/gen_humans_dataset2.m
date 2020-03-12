clear
% close all
% load('caesar_VT.mat')

tuple1 = {{[6, 7, 8, 9, 10],[1],[2],[3],[4],[5]},{[1, 7, 3, 4, 9],[2],[5],[6],[8],[10]...
    },{[1, 8, 4, 3, 7],[2],[5],[6],[9],[10]},{[2, 5, 4, 7, 6],[1],[3],[8],[9],[10]},{[2, 10, 9, 4, 5],...
    [1],[3],[6],[7],[8]}};
tuple3 = {{[1,2]},{[1,2]},{[1,2]},{[1,2]}};
tuple = tuple1;


for    jj=[1,2,3,7]
    [V,T] = read_ply(sprintf('jitter_%03d.ply', jj));
    %%
    normF = 1;
    V = V / normF;
    V = V - mean(V, 1);
    k = length(tuple);
    params.n = 5;
    params.nFarthest = k;
    params.doplot = 0;
%     [cones, AGD] = getPointTripletsByFarthestPointSampling(V, T, params);
%     %%
%     [cones, min_agd_point] = sort_cones_in_plane(V,T, cones, AGD);
    cones = [2447,3300,6700,5906,412];
%     cones = circshift(cones,randi(5));
    min_agd_point = 1785;
%     visualize_with_lms(T,V,[], [cones, min_agd_point])

    %%
    gluer = Gluer(V,T,cones,tuple,min_agd_point);
    %%
    target = [0.00282030927717228,0.334446865442616,-0.0923488374371154];
    [~, iloc] = min(vecnorm((V-repmat(target,[size(V,1) 1]))'));
    idx = find(gluer.torus_to_sphere == iloc);
    flattenerO = Torus_Flattener(gluer.V_torus,gluer.T_torus, idx);
    cm = CutMesh(gluer, flattenerO, 1);
    tr = triangulation(T,V);

    params.sz = 128;
    [pushed_function] = push_functions_to_flattening_edit(cm, V, params);
    %%
    processed = transform_flat_torus(pushed_function, [0;0]/512, 0);
    processed = pushed_function;
    figure 
    imagesc((processed + 1)*.7)
    
end
    %%
%     jj=2;
%     [V2,T2] = read_ply(sprintf('jitter_%03d.ply', jj));
%     
%     normF = 1;
%     V2 = V2 / normF;
%     V2 = V2 - mean(V2, 1);
%     k = length(tuple);
%     params.n = 5;
%     params.nFarthest = k;
%     params.doplot = 0;
%     [cones2, AGD2] = getPointTripletsByFarthestPointSampling(V2, T2, params);
%     %%
%     [cones2, min_agd_point2] = sort_cones_in_plane(V2,T2, cones2, AGD2);
%     visualize_with_lms(T2,V2,[], [cones2, min_agd_point2])
% 
%     %%
%     gluer2 = Gluer(V2,T2,cones2,tuple,min_agd_point2);
%     %%
%     idx2 = find(gluer2.torus_to_sphere == max(cones2));
%     idx2 = find(gluer2.torus_to_sphere == 3000);
%     flattenerO2 = Torus_Flattener(gluer2.V_torus,gluer2.T_torus, idx2);
%     cm2 = CutMesh(gluer2, flattenerO2, 1);
%     tr = triangulation(T2,V2);
% 
%     params.sz = 512;
%     [pushed_function2] = push_functions_to_flattening_edit(cm2, V2, params);
%     %%
%     processed = transform_flat_torus(pushed_function2, [0;0]/512, 0);
%     processed = pushed_function2;
%     figure 
%     imagesc((processed + 1)*.7)

% 
% first_cone = cm2.V2(cm2.inds_mesh_divided_to_inds_plane{cones2(1)},:)

%%

figure
scatter(mod(cm.V(:,1),1),mod(cm.V(:,2),1), '.')

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
