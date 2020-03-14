clear
% close all

tuple = {{[6, 7, 8, 9, 10],[1],[2],[3],[4],[5]},
    {[1, 7, 3, 4, 9],[2],[5],[6],[8],[10]...
    },
    {[1, 8, 4, 3, 7],[2],[5],[6],[9],[10]},
    {[2, 5, 4, 7, 6],[1],[3],[8],[9],[10]},
    {[2, 10, 9, 4, 5],...
    [1],[3],[6],[7],[8]}
    };
v1=[0.918936237897646,0.350783089827239,0.180298126916992];
v2=[0.163472145899359,-0.908418397022596,-0.384776394112160];
for m = 0:5
    m
    for ii = 0:0
        [V,T] = read_ply(sprintf('tr_%d_%03d.ply',m, ii));
        V = V - mean(V, 1);
        k = length(tuple);
        params.n = 5;
        params.nFarthest = k;
        [cones, AGD] = getPointTripletsByFarthestPointSampling(V, T, params);
        [cones, min_agd_point] = sort_cones_in_plane(V,T, cones, AGD, [0 1 0]);
        if m == 10
            cones = flip(cones);
        end
        gluer = Gluer(V,T,cones,tuple,min_agd_point);
        %%
        flattenerO = Torus_Flattener(gluer.V_torus,gluer.T_torus);
        cm = CutMesh(gluer, flattenerO, 1);
        %%
        maps = sort(cm.inds_mesh_divided_to_inds_plane{cones(1)});
        [loop{1}, lps{1}] = trace_loops_flat_torus(flattenerO, gluer, maps(1), 1,[]);
        [loop{2}, lps{2}] = trace_loops_flat_torus(flattenerO, gluer, maps(1), 2, loop{1});
        if ii == 0 && m==0
            v1 = normr(sum(normr(gluer.V_torus(loop{1}(2:4),:) - repmat(gluer.V_torus(loop{1}(1),:), [3 1])),1));
            v2 = normr(sum(normr(gluer.V_torus(loop{2}(2:4),:) - repmat(gluer.V_torus(loop{2}(1),:), [3 1])),1));
        else
            w1 = normr(sum(normr(gluer.V_torus(loop{1}(2:4),:) - repmat(gluer.V_torus(loop{1}(1),:), [3 1])),1));
            w2 = normr(sum(normr(gluer.V_torus(loop{2}(2:4),:) - repmat(gluer.V_torus(loop{2}(1),:), [3 1])),1));
            if abs(v1*w2') > abs(v1*w1')
                tmp = loop{1};
                loop{1} = loop{2};
                loop{2} = tmp;
                tmp = lps{1};
                lps{1} = lps{2};
                lps{2} = tmp;
            end
        end
        %%
        figure
        quickscatter2d(cm.V,0, '.')
        hold on
        quickscatter2d(cm.V(lps{1},:),0, 'o')
        quickscatter2d(cm.V(lps{2},:),0, 'o')
        
        visualize_with_lms(gluer.T_torus, gluer.V_torus, [], {loop{1} loop{2} find(gluer.torus_to_sphere == min_agd_point)})
        
        flattenerP = Torus_Flattener(gluer.V_torus,gluer.T_torus, loop);
        cm = CutMesh(gluer, flattenerP, 1);
        %%
        params.sz = 128;
        [pushed_function] = push_functions_to_flattening_edit(cm, V, params);
       
        figure 
        imagesc((pushed_function + 1)*.7)
        flattenerP.cut_vertex
    end
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
