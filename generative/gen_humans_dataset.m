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

normF = 1000;
V = points / normF;
V = V - mean(V, 1);


inds = [inds(1) inds(3) inds(5) inds(4) inds(2)];
gluer = Gluer(V,faces,inds,tuple1,minAGD);
% flattenerO = Torus_Flattener(gluer.V_torus,gluer.T_torus);
% cutMeshO = CutMesh(gluer, flattenerO);
flattener = Torus_Flattener(gluer.V_torus,gluer.T_torus);
cutMesh = CutMesh(gluer, flattener, 1);
for ii=1:3
    ii
    weights= ((rand(4300,1)*2 -1).*sqrt(evalues)') * scalef;
    variation= squeeze(sum(weights.*evectors,1)) / normF;
    newpoints = V + variation;
%     visualize_with_lms(faces, newpoints, [])
    
%     flattener = Torus_Flattener(gluer.V_torus,gluer.T_torus);
%     cutMesh = CutMesh(gluer, flattener, 1);

    params = struct();
    params.sz = 512;
    [pushed_function] = push_functions_to_flattening_edit(cutMesh, newpoints, params);
%     save(strcat('512h', num2str(ii)), 'pushed_function')
end
% profile viewer

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
