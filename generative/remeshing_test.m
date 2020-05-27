clear
close all

tuple = {{[6, 7, 8, 9, 10],[1],[2],[3],[4],[5]},
    {[1, 7, 3, 4, 9],[2],[5],[6],[8],[10]...
    },
    {[1, 8, 4, 3, 7],[2],[5],[6],[9],[10]},
    {[2, 5, 4, 7, 6],[1],[3],[8],[9],[10]},
    {[2, 10, 9, 4, 5],...
    [1],[3],[6],[7],[8]}
    };

cone_vec = [-1 0 0];
load('meanShape.mat')
load('facesShapeModel.mat')
V = points / 1000; V = V - mean(V, 1); T = faces; V(:,2) = -V(:,2);
params.n = 5;
params.nFarthest = 5;
[cones, AGD] = getPointTripletsByFarthestPointSampling(V, T, params);
[cones, min_agd_point] = sort_cones_in_plane(V,T, cones, AGD, cone_vec);
cones = flip(cones);
% visualize_with_lms(T,V, [], [cones, min_agd_point])
ivert = cones(1);
[cmO, gluer, flattener, v1, v2, lms_plane] = get_consistent_mapping(tuple, V, T, cones, min_agd_point, ivert);

params = struct();
params.sz = 512;
[pushed_function_og] = push_functions_to_flattening_edit(cmO, V, params);


%%
figure 
quickscatter2d(cmO.V, 1, '.',10)
hold on

for m = 0:1
    m
    for ii = 1:1
        [V,T] = read_ply(sprintf('tr_%d_%03d.ply',m, ii));
        if ii == 1
            [V,T] = read_ply(sprintf('tr_reg_%03d.ply',m*10));
        end
        V = [V(:,1), V(:,3), V(:,2)];
        V = V - mean(V, 1);
        params.n = 5;
        params.nFarthest = 5;
        [cones, AGD] = getPointTripletsByFarthestPointSampling(V, T, params);
        [cones, min_agd_point] = sort_cones_in_plane(V,T, cones, AGD, cone_vec);
        cones = flip(cones);
%         visualize_with_lms(T,V, [], [cones, min_agd_point])
        ivert = cones(1);
        
        [cm, gluer, flattener] = get_consistent_mapping(tuple, V, T, cones, min_agd_point, ivert, v1, v2, lms_plane);
        recon = reconstruct_function(cm, pushed_function_og);
%         
%         visualize_with_lms(T, recon{3}, [], {}, 0,30,0)
%         title(sprintf('tr m%i i%i', m, ii))
    end
end
