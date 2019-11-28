clear
close all
load('names.mat');
load('rotmat.mat');
load('minAGD_points.mat');
load('landmarks.mat');
n = length(names_cell);

tuple1 = {{[6, 7, 8, 9, 10],[1],[2],[3],[4],[5]},{[1, 7, 3, 4, 9],[2],[5],[6],[8],[10]...
    },{[1, 8, 4, 3, 7],[2],[5],[6],[9],[10]},{[2, 5, 4, 7, 6],[1],[3],[8],[9],[10]},{[2, 10, 9, 4, 5],...
    [1],[3],[6],[7],[8]}};


for ii = 1:n
    ii
    fname = names_cell{ii};
    
    [V,F] = read_ply(strcat('meshes/all/',fname));
    [V,F] = preprocess(V,F,R,ii);
    lms = landmarks(ii,:);
    minAGD = minAGD_points(ii);
    [cutMesh, ~, ~] = flatten_sphere(V ,F, lms, minAGD, tuple1);
    f = V;
    params.sz = 100;
    pushed_function = push_functions_to_flattening_AE(cutMesh, f, params);
    save(strcat('generative/meshes/', num2str(ii)), 'pushed_function')
%     newF = zeros((params.sz-1)^2*2, 3);
%     jj = 1;
%     for xx = 1:99
%         for yy = 1:99 
%             newF(jj*2-1,:) = [(xx + (yy-1)*100) (xx +1+ (yy-1)*100) (xx + (yy)*100)];
%             newF(jj*2,:) = [(xx +1+ (yy-1)*100) (xx +1+ (yy)*100) (xx + (yy)*100)];
%             
%             jj = jj+1;
%         end
%     end
%     trisurf(triangulation(newF, reshape(pushed_function, [100*100 3])));
    
end

%---------------
% tuple1 = {{[6, 7, 8, 9, 10],[1],[2],[3],[4],[5]},{[1, 7, 3, 4, 9],[2],[5],[6],[8],[10]...
%     },{[1, 8, 4, 3, 7],[2],[5],[6],[9],[10]},{[2, 5, 4, 7, 6],[1],[3],[8],[9],[10]},{[2, 10, 9, 4, 5],...
%     [1],[3],[6],[7],[8]}};
% 
% 
% 
% k = length(tuple1);
% 
% params.n = 5;
% params.nFarthest = k;
% params.doplot = 1;
% params.T = 0.2;
% params.constraint = @pos;
% 
% [cones, AGD] = getPointTripletsByFarthestPointSampling(V, T, params);
% [~, min_agd_point] = min(AGD);
% 
% [cutMesh, ~, ~] = flatten_sphere(V ,T, cones, min_agd_point, tuple1);
% 
% pushed_function = push_functions_to_flattening_AE(cutMesh, f);
% 
% 
% figure
% image((P+1.0)./2.0)

function [Vr, Fr] = preprocess(V, F, R, ii)
    Vr = V;
    Vr(:,3) = -Vr(:,3);
    [Vr, Fr] = delaunayize(Vr,F);
    m = mean(Vr, 1);
    Vr = Vr-m;
    Vr =(R(:,:,ii)*Vr')';
end


