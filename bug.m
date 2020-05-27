clear
load('bug')

V = obj.V_cut_sphere(1:obj.l,:);

nc = [obj.seams{1}(end,1:2) obj.seams{2}(end,1:2) obj.seams{3}(end,1:2) obj.seams{4}(end,1:2) obj.seams{5}(end,1:2) ];
nc = unique(nc);
F = obj.T_cut_sphere(1:12894,:);

avgs = [];
for i=1:length(fn)
    [row, col]  = find(F==fn(i));
    ids = F(row,:);
    avgs(i,:) = mean(V(ids(:),:), 1);
    V(fn(i),:) = avgs(i,:);
%     V(fn(idx)-3:fn(idx)+3,:)= V(fn(idx)-3:fn(idx)+3,:) + 0.05*(repmat(avgs(i,:), [7 1])-V(fn(idx)-3:fn(idx)+3,:));
end

F2=[]
for i = 1:length(F)
    if length(setdiff(F(i,:), fn))<3
        F2 = [F2; F(i,:)];
    end
end

visualize_with_lms(F2,V, fn);



