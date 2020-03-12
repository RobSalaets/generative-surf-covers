function visualize_with_lms(faces, newpoints, varargin)
    
    nargs = length(varargin);
    V = newpoints;
    lms_idx = [];
    if nargs > 1
        lms_idx = varargin{2};
    end
    
    if nargs > 2 && varargin{3}==1
        avgs = [];
        for i=1:length(varargin())
            [row, ~]  = find(faces==lms_idx(i));
            ids = faces(row,:);
            avgs(i,:) = mean(V(ids(:),:), 1);
            V(lms_idx(i),:) = avgs(i,:);
        end
    end

    figure
    if nargs > 0 && ~isempty(varargin{1})
        trisurf(faces, V(:,1),V(:,2),V(:,3), 'FaceColor','interp', 'facevertexcdata', varargin{1})
    else
        trisurf(triangulation(faces, V));
    end
    axis vis3d
    axis equal
    hold on
    
    
    
    maxdim = max(max(abs(V)));
    xlim([-maxdim maxdim]);
    ylim([-maxdim maxdim]);
    zlim([-maxdim maxdim]);
    
    if iscell(lms_idx)
        for ii=1:length(lms_idx)
            lms = V(lms_idx{ii},:);
            scatter3(lms(:,1), lms(:,2),lms(:,3), 'filled')
        end
    else
        lms = V(lms_idx,:);
        scatter3(lms(:,1), lms(:,2),lms(:,3), 'filled', 'red')
        if length(lms_idx) < 50
            text(lms(:,1), lms(:,2), lms(:,3), cellstr(num2str(lms_idx(:))),'Color','green')
        end
    end
end