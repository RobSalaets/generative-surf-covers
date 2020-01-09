function visualize_with_lms(faces, newpoints, lms_idx, doAvg)

    V = newpoints;
    if doAvg==1
        avgs = [];
        for i=1:length(lms_idx)
            [row, ~]  = find(faces==lms_idx(i));
            ids = faces(row,:);
            avgs(i,:) = mean(V(ids(:),:), 1);
            V(lms_idx(i),:) = avgs(i,:);
        end
    end

    figure
    trisurf(triangulation(faces, V));
    axis vis3d
    axis equal
    hold on
    lms = V(lms_idx,:);
    xlim([-0.5 0.5]);
    ylim([-0.5 0.5]);
    zlim([-0.5 0.5]);
    scatter3(lms(:,1), lms(:,2),lms(:,3), 'filled')

    
    if length(lms_idx) < 250
    
        text(lms(:,1), lms(:,2), lms(:,3), cellstr(num2str(lms_idx(:))))
    end
end