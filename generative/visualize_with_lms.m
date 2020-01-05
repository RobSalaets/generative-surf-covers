function visualize_with_lms(faces, newpoints, lms_idx)
    figure
    trisurf(triangulation(faces, newpoints));
    axis vis3d
    axis equal
    hold on
    lms = newpoints(lms_idx,:);
    xlim([-0.5 0.5]);
    ylim([-0.5 0.5]);
    zlim([-0.5 0.5]);
    scatter3(lms(:,1), lms(:,2),lms(:,3), 'filled')
end