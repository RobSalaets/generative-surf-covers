function visualize_with_lms(faces, newpoints, lms_idx)
    figure
    trisurf(triangulation(faces, newpoints));
    axis vis3d
    axis equal
    hold on
    lms = newpoints(lms_idx,:);
    scatter3(lms(:,1), lms(:,2),lms(:,3), 'filled')
end