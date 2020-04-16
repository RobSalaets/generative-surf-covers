function recon = reconstruct_function(cm, pushed_function, varargin)
    if length(varargin) > 0
        pushed_function = pushed_function + varargin{1};
    end
    params.sz = size(pushed_function, 2);

    all_scales = zeros(size(cm.V_divided,1),10);
    for ii = 1:size(cm.V_divided,1)
        ixs = cm.inds_mesh_divided_to_inds_plane{ii};
        scales = cm.vertex_scale(ixs);
        scales = scales./sum(scales);
        all_scales(ii,1:length(scales)) = scales';
        [maxS(ii), mloc] = max(scales);
        max_scale_idx(ii) = ixs(mloc);
    end
    [X,Y] = meshgrid(linspace(0,1,params.sz+1),linspace(0,1, params.sz+1));
    u = mod(cm.V(:,1),1.0);
    v = mod(cm.V(:,2),1.0);
    aug_pf = [pushed_function pushed_function(:,1,:); pushed_function(1,:,:) pushed_function(1,1,:) ];
    interpV(:,1) = interp2(X,Y,aug_pf(:,:,1),u,v, 'spline');
    interpV(:,2) = interp2(X,Y,aug_pf(:,:,2),u,v, 'spline');
    interpV(:,3) = interp2(X,Y,aug_pf(:,:,3),u,v, 'spline');
    recon{1} = interpV;
    %Max
    recon{2} = interpV(max_scale_idx,:);
    reconVavg = interpV(max_scale_idx,:);
    %AVG
    for ii = 1:size(cm.V_divided,1)
        ixs = cm.inds_mesh_divided_to_inds_plane{ii};
        reconVavg(ii,:) = all_scales(ii, 1:length(ixs)) * interpV(ixs,:); 
    end
    recon{3} = reconVavg;
end