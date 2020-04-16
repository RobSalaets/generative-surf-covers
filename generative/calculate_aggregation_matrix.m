function [M] = calculate_aggregation_matrix(cutMesh, sze, avg_mode, varargin)
%Calculate the aggregation matrix based on the information in cutMesh.
%This assumes a tiling of 1 is sufficient, i.e. mod(V,1) has no overlap.
    M = sparse(sze^2,sze^2);
    do_threshold = 0;
    threshold = 0;
    if ~isempty(varargin)
        threshold = varargin{1};
        do_threshold = 1;
    end
    for v = 1:size(cutMesh.V_divided,1)
        v
        copy_indices = cutMesh.inds_mesh_divided_to_inds_plane{v};
        copy_pxs = floor(mod(cutMesh.V(copy_indices,:),1)*sze);
        copy_pxs = copy_pxs(:,1) + copy_pxs(:,2)*sze + 1;
        if avg_mode == 0
            for ii=1:length(copy_pxs)-1
                if do_threshold == 1 && cutMesh.vertex_scale(copy_indices(ii)) < threshold
                    continue
                end
                for jj = ii+1:length(copy_pxs)
                    if do_threshold == 1 && cutMesh.vertex_scale(copy_indices(jj)) < threshold
                        continue
                    end
                    M(copy_pxs(ii),copy_pxs(jj)) = 1;
                    M(copy_pxs(jj),copy_pxs(ii)) = 1;
                end
            end
        else
            for ii=1:length(copy_pxs)-1
                if do_threshold == 1 && cutMesh.vertex_scale(copy_indices(ii)) < threshold
                    continue
                end
                for jj = ii+1:length(copy_pxs)
                    if do_threshold == 1 && cutMesh.vertex_scale(copy_indices(jj)) < threshold
                        continue
                    end
                    cpi = copy_pxs(ii);
                    cpj = copy_pxs(jj);
                    M(cpi, cpj) = max(cutMesh.vertex_scale(copy_indices(jj)), M(cpi,cpj));
                    M(cpj, cpi) = max(cutMesh.vertex_scale(copy_indices(ii)), M(cpj,cpi));
                end
            end
        end
    end
end