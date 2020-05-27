function [cm, gluer, flattener, varargout] = get_consistent_mapping(tuple, V, T, cones, min_agd_point, ivertex, varargin)
    
    gluer = Gluer(V,T,cones,tuple,min_agd_point);
    flattenerO = Torus_Flattener(gluer.V_torus,gluer.T_torus);
    cmO = CutMesh(gluer, flattenerO, 0);
    
    maps = sort(cmO.inds_mesh_divided_to_inds_plane{ivertex});
    %pick first occurence of ivert in plane as intersection point
    [loops{1}, lps{1}] = trace_loops_flat_torus(flattenerO, gluer, maps(1), 1,[]);
    [loops{2}, lps{2}] = trace_loops_flat_torus(flattenerO, gluer, maps(1), 2, loops{1});
    w1 = normr(gluer.V_torus(loops{1}(3),:)-gluer.V_torus(loops{1}(end-2),:));
    w2 = normr(gluer.V_torus(loops{2}(3),:)-gluer.V_torus(loops{2}(end-2),:));
    v1 = w1;
    v2 = w2;
    if ~isempty(varargin)
        v1 = varargin{1};
        v2 = varargin{2};
    end
    if abs(v1*w2') > abs(v1*w1')
        tmp = loops{1};
        loops{1} = loops{2};
        loops{2} = tmp;
        tmp = lps{1};
        lps{1} = lps{2};
        lps{2} = tmp;
    end
    if length(varargin) > 2
        lms_plane = varargin{3};
        for ii = 1:length(lms_plane)
            c = cones(lms_plane{ii}{2});
            lms_plane{ii}{3} = find(gluer.torus_to_sphere == c);
        end
        flattener = Torus_Flattener(gluer.V_torus,gluer.T_torus, loops, v1, v2, lms_plane);
    else
        flattener = Torus_Flattener(gluer.V_torus,gluer.T_torus, loops, v1, v2);
    end
    cm = CutMesh(gluer, flattener, 1);
    if nargout > 3
        varargout{1} = v1;
        varargout{2} = v2;
    end
    if nargout > 5
        varargout{3} = {};
        for ii=1:length(cones)
            icmV = cm.inds_mesh_divided_to_inds_plane_unsorted{cones(ii)};
            varargout{3}{ii} = {cm.V(icmV,:), ii};
        end
    end
end