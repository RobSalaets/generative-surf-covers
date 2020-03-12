function [G,W,I,ds_cone_to_generic_vertex] = subdivide_and_cut2(obj, genertic_vertex)
    %%
    %== takes the mesh and cuts it to a disc topology along a star whos center
    %== is a specified generic vertex and whos enpoints are the
    %== cones. If no such paths exists, locally subdivide the original mesh. 
    %== we cut in stages. At step i there will be i-1 copies of the
    %== generic vertex on the boundary of the cut mesh. We cut along a shortest path from the
    %== i'th cone to the first copy of the genric vertex in the
    %== boundary 

    % at first set the divided mesh to be equal to the original mesh
    obj.dividedTs2Ts = (1:length(obj.T_orig))';
    obj.V_divided_mesh = obj.V_orig;
    obj.T_divided_mesh = obj.T_orig;
    while true
        %initialize the array holding the distances between the
        %cones and the genric point
        ds_cone_to_generic_vertex = zeros(length(obj.cones),1);
        %create a local copy of cones indices
        cones = obj.cones;
        %update the adjacenct matrix of the divided mesh
        AA = obj.A;
        [Ai, Aj] = find(obj.A);
        AAdist = vecnorm((obj.V_divided_mesh(Ai,:)-obj.V_divided_mesh(Aj,:))');
        AAweights = sparse(Ai,Aj,AAdist);
        % find a path from the generic point to the first cone that
        % intersect none of the other cones
        cones_t = setdiff(cones,cones(1));
        AA(cones_t,:) = 0;
        AA(:,cones_t) = 0;
        for cti = 1:length(cones_t)
            AAweights(cones_t(cti), Aj(Ai==cones_t(cti))) = inf;
            AAweights(Ai(Aj==cones_t(cti)),cones_t(cti)) = inf;
        end
%         [ds_cone_to_generic_vertex(1),p] = graphshortestpath(AA,cones(1),genertic_vertex);
        [~,p] = graphshortestpath(AAweights,cones(1),genertic_vertex);
        ds_cone_to_generic_vertex(1) = length(p)-1;
        % get the edges in the path
        E = [p(1:end-1)' p(2:end)'];
        % cut along the first path
        [G,I] = cut_edges(obj.T_divided_mesh,E);
        I1=I;
        %get the indices of the cones and the generic vertex in the
        %mesh after cutting once
        [~,IA] = unique(I1);
        cones = IA(cones);
        none_cones = IA(genertic_vertex);
        W = obj.V_divided_mesh(I,:);
        %try cutting cones 2 throgh k
        for i = 2:length(obj.cones)
            %locate the correct copy of the none cone point on the
            %bounday. We always take the first copy
            t = triangulation(G,W);
            % get f, the boundary of the cut mesh
            f = t.freeBoundary();
            f = f(:,1);
            % find where the first cone is on the buondary
            ind1 = find(f == cones(1));
            % circle f so that the first cone will be the first
            % index on the boundary of the cut mesh
            f = circshift(f,1-ind1);
            %find the first copy of the generic point
            non_cone = f(1+ds_cone_to_generic_vertex(1));
            % try and find a path from the current cone to the
            % first copy of the generic vertex that does not
            % intersect the boundary of the cut mesh or any of the
            % other cones
            ff = setdiff(f,non_cone);
            AA = adjacency_matrix(G);   
            [Ai,Aj] = find(adjacency_matrix(G));
            AAdist = vecnorm((W(Ai,:)-W(Aj,:))');
            AAweights = sparse(Ai,Aj,AAdist);
            AA(ff,:) = 0;
            AA(:,ff) = 0;
            for cti = 1:length(ff)
                AAweights(ff(cti), Aj(Ai==ff(cti))) = inf;
                AAweights(Ai(Aj==ff(cti)),ff(cti)) = inf;
            end
            cones_t = setdiff(cones,cones(i));
            AA(cones_t,:) = 0;
            AA(:,cones_t) = 0;   
            for cti = 1:length(cones_t)
                AAweights(cones_t(cti), Aj(Ai==cones_t(cti))) = inf;
                AAweights(Ai(Aj==cones_t(cti)),cones_t(cti)) = inf;
            end
            [~,p] = graphshortestpath(AAweights,cones(i),non_cone);
            ds_cone_to_generic_vertex(i) = length(p)-1;
            %if we could find no path then subdivide the mesh and
            %try again
            if ds_cone_to_generic_vertex(i) == -1
                %subdivide the mesh
                disp 'try different order of the cones'
                [obj.V_divided_mesh, obj.T_divided_mesh, cutE, J] = subdivide_mesh_along_line(obj,G,W,non_cone,cones(i),I);
                % update the map from the divided triangles to the
                % undivided triangles
                obj.dividedTs2Ts = obj.dividedTs2Ts(J);
                % add a cell containing the edges cut in the
                % current cutting attempt
                obj.divided_edges{end+1} = cutE;
                %update the adjecency matrix
                obj.A = adjacency_matrix(obj.T_divided_mesh);
                break
            end
            %if we could find a path, get the edges and continue
            %cutting
            E = [p(1:end-1)' p(2:end)'];
            [G,I1] = cut_edges(G,E);
            [~,IA] = unique(I1);
            cones = IA(cones);
            %%find the indices of the copies of the generic vertex
            %%in the boundary of the cut mesh
            none_cones_new = [];
            for j =1:length(none_cones)
                fin = find(I1==none_cones(j));
                none_cones_new = [none_cones_new; fin];
            end
            none_cones = none_cones_new;
            %update I the map from cut vertices to uncut vertices
            I = I(I1);
            W = W(I1,:);
        end
        % if at the last cut there exists a path from the correct
        % copy of the generic point to the last point, then we are
        % done.
        if i==length(obj.cones) && ds_cone_to_generic_vertex(i)~=-1
            break;
        end
    end
end