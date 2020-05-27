classdef Torus_Flattener < handle
    properties  
        % the locations of the flattened vertices on the plane
        V_plane; 
        % the vertices after being cut
        V_cut_torus;  
        % the triangles after being cut
        T_cut_torus;
        % triangulation object of the cut torus
        tri_cut_torus; 
        %  array where I(i) is the index of the original on the torus
        % of vertex i on the cut torus
        I_cut_to_uncut;
        % a 2 by 1 cell array where cut_first_circle{i} is an array of the
        % indices of i'th copy of the vertices in the homology cycle in the cut torus.
        cut_first_circle;
        % a 2 by 1 cell array where cut_second_circle{i} is an array of the
        % indices of i'th copy of the vertices in the homology cycle in the cut torus.
        cut_second_circle;
        % the index of the vertex where the homology cycles intersect
        vertex;
        % Homology cycles
        hb;
        % the four copies of vertex in the cut torus
        cut_vertex;
        % Laplacian matrix
        L;
        L2;
        % the vector such that L*V_plane = b
        b;
        %the affine transforatmion of each triangle
        As;
        %determinant of linear trans of each triangle
        dets;
        %frobenius norm of linear trans of each tri
        frobenius;
        %smaller singular value of linear trans of each tri
        smin;
        %larger singular value of linear trans of each tri
        smax;
        %condition_numbers of linear trans of each tri
        condition_numbers
        %area of each tri
        areas;
        %location of second vertex of parallelogram, the first being [1 0]
        axy;
        % the adjacency matrix of the torus
        adjacency_matrix_of_torus;
    end
    
    methods
        function obj=Torus_Flattener(V,T,varargin)
            % hb is a 2 by 1 cell array, where hb{i} is the indices of
            %the vertices in the it'h homology cycle
            if isempty(varargin)
                hb = homology_of_torus(T,V);
            else
                hb = varargin{1};
            end
            % make sure the cycle intersect at a single vertex
            assert(length(intersect(hb{1},hb{2}))==1,'The homology cycles intersect in an inappropriate way.');
            obj.vertex = intersect(hb{1},hb{2});
            obj.hb = hb;
            % create the edges of the homology cycles
            first_circle = [hb{1}(1:end -1) hb{1}(2:end)];
            second_circle = [hb{2}(1:end -1) hb{2}(2:end)];
            % creat E, the edges we cut the torus along
            E = [first_circle; second_circle];
            %cut the torus to disc
            [obj.T_cut_torus,obj.I_cut_to_uncut] = cut_edges(T,E);
            obj.V_cut_torus = V(obj.I_cut_to_uncut,:);
            % create the triangulation object and adjacency matrix of the cut torus
            obj.tri_cut_torus = triangulation(obj.T_cut_torus,obj.V_cut_torus);
            obj.adjacency_matrix_of_torus = adjacency_matrix(T);
            % creat the boundary of the cut torus
            f = obj.tri_cut_torus.freeBoundary;
            f = f(:,1);
            % Flip the fb for consistency
            f = flip_free_boundary(obj, f);
            % Shift the indices of f for consitency
%             visualize_with_lms(T, V, [], {hb{1}, hb{2}})
            if length(varargin) > 1
                [f, cut_vertex_index] = shift_to_consistent_gradient(obj,f, varargin{2}, varargin{3});
            else
                % find where the cut vertex is on the boundary of the cut
                % torus
                cut_vertex_index = find(obj.I_cut_to_uncut(f)== obj.vertex);
                assert(length(cut_vertex_index) == 4);
                % circle the boundary so the it begins with the cut vertex
                f = circshift(f,1-cut_vertex_index(1));
                cut_vertex_index = cut_vertex_index- cut_vertex_index(1) +1;
            end
            %find the indices of the copies of the two homology cucles on
            %the cut torus
            find_cycels_on_cut_torus(obj,f, cut_vertex_index);
            % make sure the two copies of each homology cycle come from the
            % same homology cycle in the uncut torus.
            assert(all(obj.I_cut_to_uncut(obj.cut_first_circle{1}) == obj.I_cut_to_uncut(obj.cut_first_circle{2})));
            assert(all(obj.I_cut_to_uncut(obj.cut_second_circle{1}) == obj.I_cut_to_uncut(obj.cut_second_circle{2})));
            %create the laplacian matrix
            obj.L = cotmatrix(obj.V_cut_torus,obj.T_cut_torus);
            %clamp negative weights
            clamp(obj);
            obj.L2 = obj.L;
            % create the toric prbifold boundary conditions
            set_boundary_conditions(obj);
            obj.V_plane = solve_tutte(obj);
            if length(varargin) > 3
                obj.V_plane = solve_tutte_with_lms(obj, varargin{4});
            end
            % create the flattening
            min_conformal_distorion_unit_squre(obj);
        end
        function f2 = flip_free_boundary(obj, f)
            %flip f if not counter clockwise
            [f1i, f1j] = find(obj.T_cut_torus == f(1));
            T1 = obj.T_cut_torus(f1i,:);
            [f2i, f2j] = find(T1 == f(2));
            a = [f2j f1j(f2i)];
            if isequal(a,[1 2]) || isequal(a,[2 3]) || isequal(a,[3 1])
                f2 = flip(f);
                disp 'Flipped the free boundary'
            else
                f2 = f;
            end
        end
        function [f2, index2] = shift_to_consistent_gradient(obj,f,v1,v2)
            index = find(obj.I_cut_to_uncut(f)== obj.vertex);
            assert(length(index) == 4);
            c = -1;
            for ii=1:4
                i1 = mod(index(ii)+2, length(f)); i1(i1 == 0) = length(f);
                w1 = normr(obj.V_cut_torus(f(i1),:)-obj.V_cut_torus(f(index(ii)),:));
                i2 = mod(index(ii)-2, length(f)); i2(i2 == 0) = length(f);
                w2 = normr(obj.V_cut_torus(f(index(ii)),:)-obj.V_cut_torus(f(i2),:));
                dots = [v1*w1'; v2*w1'; v1*w2'; v2*w2'];
                [~, m1] = max(abs(dots(1:2)));
                [~, m2] = max(abs(dots(3:4)));
                dot1 = dots(m1); dot2 = dots(2+m2);
%                 assert(abs(dot1) > 0.5 && abs(dot2) > 0.5, 'gradients do not align');
                dot1 = round(dot1 / abs(dot1));
                dot2 = round(dot2 / abs(dot2));
                if dot1 == 1 && dot2 == 1
                    c = ii;
                end
            end
            assert(c > 0, 'No consistent index found');
            f2 = circshift(f,1-index(c));
            index2 = index- index(c) +1;
            index2(index2 < 1) = index2(index2 < 1) + length(f);
            index2 = circshift(index2, 1-c);
        end
        function find_cycels_on_cut_torus(obj,f, index)
            % find the copies of the homology cycles in the boundary of the
            % cut torus.
            obj.cut_first_circle = cell(1,2);
            % the two copies of the first cycle on the cut torus
            % a copy starting from the second appearence of the cut vertex
            % and ending in the third appearence
            obj.cut_first_circle{1} = f(index(2)+1:index(3)-1);
            % a copy starting from the end of the boundary and going backwards
            % until the fourth appearence of the cut vertex
            obj.cut_first_circle{2} = f(end: -1:index(4)+1);
            obj.cut_second_circle = cell(1,2);
            % a copy starting from the first appearence of the cut vertex
            % and ending in the second appearence
            obj.cut_second_circle{1} = f(index(1)+1:index(2)-1);
            % a copy starting from the fourth apperence of the cut vertex
            % and going backwards until the third appearence of the cut vertex
            obj.cut_second_circle{2} = f(index(4)-1:-1:index(3)+1);
            % create the four copies of the intersection of the hoology
            % cycle in the cut torus
            obj.cut_vertex = f(index);
        end
        function clamp(obj)
            m=min(min(tril(obj.L,-1)));
            if m<0
                warning('Mesh is not Delaunay!!');
%                 clamp=1e-2;
                clamp=0.0001;
                fprintf('clamping negative weights to %f\n',clamp);
                obj.L(obj.L<0)=clamp;
                % now fix the laplacian
                inds=sub2ind(size(obj.L),1:length(obj.L),1:length(obj.L));
                obj.L(inds)=0;
                l = linspace(1,length(obj.L),length(obj.L));
                obj.L = obj.L + sparse(l,l,-sum(obj.L(l,:)));
            end
        end
        function set_boundary_conditions(obj)
            
            % create b, the vector such that L*u=b
            obj.b = sparse(length(obj.L),2);
            %fix the four copies of the interesection of the homology cycle
            %to the vertices of unit square.
            obj.L(obj.cut_vertex,:) = 0;
            % set L(cut_vertex(j),i) to be 1 in delta(i,cut_vertex(j))
            obj.L = obj.L  + sparse(obj.cut_vertex,obj.cut_vertex,ones(4,1),...
                length(obj.L),length(obj.L));
            obj.b(obj.cut_vertex,:) = [ 0 0 ;1 0; 1 1 ; 0 1];
            %set every vertex in the first copy of each cycle to be a
            %convex combination of both its neighbors and the neighbors of
            %its "twin" on the second copy.
            set_ccm_proprty_on_cycle(obj,obj.cut_first_circle,1);
            set_ccm_proprty_on_cycle(obj,obj.cut_second_circle,2);
            %set every vertex in the second copy of each cycle to be a
            %a translation of its twin in the first copy
            set_cycles_idendical_up_to_translation(obj,obj.cut_first_circle,1);
            set_cycles_idendical_up_to_translation(obj,obj.cut_second_circle,2);

        end
        function u= solve_tutte_with_lms(obj, lms_plane)
            obj.L = obj.L2;
            w = obj.vertexScale();
            % create b, the vector such that L*u=b
            obj.b = sparse(length(obj.L),2);
            
            %set every vertex in the first copy of each cycle to be a
            %convex combination of both its neighbors and the neighbors of
            %its "twin" on the second copy.
            set_ccm_proprty_on_cycle(obj,obj.cut_first_circle,1);
            set_ccm_proprty_on_cycle(obj,obj.cut_second_circle,2);
            %set every vertex in the second copy of each cycle to be a
            %a translation of its twin in the first copy
            set_cycles_idendical_up_to_translation(obj,obj.cut_first_circle,1);
            set_cycles_idendical_up_to_translation(obj,obj.cut_second_circle,2);
            
            fixed_vars_idx = [];
            V_fixed = [];
            for ii = 1:length(lms_plane)
                torus_idx = lms_plane{ii}{3}; % torus ids zitten in dezelfde kopie en volgorde
                Vlms= lms_plane{ii}{1};
                if ii == 1
                    Vlms = Vlms(3:7,:);
                end
                [ft_idx, ~] = find(obj.I_cut_to_uncut == torus_idx');
                ft_idx = sort(ft_idx);
                [ft_idx, iD] = setdiff(ft_idx, obj.cut_vertex);
                
                fixed_vars_idx = [fixed_vars_idx; ft_idx];
                V_fixed = [V_fixed; Vlms];
            end
            selection = vecnorm((obj.V_plane(fixed_vars_idx,:) - V_fixed)') < 0.4;
            
            V_fixed = V_fixed(selection,:);
            fixed_vars_idx = fixed_vars_idx(selection);
            
            V_fixed = [V_fixed; [ 0 0 ;1 0; 1 1 ; 0 1]];
            fixed_vars_idx = [fixed_vars_idx; obj.cut_vertex];
            

            free_vars_idx = setdiff(1:length(obj.L), fixed_vars_idx);
            obj.b = obj.b - obj.L(:,fixed_vars_idx) * V_fixed;
            
            u =  zeros(length(obj.L),2);
            %%
            exp = -0.5;
            res = 1;
            ii = 0;
            while res > 1e-2 && ii < 5
            wmapped = real(abs(w).^(exp));
            
            W = sparse(1:length(w), 1:length(w), wmapped, length(w), length(w));
            
            u(free_vars_idx,1) = (W*obj.L(:,free_vars_idx))\(W*obj.b(:,1));
            u(free_vars_idx,2) = (W*obj.L(:,free_vars_idx))\(W*obj.b(:,2));
            u(fixed_vars_idx,1) = V_fixed(:,1);
            u(fixed_vars_idx,2) = V_fixed(:,2);
            
            residu_w = vecnorm(W*obj.L(:,free_vars_idx)*u(free_vars_idx,1)-W*obj.b(:,1));
            res = vecnorm(obj.L(:,free_vars_idx)*u(free_vars_idx,1)-obj.b(:,1));
            exp = exp*0.7;
            ii = ii + 1;
            end
            disp(sprintf('Took %i iterations to find weights WLS', ii))
        end
        function set_cycles_idendical_up_to_translation(obj,cycle,cycle_number)
            % set the copies of the first cycle identical up to translation
            % by [0 1]. The copies of the first cycle identical up to translation
            % by [1 0].
            obj.L(cycle{2},:) = 0;
            I = [cycle{2} ; cycle{2}];
            J = [cycle{2} ; cycle{1}];
            obj.L = obj.L + sparse(I,J,[ones(length(I)/2,1);-ones(length(I)/2,1)],length(obj.L),length(obj.L));
            if cycle_number == 1
                obj.b(cycle{2},1) = -1;
            else
                obj.b(cycle{2},2) = 1;
            end
                
        end
        function set_ccm_proprty_on_cycle(obj,cycle,cycle_number)
            % input: a 2 by cell array where each cell is the indices of a
            % copy of cycle on the cut torus.
            % the function sets obj.L and obj.b such that the vertices on
            % the cycle will satisfy the convex combination property
            I =cycle{1};
            J=cycle{2};
            %for every vertex i in the first copy of the cycle add the rwo
            %j, corresponding to its "twin", to row i
            obj.L(I,:) = obj.L(I,:)+obj.L(J,:) ;
            % now L(i,i) appears once where it should appear twice and
            % L(i,j) appears once where it should not appear at all.
            ind = sub2ind(size(obj.L),J,J);
            obj.L = obj.L + sparse(I,I,obj.L(ind),length(obj.L),length(obj.L));
            obj.L = obj.L - sparse(I,J,obj.L(ind),length(obj.L),length(obj.L));
            %set b in the i'th vertices to be the sum of the weight's in
            %the j'th row, where j is the twin of i
            ind = sub2ind(size(obj.L),J,J);
            %in the first cycle we  set the x coordinate 
            if cycle_number == 1
                obj.b(I,1) = obj.L(ind);
            % in the second cycle we set the y coordinate with minus sign
            else
                obj.b(I,2) = -obj.L(ind);
            end
            
        end
        function u = solve_tutte(obj)
            u =  zeros(length(obj.L),2);
            u(:,1) = obj.L\obj.b(:,1);
            u(:,2) = obj.L\obj.b(:,2);
        end
        function [ax, ay] = min_conform_dist(obj,u)
            %takes u that is a flattening of the torus to the square with
            %vertices 0,1 and 1,0 and moves to the parallelogram 1,0 a,b
            % where a,b are such that the conformal energy is minimized.
            E = obj.tri_cut_torus.edges;
            I = E(:,1);
            J = E(:,2);
            x_i = u(I,1);
            y_i = u(I,2);
            x_j = u(J,1);
            y_j = u(J,2);
            ind = sub2ind(size(obj.L),I,J);
            WW = 0.5*full(obj.L(ind));
            %the conformal energy is sum(WW(ind).*(x_i-x_j).^2+(y_i-y_j).^2)- area
            %We compurw ax,ay the minimize the energy of V*[1 0; ax by]
            ax = -sum(WW.*(y_i-y_j).*(x_i-x_j))/sum(WW.*((y_i-y_j).^2));
            ay = 1/(2*(sum(WW.*(y_i-y_j).^2)));
            U = [1 0 ; ax ay];
            obj.axy = [ ax ay];
            u = u*U;
            obj.V_plane = u;
        end
        function draw(obj,u)
            figure(10);
            clf;
            patch('faces',obj.T_cut_torus,'vertices',u,'facecolor','red','facealpha',0.1);
            axis equal;
        end
        function computeAs(obj)
            % calculate map between 2d vertices to matrix mapping each
            % triangle induced by the placement of 2d vertices.
            [V2A, obj.areas] = getFlatteningDiffCoefMatrix(obj.V_cut_torus,obj.T_cut_torus);
            if ispc || ismac
                disp('using mex file. the files under mex folder need to be compiled with mex');
                obj.As = reshape(V2A*obj.V_plane(:),2,2,[]);
            else
                obj.As = reshape((V2A*obj.V_plane(:))',2,2,[]);
            end
        end
        function computeDistortion(obj,force)
            if nargin<2
                force=false;
            end
            if ~force && ~isempty(obj.dets)
                return;
            end
            if isempty(obj.As)||force
                obj.computeAs();
            end
            fprintf('==== Computing distortion etc. ===\n');
            As=obj.As;
            %the entries of A
            a = squeeze(As(1,1,:))';
            b = squeeze(As(1,2,:))';
            c = squeeze(As(2,1,:))';
            d = squeeze(As(2,2,:))';
            % compute area scale of each face
            obj.dets=(a.*d-b.*c)';
            obj.frobenius=sqrt(a.^2+b.^2+c.^2+d.^2)';
            % compute singular vals of A matrices
            alpha = [a+d;b-c]/2; 
            beta = [a-d;b+c]/2;
            alpha=sqrt(sum(alpha.^2));
            beta=sqrt(sum(beta.^2));
            obj.smax=(alpha+beta)';
            obj.smin=abs((alpha-beta)');
            % compute angle distortion of each face
            obj.condition_numbers = obj.smax./obj.smin;
        end
        function vc=vertexScale(obj)
            %return the conformal distortion per vertex, averaged over the
            %adjacent faces according to their area
            %computing the distortion metrics per face
            if isempty(obj.dets)
                obj.computeDistortion();
            end
            % per face - area* det
            c=obj.areas.*obj.dets;
            %accumulate c per vertex
            vals=repmat(c,1,3);
            %v=accumarray(obj.flat_T(:),vals(:));
            v=accumarray(obj.I_cut_to_uncut(obj.T_cut_torus(:)),vals(:));
            %accumulate areas per vertex
            areas=repmat(obj.areas,1,3);
            %vareas=accumarray(obj.flat_T(:),areas(:));
            vareas=accumarray(obj.I_cut_to_uncut(obj.T_cut_torus(:)),areas(:));
            %average the scale by areas
            vc=v./vareas;
            vc = vc(obj.I_cut_to_uncut);
        end
        function conformal_energy=compute_conformal_energy(obj,u)
            %the conformal energy is sum(WW(ind).*(x_i-x_j).^2+(y_i-y_j).^2)- area
            E = obj.tri_cut_torus.edges;
            E = [E; E(:,[2,1])];
            I = E(:,1);
            J = E(:,2);
            x_i = u(I,1);
            y_i = u(I,2);
            x_j = u(J,1);
            y_j = u(J,2);
            ind = sub2ind(size(obj.L),I,J);
            WW = (0.5)*full(obj.L(ind));
            target_edge_lengths = (x_i-x_j).^2+(y_i-y_j).^2;
            f=WW.*target_edge_lengths;
            conformal_energy=sum(f);
        end
        function min_conformal_distorion_unit_squre(obj)
            % set the V to minimize the conformal energy
            obj.min_conform_dist(obj.V_plane);
            % correct by a PSL2Z tranformation to make the lattice base closer to being a square
            % give the base of our current lattice as input
            v = obj.axy;
            u =  [1 0];
            A = correcting_PSL2Z_matrix(u,v);
            % B is the matrix taking our corrected lattice base to [0 1], [1 0]
            B = inv(A*[u;v]);
            c = [obj.axy(1)+1, obj.axy(2)] * B;
            % (double(c<0) is a heuristic to shift back to unit square if needed)
            obj.V_plane = obj.V_plane * B + double(c<0);
%             obj.axy = [1 0];
        end
    end
end

