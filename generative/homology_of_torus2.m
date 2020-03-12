function hb = homology_of_torus2(T,V, idx)
    % hb is a 2 by 1 cell array with each cell containig a cycle
    hb = cell(2,1);
    E = cut_graph2(V,T);
    % make an adjacency matrix from E
    Ea = sparse(E(:,1),E(:,2),ones(size(E,1),1),size(V,1),size(V,1));
    Ea = Ea + Ea';
    %find the first cycle on the cut graph
    %Why does looking for loops make sense? MST will create a disk
    %connectivity because of minimal representation. Genus one will have 2
    %loops. Try to force first loop to pass certain vertex
    hb{1} = find_longest_loop_in_graph(Ea);
    Adj = adjacency_matrix(T);
    %Find closest on loop
    mini = [];
    mindist = inf;
    minp = [];
    for jj = 1:length(idx)
        for ii = 1:round(length(hb{1})/5):length(hb{1})
            AT = Adj;
            others = [hb{1}(1:ii-1) hb{1}(ii+1:end)];
            AT(others,:)=0;
            AT(:,others)=0;
            [~,p] = graphshortestpath(AT, idx(jj), hb{1}(ii));
            if length(find(hb{1} == idx(jj))) > 0
                minp = [];
                mindist = 0;
                mini(:) = find(hb{1}==idx(jj));
                mini(2) = jj;
                break
            end
            if length(p) < mindist && length(p) > 0
                mindist = length(p);
                minp = p;
                mini(1) = ii;
                mini(2) = jj;
            end
        end
        if mindist == 0
            break
        end
    end
    
    Adj(minp(2:end),:)=0;
    Adj(:,minp(2:end))=0;
    range = length(minp)*2;
    for ii=1:range
        AT = Adj;
        loc1 = mod(mini(1)-ii, length(hb{1})-1);
        if loc1 == 0
            loc1 = loc1 + length(hb{1})-1;
        end
        others = [hb{1}(1:loc1-1) hb{1}(loc1+1:end)];
        AT(others,:)=0;
        AT(:,others)=0;
        [~,p1] = graphshortestpath(AT,  hb{1}(loc1),idx(mini(2)));
        if ~isempty(p1) && length(p1) < 1.8*length(minp)
            break
        end
    end
    for ii=1:range
        AT = Adj;
        loc2 = mod(mini(1)+ii, length(hb{1})-1);
        if loc2 == 0
            loc2 = loc2 + length(hb{1})-1;
        end
        others = [hb{1}(1:loc2-1) hb{1}(loc2+1:end)];
        AT(others,:)=0;
        AT(:,others)=0;
        [~,p2] = graphshortestpath(AT, idx(mini(2)), hb{1}(loc2));
        if ~isempty(p2) && length(p2) < 1.8*length(minp)
            break
        end
    end
    assert(~isempty(p1) && ~isempty(p2));
    
    hb{1} = [hb{1}(1:loc1-1) p1(1:end-1) p2(1:end-1) hb{1}(loc2:end)];
    c = [hb{1}(1:end-1); hb{1}(2:end)]';
    hb{1}  = hb{1}';
    %cut along the  first cycle
    [G,I] = cut_edges(T,c);
    W = V(I,:);
    %pick a point on the boundary and compute the shortest path to it's
    %"twin"
    t = triangulation(G,W);
    f = freeBoundary(t);
    f= f(:,1);
    A = adjacency_matrix(G);
    %   for every boundary vertex v if there is a loop from  v to it's twin 
    %   in the boundary of cut cylinder, disjoint from the boundary
    %   that is our second loop
%     for j  = 1: size(f,1)
    half_loop = floor(length(hb{1})/2)-1;
    visualize_with_lms(T, V, [], {hb{1}, minp, p1, p2})
    indices = find(hb{1} == idx(mini(2))) + floor(0.5:0.5:half_loop+.5).*(cos(pi*(0:2*half_loop)));
    for j = indices
        %S are the copies of v in W (the torus cut to a cylinder)
        S = find(I==hb{1}(j));
        % if there is a path not intersecting the boundary from f(j) to its
        % twin then we are done. otherwise keep searching.
        A_temp = A;
        f_temp = setdiff(f,S);
        A_temp(f_temp,:) = 0;
        A_temp(:,f_temp) = 0;
        [~,p]= graphshortestpath(A_temp,S(1),S(2));
        if ~isempty(p)
            break
        end
    end
    hb{2} = I(p);

end