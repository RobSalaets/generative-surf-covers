function E = cut_graph2(V,T)
%%= the output are edges E in the mesh such that (T,V)/E is a disc
    %G is an adjacency matrix of the dual graph with the i'th row of 
    %T being the i'th vertex of the dual graph
    %Add weights wrt idx to have some consistency
    tri  = triangulation(T,V);
    % create the dual graph
    E = tri.edges();
    dual_edges = tri.edgeAttachments(E);  %returns triangles lying l/r from edge, therefore this constructs dual edges, as faces become verts
    dual_edges = cell2mat(dual_edges);
    % create the dual graph adjacency matrix
    A = sparse(dual_edges(:,1),dual_edges(:,2),ones(length(dual_edges),1),length(T),length(T));
    A = A+ A';    
    
    [Ai, Aj] = find(A);
    centroids = tri.incenter();
    minX = min(centroids(:,1)) - 0.1*abs(min(centroids(:,1)));
    maxY = max(abs(centroids(:,2)))*1.1;
    weights = (minX - 0.5*(centroids(Ai,1) + centroids(Aj,1)));% + ... 
              %0.8*(maxY-abs(0.5*(centroids(Ai,1) + centroids(Aj,1))));
%     AAdist = vecnorm((centroids(Ai,:)-centroids(Aj,:))');
    A = sparse(Ai,Aj,weights,length(T),length(T));
%     Ag = graph(A);
%     
%     [idxe, ~] = find(E==idx);
%     fidx = dual_edges(idxe(1), 1);
%     node_dist = Ag.distances(fidx);
%     edge_weight = (node_dist(Ai) + node_dist(Aj))/2;
%     A = sparse(Ai,Aj,1./(edge_weight),length(T),length(T));
    gr = graph(A);    %graph with orig faces as nodes and bidirectional
    
    %Tree is a matrix of the tree in G
    Tree = gr.minspantree;   
    dual_tree_edges = Tree.Edges;
    %convert the edges of the co-spanningtree to a matrix
    dual_tree_edges = table2array(dual_tree_edges);
    dual_tree_edges = dual_tree_edges(:,[1,2]);
    dual_edges = [ dual_edges ; dual_edges(:,[2,1])];
    % convert the coedges to edges
    [~,indx]=ismember(dual_tree_edges(:,:),dual_edges,'rows');
    indx = mod(indx,length(dual_edges)/2);
    indx(indx==0) = length(dual_edges)/2;
    % TE are the edges in the original matrix corresponding to the edges of
    % the minimal spanning tree of the dual graph
    TE = E(indx,:);
    % the cut graph is the edges of the mesh\TE
    E = setdiff(E,TE,'rows');
end