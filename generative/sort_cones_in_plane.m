function [sorted_cones, varargout] = sort_cones_in_plane(V,T, cones, varargin)
    %%%
    % FRONT SIZE SHOULD BE IN POSITIVE AXIS DIR
    %%%
    conesV = V(cones,:);
    mConesV = mean(conesV, 1);
    conesV = conesV-repmat(mConesV, [length(cones), 1]);
    [~, ~, sV] = svd(conesV'*conesV);
    perp = sV(:,end);
    [~, maxi] = max(abs(perp));
    perp = perp * sign(perp(maxi));
    first = conesV(1,:);
    angles(1) = 0;
    for ii=2:size(conesV, 1)
        d = conesV(ii,:)*first'/(norm(conesV(ii,:))*norm(first));
        angles(ii) = sign(conesV(ii,:)*cross(perp, first)')*acos(d);
    end
    [~, perm] = sort(angles, 'descend');
    sorted_cones = cones(perm);
    
    if ~isempty(varargin)
        [sAGD, I] = sort(varargin{1});
        topnum = round(size(varargin{1},1)/80);
        top = I(1:topnum);
        adj = adjacency_matrix(T);
        unseen = top;
        nb_c = 0;
        while ~isempty(unseen)
            nb_c = nb_c + 1;
            clust{nb_c} = unseen(1);
            %add boys
            possibrel = find(adj(clust{nb_c},:));
            while ~isempty(possibrel)
                clust{nb_c} = [clust{nb_c} intersect(top, possibrel)'];
                [~, possibrel] = find(adj(clust{nb_c},:));
                possibrel = unique(possibrel);
                possibrel = intersect(setdiff(possibrel, clust{nb_c}), top)';
            end
            unseen = setdiff(unseen, clust{nb_c});
        end
        for ii = 1:length(clust)
            cMeans(ii,:) = mean(V(clust{ii},:), 1);
        end
        dp = (cMeans-repmat(mConesV, [nb_c 1]))*perp;
        [~, mdp] = min(dp);
        ids = clust{mdp};
        
        [sAGD, I] = sort(varargin{1}(ids));
        t5 = ids(I(1:5));
        nn = zeros(length(t5),1);
        for i=1:length(t5)
            nn(i) = length(find(adj(t5(i),:)));
        end
        [~, mnn] = max(nn);
        
        varargout{1} = t5(mnn);
        if length(varargin) >= 2
            dir = varargin{2};
            conesV = V(sorted_cones,:);
            mConesV = mean(conesV, 1);
            conesV = conesV-repmat(mConesV, [length(cones), 1]);
            assert(all(size(dir) == [1 3]));
            dps = conesV * dir';
            [~, iloc] = max(dps);
            sorted_cones = circshift(sorted_cones, -(iloc-1));
        end
    end
end