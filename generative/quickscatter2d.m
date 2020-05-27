function quickscatter2d(V, m, c, varargin)
    if m == 1
        scatter(mod(V(:,1),1), mod(V(:,2),1), 5, c)
    else
        scatter(V(:,1),V(:,2), 5,c)
    end
end