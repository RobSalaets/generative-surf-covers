function quickscatter2d(V, m, c, varargin)
    size = 5;
    if ~isempty(varargin)
        size=varargin{1};
    end
    if m == 1
        scatter(mod(V(:,1),1), mod(V(:,2),1), size, c)
    else
        scatter(V(:,1),V(:,2), size,c)
    end
end