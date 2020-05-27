clear
close all
load('tmp.mat')
V = cmO.V;
T = cmO.T;
sz = 128;
triangle_map = cell(length(cmO.T_orig),1);
for iit=1:length(T)
    if mod(iit,100)==0
        iit
    end
    triad = T(iit,:);
    iito = cmO.Ts_plane_to_Ts_orig_mesh(iit);
    tV = V(triad,:);
    [b] = bounds(tV);
    nV = mod(tV, 1.0);
    [pbs] = snap(b,sz);
    %check if pxl in triangle
    if pbs(1) == pbs(2) && pbs(3) == pbs(4)
        weight = polyshape(nV*sz).area;
        triangle_map{iito} = [triangle_map{iito}; pbs(1) pbs(3) weight];
    else
        range1 = mod(pbs(1):(pbs(2)),sz);
        range2 = mod(pbs(3):(pbs(4)),sz);
        t = polyshape(nV*sz);
        if pbs(2) < pbs(1)
            range1 = [0:pbs(2) pbs(1):(sz-1)];
            dir = -sign(min(tV(:,1)));
            t = [polyshape(tV*sz); polyshape((tV + repmat([dir 0], [3 1]))*sz)];
        end
        if pbs(4) < pbs(3)
            range2 = [0:pbs(4) pbs(3):(sz-1)]; %Niet af
            dir = -sign(min(tV(:,2)));
            t = [polyshape(tV*sz); polyshape((tV + repmat([0 dir], [3 1]))*sz)];
        end
        for ix = range1
            for iy = range2
                sq = polyshape([ix ix ix+1 ix+1], [iy+1 iy iy iy+1]);
                weight = sum(sq.intersect(t).area);
                triangle_map{iito} = [triangle_map{iito}; ix iy weight];
%                 if inside([ix; iy], nV(1,:)*sz, nV(2,:)*sz, nV(3,:)*sz)
%                     triangle_map{ix+1,iy+1} = [triangle_map{ix+1,iy+1} iit];
%                 end
            end
        end
    end
end
%%
load('generative/aggr_matrix/triangle_map_128')
sz = 128;
M = sparse(sz^2,sz^2);
threshold = 0.1;
tic;
for tt = 1:length(cmO.T_orig)
    if mod(tt,100) == 0
        tt/length(cmO.T_orig)
    end
    tdata = triangle_map{tt}; % [px py weight]
    tdata = tdata(tdata(:,3) > threshold,:);
    for ii = 1:size(tdata,1)-1
        ilin = tdata(ii,1) + tdata(ii,2)*sz+1;
        others = ii+1:size(tdata,1);
        jlins = tdata(others,1) + tdata(others,2)*sz + 1;
        ilins = ones(length(jlins),1)*ilin;
        Md = sparse([ilins; jlins], [jlins; ilins], repmat(tdata(ii,3) * tdata(others,3), [2 1]), length(M), length(M));
        M = M + Md;
    end
end
toc

function [bnds] = bounds(tV)
    bnds = [min(tV(:,1));
    max(tV(:,1));
    min(tV(:,2));
    max(tV(:,2))];
end

function [pbnds] = snap(bnds,sz)
     pbnds = mod(floor(bnds*sz),sz);
end

function res=  inside(s, a, b, c)
    as_x = s(1)-a(1);
    as_y = s(2)-a(2);

    s_ab = (b(1)-a(1))*as_y-(b(2)-a(2))*as_x > 0;

    if ((c(1)-a(1))*as_y-(c(2)-a(2))*as_x > 0) == s_ab
        res = 0;
        return
    end

    if ((c(1)-b(1))*(s(2)-b(2))-(c(2)-b(2))*(s(1)-b(1)) > 0) ~= s_ab
        res = 0;
        return
    end

    res = 1;
end