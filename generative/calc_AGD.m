function [AGD_raw, minAGDidx] = calc_AGD(V,F)
    adj = triangulation2adjacency_change(F,V');
    dist = graphallshortestpaths(adj,'directed',false);
    faceAreas = computeSurfAreas(V,F);
    nV = length(V);
    oneRingAreas = zeros(nV,1);
    for ii = 1:nV
        ff = any(F'==ii); % indices of faces in the ii'th vertex 1-ring
        oneRingAreas(ii) = (1/3)*sum(faceAreas(ff));
    end
    AGD_raw = dist*oneRingAreas;
    [~, minAGDidx] = min(AGD_raw);
end