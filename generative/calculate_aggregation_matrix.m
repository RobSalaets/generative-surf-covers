function [groups] = calculate_aggregation_matrix(avg_pushed_function, divide_factor)
%Calculate the aggregation matrix based on the information in cutMesh.
%This assumes a tiling of 1 is sufficient, i.e. mod(V,1) has no overlap.
    
    apf = avg_pushed_function;
    xX = [min(apf(:,:,1), [], 'all') max(apf(:,:,1), [], 'all')];
    xY = [min(apf(:,:,2), [], 'all') max(apf(:,:,2), [], 'all')];
    xZ = [min(apf(:,:,3), [], 'all') max(apf(:,:,3), [], 'all')];
    range = [xX(2)-xX(1); xY(2)-xY(1); xZ(2)-xZ(1)];
    bSize = min(range) / divide_factor;
    offset = [abs(floor(xX(1)/bSize)) abs(floor(xY(1)/bSize)) abs(floor(xZ(1)/bSize))];
    groups = containers.Map('KeyType', 'double', 'ValueType', 'any');
    ymul = ceil(range(2))+2;
    zmul = ymul*(ceil(range(3))+2);
    for px = 1:size(apf, 1)
        for py = 1:size(apf,2)
            p = apf(px, py,:);
            iloc = floor(p / bSize)+offset;
            ilin = iloc(1) + iloc(2)*ymul + iloc(3)*zmul;
            if groups.isKey(ilin)
                groups(ilin) = [groups(ilin) [px; py]];
            else
                groups(ilin) = [px;py];
            end
        end
    end
%     figure
%     hold on
%     vals = groups.values;
%     p = randperm(groups.Count);
%     for ii = 1:groups.Count
%         pixels = vals{ii};
%         scatter(pixels(1,:), pixels(2,:), 2, ones(length(pixels(1,:)),1)*double(p(ii))/double(groups.Count), 'filled')
%     end
    
end