clear
close all
load('landmark_plane_coords.mat')
N = 81;
for ii = 1:N
    points(ii,:,:) = reshape(lm2D(ii, :,:,:), [25 2]);
end
idxs = [1 2 3 6 7 8]; % Indices of center 'green/yellow structure'
Dall = zeros(N,3,2);

for ii = 1:N
    ii
    avg = squeeze(sum(points(ii, idxs, :), 2)./6)';
    for jj=1:6
        if norm(avg-squeeze(points(ii, idxs(jj),:))')>= 0.1
            disp 'do manual select'
            avg45 = [0.0020    0.3032];
            avg = avg45;
            break
        end
    end
    
    D = [];
    D(2,:) = avg;
    %purpple
    [pp1, iloc] = getClosestPointTo(avg, squeeze(points(ii, 4:5:24,:)));
    [pp2, ~] = getClosestPointTo(pp1, squeeze(points(ii, setdiff(4:5:24, 4+(iloc-1)*5),:)));
    D(1,:) = (pp1 + pp2) * 0.5;
    %green
    [pg1, iloc] = getClosestPointTo(avg, squeeze(points(ii, 5:5:25,:)));
    a = D(1,:)-avg;
    b = pg1 - avg;
    pg2pred = pg1 - 2*(b-(a*b')/(a*a')*a);
    [pg2, ~] = getClosestPointTo(pg2pred, squeeze(points(ii, setdiff(5:5:25, iloc*5),:)));
    D(3,:) = (pg1 + pg2) * 0.5;
    Dall(ii, :,:) = D;
    
%     load(num2str(ii))
%     figure
%     image(pushed_function*10)
%     hold on
%     scatter(D(:,1)*100,D(:,2)*100)
end

transl = zeros(2,N);
transl(:,2) = -mean(squeeze(Dall(2,:,:))',2);
rotations = zeros(N, 2,2);
rotations(2,:,:) = eye(2);
angles = zeros(N,1);
% Calc slope for rot alignment
A = [squeeze(Dall(2,:,1))' ones(3,1)]; b = squeeze(Dall(2,:,2))';
rc2 = A \ b;
rc2 = rc2(1);
angle_r = atan2(rc2,1);
angles(2) = angle_r;

% figure
% hold on
for ii = setdiff(1:N,2)
    A = [squeeze(Dall(ii,:,1))' ones(3,1)]; b = squeeze(Dall(ii,:,2))';
    rc = A \ b;
    rc = rc(1);
    if (A(1,1) < A(2,1))
        angle = atan2(rc,1);
    else
        angle = atan2(-rc,-1);
    end
    dangle = angle;
    angles(ii)=dangle;
    rotations(ii,:,:) = [cos(dangle) -sin(dangle);
           sin(dangle) cos(dangle)];
    transl(:,ii) = -mean(squeeze(Dall(ii,:,:))',2);
    
    T = squeeze(rotations(ii,:,:)) * (squeeze(Dall(ii,:,:))' + [transl(:,ii) transl(:,ii) transl(:,ii)]);
    
    x = T(1,:);
    y = T(2,:);
%     axis([min(x)-1, max(x)+1, min(y)-1, max(y)+1])
%     text(x(1),y(1),num2str(ii),'Color','red')
%     text(x(2),y(2),num2str(ii), 'Color','green')
%     text(x(3),y(3),num2str(ii), 'Color','blue')
    
end
save('generative/toricTransforms.mat', 'angles','rotations','transl')

function [p, iloc] = getClosestPointTo(target, points)
    assert(size(target,2)== 2);
    assert(size(points,2) == 2);
    assert(length(size(points))==2);
    pp_closest = zeros(length(points),2);
    for jj=1:length(points)
        c=1;
        %shift and select closest to target
        loc9 = zeros(9,2);
        for s1 = -1:1
            for s2 = -1:1
                shift = [s1 s2];
                loc9(c, :) = points(jj,:) + shift;
                c = c +1;
            end
        end
        diff = loc9 - repmat(target , [9 1]);
        normsq = diff(:,1).^2 + diff(:,2).^2;
        [~, mloc] = min(normsq);
        pp_closest(jj,:) = loc9(mloc, :);
    end
    
    diff = pp_closest - repmat(target , [length(points) 1]);
    normsq = diff(:,1).^2 + diff(:,2).^2;
    [~, iloc] = min(normsq);
    p = pp_closest(iloc,:);
end
