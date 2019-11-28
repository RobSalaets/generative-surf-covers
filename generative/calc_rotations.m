% clear
close all


landmarks = load('landmarks.mat'); landmarks = landmarks.landmarks;
names_cell = load('names.mat'); names_cell = names_cell.names_cell;
n = length(names_cell);
lms_V = zeros(5,3,n);
for ii = 1:n
    ii
    fname = names_cell{ii};
    
    [V,F] = read_ply(strcat('meshes/all/',fname));
    V(:,3) = -V(:,3);
    [V, F] = delaunayize(V,F);
    m = mean(V, 1);
    V = V-m;
    lms_V(:,:,ii) = V(landmarks(ii,:),:);
end
% lms_N = zeros(5,3,n);
% for ii = 1:n
%     norms = [norm(lms_V(1,:,ii));norm(lms_V(2,:,ii));norm(lms_V(3,:,ii));norm(lms_V(4,:,ii));norm(lms_V(5,:,ii)) ];
%     avg_norm = mean(norms);
%     lms_N(:,:,ii) = lms_V(:,:,ii) / avg_norm;
% end
ref = lms_V(:,:,1)';
R = zeros(3,3, n);
R(:,:,1) = eye(3);
for ii=2:n
    ii
    L = lms_V(:,:,ii)';
    L(2,:) = -L(2,:);
    R(:,:,ii) = calc_rot(L, ref);
end
% save('generative\rotmat.mat', 'R');

figure
hold on
for ii = 1:n
    pcshow((R(:,:,ii)*lms_V(:,:,ii)')',.2:0.2:1.0,'MarkerSize',100);
end

figure
hold on
for ii = 1:n
    pcshow(lms_V(:,:,ii),.2:0.2:1.0,'MarkerSize',100);
end

function [Rr] = calc_rot(L, ref)
    A = zeros(15, 9);
    b = zeros(15, 1);
    c = 1;
    for ii = 1:3
        for jj = 1:5
            A(c,(1+3*(ii-1)):(3+3*(ii-1))) = L(1:3,jj)';
            c=c+1;
        end
    end
    b(1:5,1) = ref(1,:);
    b(6:10,1) = ref(2,:);
    b(11:15,1) = ref(3,:);


    options = optimoptions('fmincon','SpecifyObjectiveGradient',true);
    rot = fmincon(@(x)LS(x,A,b), A\b, [], [], [], [], [], [], @constr, options);

    Rr = [rot(1:3)'; rot(4:6)'; rot(7:9)'];
end

function [f,G] = LS(x,A,b)
    f = 0.5*(A*x-b)'*(A*x-b);
    G = A'*A*x -A'*b;
end

function [c, ceq] = constr(x)
    c = [];
    ceq = zeros(5,1);
    
    
    
    for ii = 1:3
        for jj = 1:3
            v = x((1 + (ii-1)*3):(3 + (ii-1)*3))' * x((1+(jj-1)):3:(7+(jj-1)));
            if ii ~= jj
                ceq(jj + (ii-1)*3) = v;
            end
        end
    end
    v = x(1:3:7)' * x(1:3:7) - x(2:3:8)' * x(2:3:8);
    ceq(1) = v;
    v = x(1:3:7)' * x(1:3:7) - x(3:3:9)' * x(3:3:9);
    ceq(2) = v;
    v = x(1:3:7)' * x(2:3:8);
    ceq(3) = v;
    v = x(1:3:7)' * x(3:3:9);
    ceq(4) = v;
    v = x(3:3:9)' * x(2:3:8);
    ceq(5) = v;
end
