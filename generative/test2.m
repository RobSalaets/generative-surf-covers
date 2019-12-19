

clear
close all
load('evalues.mat')
load('evectors.mat')
load('meanShape.mat')
load('facesShapeModel.mat')
eigv100 = evalues(:);
eigvec100 = evectors(:,:);
eigvec100 = reshape(eigvec100, [4300 6449 3]);


fig = uifigure('Position',[100 100 350 275]);
tr = trisurf(triangulation(faces, points));
axis vis3d
axis equal
ev = test3(1);

sld = uislider(fig,...
    'Position',[100 75 120 3],...
    'ValueChangedFcn',@(sld,event) updateP(sld,faces,points,eigvec100,eigv100,ev));
sld2 = uislider(fig,...
    'Position',[100 50 120 3],...
    'ValueChangedFcn',@(sld,event) updateVar(sld, ev));

function updateP(sld, faces,points,eigvec100,eigv100,ev)
    scalef =0.01*1e3;
    v = ev.ev;
    trisurf(triangulation(faces, points+ (sld.Value-50)/100*sqrt(eigv100(v))*scalef*squeeze(eigvec100(v,:,:))));
    axis vis3d
    axis equal
    drawnow limitrate
end

function updateVar(sld, ev)
    ev.ev = round(sld.Value)
end