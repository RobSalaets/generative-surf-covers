clear
close all
load('evalues.mat')
load('evectors.mat')
load('meanShape.mat')
load('facesShapeModel.mat')
load('human_lms.mat')


tuple1 = {{[6, 7, 8, 9, 10],[1],[2],[3],[4],[5]},{[1, 7, 3, 4, 9],[2],[5],[6],[8],[10]...
    },{[1, 8, 4, 3, 7],[2],[5],[6],[9],[10]},{[2, 5, 4, 7, 6],[1],[3],[8],[9],[10]},{[2, 10, 9, 4, 5],...
    [1],[3],[6],[7],[8]}};
scalef =0.0045*5e2;
minAGD = 1110;
normF = 1000;
V = points / normF;
V = V - mean(V, 1);
inds = [inds(1) inds(3) inds(5) inds(4) inds(2)];
evectors = reshape(evectors, [4300 6449 3]);

ivert = inds(1);
[cm, gluer, flattenerO, v1, v2] = get_consistent_mapping(tuple1, V, faces, inds, minAGD, ivert);

figure(1)
for kk=1:32
    load(sprintf('C:\\Users\\Rob\\Desktop\\Thesis\\Geometry\\progressive_growing_of_gans\\generated\\512-fmb2048-zm-3993\\Gs%i.mat', kk-1))
    load('hstats10k.mat')

    recon = reconstruct_function(cm, pushed_function+pf_avg);

    visualize_with_lms(faces, recon{3}, [], {}, 0,30,0)
    title('Recon Avg')
    drawnow
    pause(0.1)
    frame = getframe(1);
    im = frame2im(frame);
    [imind,qsdf] = rgb2ind(im,256);
    if kk == 1
      imwrite(imind,qsdf,'testgif2.gif','gif', 'Loopcount',inf);
    else
      imwrite(imind,qsdf,'testgif2.gif','gif','WriteMode','append');
    end
end