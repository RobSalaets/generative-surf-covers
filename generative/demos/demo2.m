clear
close all

load('landmark_plane_coords')
for ii = 1:5
    load(num2str(ii))
    if ii == 1
        figure
        pcshow(reshape(pushed_function, [100*100 3]))
    end
    
    
    
    figure
    imagesc(pushed_function*10)
    hold on
    for jj=1:5
    scatter(lm2D(ii,jj,:,1)*100, lm2D(ii,jj,:,2)*100)
    end
end