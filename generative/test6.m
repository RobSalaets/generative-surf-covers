% clear
load('hstats10k.mat')
% 
figure
subplot(1,4,1)
imagesc(pf_avg + 0.5)
title('Gemiddelde')
axis off
for i = 1:3
    
    load(sprintf('bh (%d)', randi(10000)))
    
    subplot(1,4,i+1)
    maxx =max(abs(pushed_function(:)-pf_avg(:)));
    imagesc(((pushed_function-pf_avg) + maxx)/(2*maxx))
    title(sprintf('Variatie %i', i))
    axis off
end


