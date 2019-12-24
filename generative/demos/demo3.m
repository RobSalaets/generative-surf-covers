clear
close all

for ii = 1:10
    load(strcat(num2str(ii),'p'))
    
    figure
    imagesc(tiled_rot*10)
end