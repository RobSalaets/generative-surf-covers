clear
close all

for ii = 1:25
    load(strcat('h',num2str(ii)))
    
    figure
    imagesc(pushed_function*10)
end