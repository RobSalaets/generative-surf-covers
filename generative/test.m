clear

for ii = 10:20
    figure
    load(strcat('h',num2str(ii)))
    imagesc(pushed_function*10, [min(pushed_function*10,[],'all'), max(pushed_function*10,[],'all')])
end