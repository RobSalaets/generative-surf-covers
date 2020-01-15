load('testgen128.mat')

m = min(pushed_function, [],'all');
M = max(pushed_function, [], 'all');
resc = (pushed_function-m)./(M-m);
figure
imagesc(permute(resc, [2 3 1]))