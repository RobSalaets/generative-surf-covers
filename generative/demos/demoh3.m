clear
close all
load('gen0')
pushed_function = permute(pushed_function, [2 3 1]);
figure
pcshow(reshape(pushed_function, [64^2 3]))

load('gen1')
pushed_function = permute(pushed_function, [2 3 1]);
figure
pcshow(reshape(pushed_function, [64^2 3]))

load('gen2')
pushed_function = permute(pushed_function, [2 3 1]);
figure
pcshow(reshape(pushed_function, [64^2 3]))

% NOTE overlap