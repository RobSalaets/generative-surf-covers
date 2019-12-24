clear

for ii = 11:50
    ii
    figure
    load(strcat('n',num2str(ii)))
    imagesc((pushed_function + 1)*0.25)
end
