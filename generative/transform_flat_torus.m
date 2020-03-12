function tiled_rot = transform_flat_torus(pushed_function, trans, rot)
    %%%%
    % trans: 2dim vector in 0-1 space
    % rot: rotation in deg
    %%%%
    image_size = size(pushed_function,1);
    channels = size(pushed_function, 3);
    new_pushed_function = zeros(image_size,image_size,channels);
    pxl_shift = round(trans * image_size);
    px = pxl_shift(1);
    py = pxl_shift(2);

    if px > 0
        new_pushed_function(:,px+1:end, :) = pushed_function(:,1:end - px,:);
        new_pushed_function(:,1:px, :) = pushed_function(:,end - px+1:end,:);
    else
        px = abs(px);
        new_pushed_function(:,end-px+1:end, :) = pushed_function(:,1:px,:);
        new_pushed_function(:,1:end-px, :) = pushed_function(:,1+px:end,:);
    end
    temp = new_pushed_function;
    if py > 0
        new_pushed_function(py+1:end,:, :) = temp(1:end - py,:,:);
        new_pushed_function(1:py,:, :) = temp(end - py+1:end,:,:);
    else
        py = abs(py);
        new_pushed_function(end-py+1:end,:, :) = temp(1:py,:,:);
        new_pushed_function(1:end-py,:, :) = temp(1+py:end,:,:);
    end

    %Tile for rotate
    tiled = zeros(image_size*2,image_size*2,channels);
    tiled(1:image_size, 1:image_size,:) = new_pushed_function;
    tiled(1:image_size, image_size+1:end,:) = new_pushed_function;
    tiled(image_size+1:end, 1:image_size,:) = new_pushed_function;
    tiled(image_size+1:end, image_size+1:end,:) = new_pushed_function;
    
    tiled_rot = imrotate(tiled, rot, 'bilinear','crop'); %Rotation distorts data!!! points slide along linear surface around point cloud
    tiled_rot = tiled_rot(image_size/2+1:image_size+image_size/2,image_size/2+1:image_size+image_size/2,:);
    
end