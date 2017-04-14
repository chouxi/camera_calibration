%% Camera Calibration using 2D calibration object
% Zheng Qi

%% Extract the four corners of the image
file_list = {'images2.png', 'images9.png', 'images12.png', 'images20.png'};
[r,len] = size(file_list);
%%
% bottom left->bottom right->top left->top right
real_corners = [0 270 0 270;
                0 0 210 210;
                1 1 1 1];
for i=1:len
    file_name = char(file_list(1,i))
    image = imread(file_name);
    figure
    imshow(image);
    [corners_x, corners_y] = ginput(4);
    corner_mat = [corners_x'; corners_y'; ones(size(corners_x))']
    homography2d(real_corners, corner_mat)
end