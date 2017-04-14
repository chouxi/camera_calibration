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
homo_list = [];
for i=1:len
    file_name = char(file_list(1,i));
    image = imread(file_name);
    figure
    imshow(image);
    [corners_x, corners_y] = ginput(4);
    corner_mat = [corners_x'; corners_y'; ones(size(corners_x))'];
    homo_tmp = homography2d(real_corners, corner_mat);
    homo_list = [homo_list; homo_tmp/homo_tmp(end,end)];
end
%%
% For printing file_name with homo value
for i=1:len
	file_name = char(file_list(1,i))
	homo = homo_list(i*3-2:i*3,:)
end

%% Computing B matrix
V=[];
for i=1:len
	homo = homo_list(i*3-2:i*3,:);
	V=[V; get_v(homo)];
end

[V_U, V_S, V_V] = svd(V);
b=V_V(:,end);
B=[b(1) b(2) b(4);
   b(2) b(3) b(5);
   b(4) b(5) b(6)]

%% Computing inrinsic
v_0=(B(1,1)*B(1,3)-B(1,1)*B(2,3))/(B(1,1)*B(2,2)-B(1,2)^2);
