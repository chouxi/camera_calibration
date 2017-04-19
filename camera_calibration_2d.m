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
%homo_list = [];
%for i=1:len
%    file_name = char(file_list(1,i));
%    image = imread(file_name);
%    figure
%    imshow(image);
%    [corners_x, corners_y] = ginput(4);
%    corner_mat = [corners_x'; corners_y'; ones(size(corners_x))'];
%    homo_tmp = homography2d(real_corners, corner_mat);
%    homo_list = [homo_list; homo_tmp/homo_tmp(end,end)];
%end
homo_list =[
    1.7715    0.1609   59.7405;
    0.0266   -1.6353  418.2215;
   -0.0000    0.0004    1.0000;
    2.2890    0.0845  121.0000;
    0.2942   -1.9913  433.0000;
    0.0011    0.0003    1.0000;
    1.1370    0.0866   98.0000;
   -0.3046   -1.4381  399.0000;
   -0.0009    0.0003    1.0000;
    1.7424    0.5691  118.0000;
   -0.0306   -0.8194  284.0000;
   -0.0000    0.0017    1.0000;
   ];
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

%% Computing intrinsic
v_0=(B(1,1)*B(1,3)-B(1,1)*B(2,3))/(B(1,1)*B(2,2)-B(1,2)^2);
lamda=B(3,3) - (B(1,3)^2 + v_0*(B(1,2)*B(1,3)-B(1,1)*B(2,3)))/B(1,1);
alpha=sqrt(lamda/B(1,1));
beta=sqrt(lamda*B(1,1)/(B(1,1)*B(2,2)-B(1,2)^2));
gamma=-B(1,2)*(alpha^2)*beta/lamda;
u_0=gamma*v_0/alpha-B(1,3)*(alpha^2)/lamda;
A=[alpha gamma u_0;
   0 beta v_0;
   0 0 1]

%% Computing rotation matrix and t
for i=1:len
	file_name = char(file_list(1,i))
	homo = homo_list(i*3-2:i*3,:);
	r_1=lamda*inv(A)* homo(:,1);
	r_2=lamda*inv(A) * homo(:,2);
	R =[r_1 r_2 cross(r_1, r_2)]
	t = lamda*inv(A) *homo(:,3)
	R_T_R = R'*R
	%%
	% get new R by SVD
	[R_U, R_S, R_V] = svd(R);
	new_R = R_U*R_V'
	new_R_T_R=new_R'*new_R
end
%% Improving Accuracy
for i=1:len
	file_name = char(file_list(1,i))
	homo = homo_list(i*3-2:i*3,:);
	%% Step 1
	p_approx = homo * real_corners;
	p_approx = p_approx ./ repmat(p_approx(3,:),size(p_approx,1), 1);
	figure
	image=imread(file_name);
	imshow(image);
	hold on
	title('Figure 1 : Projected grid corners')
	plot(p_approx(1,:), p_approx(2,:), 'o');
	% pause for plot complete context in one figure
	pause(0.5)
	%%% Step 2
	[cim, r, c, rsubp, csubp] = harris(rgb2gray(image), 2, 500, 2, 0);
	harris_corner=[csubp, rsubp];
	figure
	image=imread(file_name);
	imshow(image);
	hold on
	title('Figure 2 : Harris corners')
	plot(harris_corner(:,1), harris_corner(:,2), '+');
	pause(0.5)
	%%% Step 3
	for j=1:size(p_approx,2)
		n = dist2(harris_corner, p_approx(1:2,j)');
		[min_val row_idx] = min(n);
		row_list(j) = row_idx;
	end
	row_list
	p_correct=harris_corner(row_list(:),:);
	figure
	image=imread(file_name);
	imshow(image);
	hold on
	title('Figure 3 : grid points')
	plot(p_correct(:,1), p_correct(:,2), '+');
	pause(0.5)
end
