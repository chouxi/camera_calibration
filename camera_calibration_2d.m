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
%homo_list =[
%    1.7715    0.1609   59.7405;
%    0.0266   -1.6353  418.2215;
%   -0.0000    0.0004    1.0000;
%    2.2890    0.0845  121.0000;
%    0.2942   -1.9913  433.0000;
%    0.0011    0.0003    1.0000;
%    1.1370    0.0866   98.0000; 
%   -0.3046   -1.4381  399.0000;
%   -0.0009    0.0003    1.0000;
%    1.7424    0.5691  118.0000;
%   -0.0306   -0.8194  284.0000;
%   -0.0000    0.0017    1.0000;
%   ];
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

%% Computing intrinsic
[V_U, V_S, V_V] = svd(V);
b=V_V(:,end);
B=[b(1) b(2) b(4);
   b(2) b(3) b(5);
   b(4) b(5) b(6)]

v_0=(B(1,2)*B(1,3)-B(1,1)*B(2,3))/(B(1,1)*B(2,2)-B(1,2)^2);
lamda=B(3,3) - (B(1,3)^2 + v_0*(B(1,2)*B(1,3)-B(1,1)*B(2,3)))/B(1,1);
alpha=sqrt(lamda/B(1,1));
beta=sqrt(lamda*B(1,1)/(B(1,1)*B(2,2)-B(1,2)^2));
gamma=-B(1,2)*(alpha^2)*beta/lamda;
u_0=gamma*v_0/alpha-B(1,3)*(alpha^2)/lamda;
K=[alpha gamma u_0;
   0 beta v_0;
   0 0 1]

%% Computing rotation matrix and t
for i=1:len
	file_name = char(file_list(1,i))
	homo = homo_list(i*3-2:i*3,:);
	lambda = 1/sqrt(sum((inv(K) *homo(:,2)).^2));
	r_1=lambda*inv(K)* homo(:,1);
	r_2=lambda*inv(K) * homo(:,2);
	R =[r_1 r_2 cross(r_1, r_2)]
	t = lambda*inv(K) *homo(:,3)
	R_T_R = R'*R
	%%
	% get new R by SVD
	[R_U, R_S, R_V] = svd(R);
	new_R = R_U*R_V'
	new_R_T_R=new_R'*new_R
end
%% Improving Accuracy
H_list =[];
Acc_R_list ={};
Acc_t_list ={};
for i=1:len
	file_name = char(file_list(1,i))
	homo = homo_list(i*3-2:i*3,:);
	%% Step 1
	real_corners = [];
	ind_i = 0;
	for i=1:10
		ind_j = 0;
		for j=1:8
			real_corners = [real_corners; ind_i, ind_j];
			ind_j = ind_j + 30;
		end
		ind_i = ind_i + 30;
	end
	real_corners = [real_corners'; ones(size(real_corners, 1),1)'];
	p_approx = homo * real_corners;
	p_approx = p_approx ./ repmat(p_approx(3,:),size(p_approx,1), 1);
	figure
	image=imread(file_name);
	imshow(image);
	hold on
	title('Figure 1 : Projected grid corners')
	plot(p_approx(1,:), p_approx(2,:), 'ro');
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
	plot(harris_corner(:,1), harris_corner(:,2), 'r+');
	pause(0.5)
	%%% Step 3
	for j=1:size(p_approx,2)
		n = dist2(harris_corner, p_approx(1:2,j)');
		[min_val row_idx] = min(n);
		row_list(j) = row_idx;
	end
	p_correct=harris_corner(row_list(:),:);
	figure
	image=imread(file_name);
	imshow(image);
	hold on
	title('Figure 3 : grid points')
	plot(p_correct(:,1), p_correct(:,2), 'r+');
	pause(0.5)
	%%% Step 4
	p_correct = [p_correct'; ones(1,size(p_correct,1))];
    H = homography2d(real_corners, p_correct)
    H_list = [H_list; H/H(end,end)];
	%%% Step 6 Sum of Square Error
	image_corner = H*real_corners;
	image_corner = image_corner ./ repmat(image_corner(3,:),size(p_approx,1), 1);
	err_reprojection = sqrt(sum(sum((image_corner-p_correct).^2))) / size(image_corner,2)
end
% Step 5 We get step 6 before 5
Acc_V=[];
for i=1:len
	homo = H_list(i*3-2:i*3,:);
	Acc_V=[Acc_V; get_v(homo)];
end
[V_U, V_S, V_V] = svd(Acc_V);
b=V_V(:,end);
B=[b(1) b(2) b(4);
   b(2) b(3) b(5);
   b(4) b(5) b(6)]

v_0=(B(1,2)*B(1,3)-B(1,1)*B(2,3))/(B(1,1)*B(2,2)-B(1,2)^2);
lamda=B(3,3) - (B(1,3)^2 + v_0*(B(1,2)*B(1,3)-B(1,1)*B(2,3)))/B(1,1);
alpha=sqrt(lamda/B(1,1));
beta=sqrt(lamda*B(1,1)/(B(1,1)*B(2,2)-B(1,2)^2));
gamma=-B(1,2)*(alpha^2)*beta/lamda;
u_0=gamma*v_0/alpha-B(1,3)*(alpha^2)/lamda;
K=[alpha gamma u_0;
   0 beta v_0;
   0 0 1]
for i=1:len
	file_name = char(file_list(1,i))
	homo = H_list(i*3-2:i*3,:);
	lambda = 1/sqrt(sum((inv(K) *homo(:,2)).^2));
	r_1=lambda*inv(K)* homo(:,1);
	r_2=lambda*inv(K) * homo(:,2);
	R =[r_1 r_2 cross(r_1, r_2)]
	Acc_R_list{i} = R;
	t = lambda*inv(K) *homo(:,3)
	Acc_t_list{i} = t;
end

%% A way to find corner automatically
% From Harris Corner detector, we can get many corners including the ones we really want.
% First the pixel of the corner we want is dark, it's filter out many times.
% And the measure of the grid is bigger than the noises.
% Second we already know the distance of these 4 points.
% Although the image is distorted, we can get use a approximate distance.
% Then we can use this approximate of grid distance to find the proper grid corners.
% From these grid corners, we get the four outward corners as the corner we want.

%% Augment Reality Image


% Because last 4 digits of my RUID are 3950, so I used picture 7.png
[clipart space alpha] = imread('clipart/7.png');
[row_clip, col_clip, space]=size(clipart);
% Resize the clipart
scalar = min(270/ col_clip, 210/row_clip);
for i=1:len
	file_name = char(file_list(1,i))
	homo = H_list(i*3-2:i*3,:);
	image = imread(file_name);
	%projection = K*[Acc_R_list{i}(:,1:2) Acc_t_list{i}];
	for j=1:row_clip
		for k=1:col_clip
			if alpha(row_clip+1-j,k) ~= 0
				p = homo*[k*scalar;j*scalar;1];
				p = round(p / p(end,end));
				image(p(2),p(1),:) = clipart(row_clip+1-j,k,:);
			end
		end
	end
	figure
	imshow(image)
	pause(0.5)
end

%% Augment Reality Object

obj_points = [0 0 0;
			  90 0 0;
			  90 90 0;
			  0 90 0;
			  0 0 90;
			  90 0 90;
			  90 90 90;
			  0 90 90;
			  ];
for i=1:len
	file_name = char(file_list(1,i))
	image = imread(file_name);
	projection = K*[Acc_R_list{i} Acc_t_list{i}];
	figure
	imshow(image)
	hold on
    obj_img_pt = zeros(size(obj_points,1),2);
	for j=1:size(obj_points ,1)
		p = projection*[obj_points(j,1);obj_points(j,2); obj_points(j,3);1];
		p = round(p / p(end,end));
        obj_img_pt(j,:)=p(1:2)';
    end
    plot([obj_img_pt(1,1) obj_img_pt(2,1)], [obj_img_pt(1,2) obj_img_pt(2,2)],'r','LineWidth',2)
    plot([obj_img_pt(2,1) obj_img_pt(3,1)], [obj_img_pt(2,2) obj_img_pt(3,2)],'r','LineWidth',2)
    plot([obj_img_pt(3,1) obj_img_pt(4,1)], [obj_img_pt(3,2) obj_img_pt(4,2)],'r','LineWidth',2)
    plot([obj_img_pt(4,1) obj_img_pt(1,1)], [obj_img_pt(4,2) obj_img_pt(1,2)],'r','LineWidth',2)

    plot([obj_img_pt(5,1) obj_img_pt(6,1)], [obj_img_pt(5,2) obj_img_pt(6,2)],'r','LineWidth',2)
    plot([obj_img_pt(6,1) obj_img_pt(7,1)], [obj_img_pt(6,2) obj_img_pt(7,2)],'r','LineWidth',2)
    plot([obj_img_pt(7,1) obj_img_pt(8,1)], [obj_img_pt(7,2) obj_img_pt(8,2)],'r','LineWidth',2)
    plot([obj_img_pt(8,1) obj_img_pt(5,1)], [obj_img_pt(8,2) obj_img_pt(5,2)],'r','LineWidth',2)

    plot([obj_img_pt(1,1) obj_img_pt(5,1)], [obj_img_pt(1,2) obj_img_pt(5,2)],'r','LineWidth',2)
    plot([obj_img_pt(2,1) obj_img_pt(6,1)], [obj_img_pt(2,2) obj_img_pt(6,2)],'r','LineWidth',2)
    plot([obj_img_pt(3,1) obj_img_pt(7,1)], [obj_img_pt(3,2) obj_img_pt(7,2)],'r','LineWidth',2)
    plot([obj_img_pt(4,1) obj_img_pt(8,1)], [obj_img_pt(4,2) obj_img_pt(8,2)],'r','LineWidth',2)
	pause(0.5)
end

%% Extra Credit
% From the chapter given
% If n = 2, we can impose the skewless constraint γ = 0, i.e., [0, 1, 0, 0, 0, 0]b = 0
% That means we assume the pixels construct a rectanguler

extra_credit(homo_list, file_list)
