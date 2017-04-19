%% Camera Calibration using 3D calibration object
% Zheng Qi
%% Draw the image points

real_points = [2 2 2;
-2 2 2;
-2 2 -2;
2 2 -2;
2 -2 2;
-2 -2 2;
-2 -2 -2;
2 -2 -2];
figure(1)
plot3(real_points(:,1),real_points(:,2),real_points(:,3),'o')
grid on

pixel_points = [422 323; 178 323; 118 483; 482 483; 438 73; 162 73; 78 117; 522 117];
figure(2)
scatter(pixel_points(:,1),pixel_points(:,2),'o')
grid on

%% Generate matrix P
[x, y] = size(real_points);
P=[];
for i=1:x
    %%
    % In each iteration add 2 rows.
    P=[P; gen_mat_P(real_points(i,:),pixel_points(i,:))];
end
P

%% Calculate M by SVD
[U,S,V] = svd(P);
%%
% Got M matrix from last col of V
M=[(V(1:4,12))';(V(5:8,12))';(V(9:12,12))']

%% Calculate camera center by SVD
[m_U, m_S, m_V] = svd(M);
%%
% t is the last col of V of SVD of M,
% center got from rescaling the t with setting the last one of t to 1
t=m_V(:,end);
center=t(1:end-1)/t(end)

%% Calculate M'
M_prime = M(:,1:3)/ M(3,3)

%% Calculate R_x
cos_x=M_prime(3,3)/sqrt(M_prime(3,2)^2 + M_prime(3,3)^2);
sin_x=-M_prime(3,2)/sqrt(M_prime(3,2)^2 + M_prime(3,3)^2);
R_x=[1, 0, 0;
    0, cos_x, -sin_x;
    0, sin_x, cos_x;]
%%
% The asin() in matlab is arsin, and the result is radian.
% We need to change it to degree.
theta_x = asin(sin_x)*180/pi
N=M_prime*R_x

%% Calculate R_z
cos_z = N(2,2)/sqrt(N(2,1)^2 + N(2,2)^2);
sin_z = -N(2,1)/sqrt(N(2,1)^2 + N(2,2)^2);
R_z = [cos_z, -sin_z, 0;
        sin_z, cos_z, 0;
        0, 0, 1]
theta_z = asin(sin_z)*180/pi

%% Calculate K by M'/(R_x*R_y)
% R_y is a diagnal matrix with the values are 1 on the diagnal.
% R_z is factorized out
K = M_prime/(R_x*diag([1,1,1]));
K = K/K(3,3)
fcoal=[K(1,1),K(2,2)]
img_center=[K(1,3),K(2,3)]
