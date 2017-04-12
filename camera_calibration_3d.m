%% Camera Calibration using 3D calibration object
%% Draw the image points

real_points = [2 2 2;
-2 2 2;
-2 2 -2;
2 2 -2;
2 -2 2;
-2 -2 2;
-2 -2 -2;
2 -2 -2]
figure(1)
plot3(real_points(:,1),real_points(:,2),real_points(:,3),'o')
grid on

pixel_points = [422 323; 178 323; 118 483; 482 483; 438 73; 162 73; 78 117; 522 117]
figure(2)
scatter(pixel_points(:,1),pixel_points(:,2),'o')
grid on

%% 