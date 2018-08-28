% MAIN
%% Test the 3D TO 3D model
% test = BTL_ICRA();
% [R_final, t_final, P_MTI_tform, P_Stereo_tform] = test.Calibration_3D23D();

%% Method 1: ANN
test = BTL_ICRA();
flag_method = 2;
[pixel_test] = test.TestDesign(flag_method);
X_laser = pixel_test(:,2); 
Y_laser = pixel_test(:,1);
[net, P_MTI_2d, P_MTI_3d] = test.Calibration_ANN();
data_ANN = net([pixel_test(:,2), pixel_test(:,1)]'); 

%% Method 2: Interpolation
test = BTL_ICRA();
flag_method = 2;
[pixel_test] = test.TestDesign(flag_method);
X_laser = pixel_test(:,2); % height -- 1 
Y_laser = pixel_test(:,1); % width -- 2
[data_Interpolation, Cell_2d, Cell_3d, P_MTI_2d, P_MTI_3d] = test.Calibration_Interpolation(X_laser, Y_laser); 

%%
figure(1);clf
subplot(1,2,1)
pcshow(P_MTI_3d,'b');
hold on;
pcshow(data_Interpolation, 'r'); 
hold on; 
pcshow(data_ANN', 'g'); 
xlabel('x');ylabel('y');zlabel('z');
set(gca,'Xdir','reverse'); 
az = 0;
el = 90;
view(az, el);
subplot(1,2,2)
imshow('template.png');
hold on;
scatter(pixel_test(:,2),pixel_test(:,1), 1);

%% Write the data to the csv file 
csvwrite('xyz_points.csv', data_Interpolation);

%% Test the interpolation method 
idx = 100; % The ith value;
figure(1);clf
subplot(1,2,1)
pcshow(P_MTI_3d,'b');
hold on;
pcshow(Cell_3d{idx}, 'r', 'MarkerSize', 20);
xlabel('x');ylabel('y');zlabel('z');
set(gca,'Xdir','reverse'); 
az = 0;
el = 90;
view(az, el);
subplot(1,2,2)
imshow('template.png');
hold on;
scatter(Cell_2d{idx}(:,1), Cell_2d{idx}(:,2), 1);
