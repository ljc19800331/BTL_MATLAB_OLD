%% Method 1: ANN
test = BTL_use();
flag_method = 2;
[pixel_test] = test.TestDesign(flag_method);
X_laser = pixel_test(:,2); 
Y_laser = pixel_test(:,1);
[net, P_MTI_2d, P_MTI_3d] = test.Calibration_ANN();

% measure the max error for ANN
% [net, P_MTI_2d, P_MTI_3d] = test.Calibration_ANN();
data_ANN_test = net([P_MTI_2d(:,1), P_MTI_2d(:,2)]');

dist = data_ANN_test' - P_MTI_3d;
err = sqrt(dist(:,1).^2 + dist(:,2).^2);
max(err) * 25.4 % 2mm

figure(1);clf
pcshow(P_MTI_3d, 'b'); hold on; 
pcshow(data_ANN_test', 'r'); hold on;

data_ANN = net([pixel_test(:,2), pixel_test(:,1)]'); 

%% Method 2: Interpolation
test = BTL_use();
% flag_method = 2;
% [pixel_test] = test.TestDesign(flag_method);
X_laser = pixel_test(:,2);          % height -- 1 
Y_laser = pixel_test(:,1);          % width -- 2
[data_Interpolation, Cell_2d, Cell_3d, P_MTI_2d, P_MTI_3d] = test.Calibration_Interpolation(X_laser, Y_laser); 

%% Show the results in the same object
figure(1);clf
subplot(1,2,1)
hold on;
pcshow(P_MTI_3d,'b');
hold on;
pcshow(data_Interpolation, 'r'); 
hold on; 
pcshow(data_ANN', 'g'); 
hold on;

xlabel('x');ylabel('y');zlabel('z');
set(gca,'Xdir','reverse'); 
az = 0;
el = 90;
view(az, el);
subplot(1,2,2)
imshow('before.png');
hold on;
scatter(pixel_test(:,2),pixel_test(:,1), 1);

%% Reorder the pixels in a close region
x_max_target = max(data_Interpolation(:,1)); x_min_target = min(data_Interpolation(:,1));
y_max_target = max(data_Interpolation(:,2)); y_min_target = min(data_Interpolation(:,2));
z_max_target = max(data_Interpolation(:,3)); z_min_target = min(data_Interpolation(:,3));

x_range = x_max_target - x_min_target;
y_range = y_max_target - y_min_target;
z_range = z_max_target - z_min_target; 

% the unit is based on the spot size of the laser center
x_line = linspace(x_min_target, x_max_target, 100);
y_line = linspace(y_min_target, y_max_target, 100);
[x_grid, y_grid] = meshgrid(x_line, y_line);
vq = griddata(data_Interpolation(:,1), data_Interpolation(:,2), data_Interpolation(:,3), x_grid, y_grid);

figure(1);clf
subplot(1,2,1)
pcshow(P_MTI_3d,'b');
hold on;
% pcshow(data_Interpolation, 'r'); 
% hold on;
pcshow([x_grid(:), y_grid(:), vq(:)], 'r'); 
hold on; 
pcshow(data_ANN', 'g'); 
xlabel('x');ylabel('y');zlabel('z');
set(gca,'Xdir','reverse'); 
az = 0;
el = 90;
view(az, el);
subplot(1,2,2)
imshow('before.png');
hold on;
scatter(pixel_test(:,2),pixel_test(:,1), 1);

data_reorder = [x_grid(:), y_grid(:), vq(:)];
idx_nan = find(isnan(data_reorder(:,3))==1);
data_reorder(idx_nan, 3) = min(vq(:));

%% Write the data to the csv file 
csvwrite('xyz_points.csv', data_ANN');
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
imshow('before.png');
hold on;
scatter(Cell_2d{idx}(:,1), Cell_2d{idx}(:,2), 1);

%% 
% figure(1); clf; 
% pcshow(P_MTI_3d, 'r');
% xlabel('x');ylabel('y');zlabel('z');
% set(gca,'Xdir','reverse'); 
% az = 0;
% el = 90;
% view(az, el);
% figure(2);clf
% imshow('before.png'); hold on;
% plot(P_MTI_2d(:,1), P_MTI_2d(:,2), 'b');


