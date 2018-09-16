% Validation based on the 3D method
%% Load the data
file_before = 'C:\Users\gm143\Documents\MATLAB\BTL\Data\exp_14\P_MTI_Calibration_interpolation_3D_before.npy';
file_after = 'C:\Users\gm143\Documents\MATLAB\BTL\Data\exp_14\P_MTI_Calibration_interpolation_3D_after.npy';
P_before = readNPY(file_before);
P_after = readNPY(file_after);

figure(1);clf
pcshow(P_before, 'r'); hold on;
pcshow(P_after, 'b');
xlabel('x'); ylabel('y'); zlabel('z');

%% Preprocessing the data set

N = double(uint8(sqrt(length(P_before))));

% P_before
xmax_before = max(P_before(:,1));
xmin_before = min(P_before(:,1)); 
ymax_before = max(P_before(:,2));
ymin_before = min(P_before(:,2));

x_line = linspace(xmin_before, xmax_before, N);
y_line = linspace(ymin_before, ymax_before, N);
[xq, yq] = meshgrid(x_line, y_line);

vq = griddata(P_before(:,1), P_before(:,2), P_before(:,3), xq, yq);
colormap_before = vq;

% P_after
xmax_after = max(P_after(:,1));
xmin_after = min(P_after(:,1)); 
ymax_after = max(P_after(:,2));
ymin_after = min(P_after(:,2));

x_line = linspace(xmin_after, xmax_after, N);
y_line = linspace(ymin_after, ymax_after, N);
[xq, yq] = meshgrid(x_line, y_line);

vq = griddata(P_after(:,1), P_after(:,2), P_after(:,3), xq, yq);
colormap_after = vq;

colormap_before = imrotate(colormap_before, 180); 
colormap_after = imrotate(colormap_after, 180); 

figure(1);clf
subplot(1,2,1)
imshow(colormap_before,[]); title('before');
xlabel('x'); ylabel('y');
subplot(1,2,2)
imshow(colormap_after,[]); title('after');
xlabel('x'); ylabel('y');

%% The label centroid for the problem
color_label = colormap_before;
thres_1 = 6.95;
color_label(color_label >= thres_1) = 0; 

thres_2 = 7.00;
color_result = colormap_after; 
color_result(color_result <= thres_2) = 0; 

figure(1);clf 
subplot(1,2,1)
imshow(color_label,[]);
subplot(1,2,2)
imshow(color_result,[]);

% case 14
label = [22,18; 52,18; 81,18; 23,47; 52,47; 81,47; 23,77; 52,77; 81,77];
result = [23,19; 52,19; 82,18; 23,48; 51,48; 83,49; 22,81; 52,80; 81,79];

% case 13
label = [22,18; 52,18; 81,18; 23,48; 52,48; 81,49; 23,78; 53,78; 82,78];
result = [22,19; 52,19; 81,19; 22,47; 51,47; 81,49; 22,78; 52,78; 81,79];

% case 11
label = [24,11; 52,12; 81,12; 24,42; 52,41; 81,42; 23,71; 52,71; 81,71];
result = [23,13; 52,13; 83,16; 23,43; 52,44; 81,44; 22,72; 52,72; 81,73];

% case 10
label = [11,9; 26,9; 40,9; 11,24; 26,24; 40,24; 12,38; 26,39; 41,39];
result = [9,10; 24,9; 38,10; 9,25; 24,25; 39,25; 9,40; 24,39; 38,39];

% case 9
label = [11,8; 26,8; 40,8; 11,22; 26,22; 40,22; 11,37; 25,37; 39,38];
result = [8,8; 23,9; 38,9; 9,24; 24,23; 38,23; 9,38; 24,41; 39,40];

DIFF = label - result;
error = sqrt(sum(    (sqrt((DIFF(:,1).^2 + DIFF(:,2).^2))).^2   )./9);
error_mm = 1./N * 2.54 * 10 * 2.67; % 0.6850 mm

%% Generate the depth map
colormap_before(isnan(colormap_before)) = 0; 
colormap_after(isnan(colormap_after)) = 0;
diff = colormap_after - colormap_before;

threshold = 0.03;
diff_org = diff;
diff(diff < threshold | diff > 1) = 0;
figure(1);clf;
subplot(2,2,1)
imshow(color_label,[]); title('label');
subplot(2,2,2)
imshow(diff,[]);
subplot(2,2,3)
imshow(colormap_before,[]);

%%
diff(diff > 0) = 1;
diff(isnan(diff)) = 0;
figure(2);clf
subplot(2,2,1)
imshow(colormap_before,[]); 
xlabel('x'); ylabel('y');
title('Scanning data before ablation', 'Fontsize', 20);
subplot(2,2,2)
imshow(colormap_after,[]);
xlabel('x'); ylabel('y');
title('Scanning data after ablation', 'Fontsize', 20);
subplot(2,2,3)
imshow(diff_org,[]); 
xlabel('x'); ylabel('y');
title('Error map before processing', 'Fontsize', 20);
subplot(2,2,4)
imshow(diff,[]);
title('Error map after processing', 'Fontsize', 20);


%% Measure the error for a single point




