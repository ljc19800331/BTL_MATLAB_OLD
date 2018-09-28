% Read the data
P_before = []; 
P_after = [];
x_before = load('C:\Users\gm143\Documents\MATLAB\Brain_MATLAB\5d_image\StereoMTI\test_3\Scan_x_before');
y_before = load('C:\Users\gm143\Documents\MATLAB\Brain_MATLAB\5d_image\StereoMTI\test_3\Scan_y_before');
z_before = load('C:\Users\gm143\Documents\MATLAB\Brain_MATLAB\5d_image\StereoMTI\test_3\Scan_z_before');

x_after = load('C:\Users\gm143\Documents\MATLAB\Brain_MATLAB\5d_image\StereoMTI\test_3\Scan_x_after');
y_after = load('C:\Users\gm143\Documents\MATLAB\Brain_MATLAB\5d_image\StereoMTI\test_3\Scan_y_after');
z_after = load('C:\Users\gm143\Documents\MATLAB\Brain_MATLAB\5d_image\StereoMTI\test_3\Scan_z_after');

P_before(:,1) = x_before; P_before(:,2) = y_before; P_before(:,3) = z_before; 
P_after(:,1) = x_after; P_after(:,2) = y_after; P_after(:,3) = z_after;

figure(1);clf
pcshow(P_before, 'r'); hold on;
pcshow(P_after, 'b');

%% Normalized the data 
P_before(:,1) = P_before(:,1) - min(P_before(:,1)); 
P_before(:,2) = P_before(:,2) - min(P_before(:,2)); 
P_before(:,3) = P_before(:,3) - min(P_before(:,3)); 

P_after(:,1) = P_after(:,1) - min(P_after(:,1));
P_after(:,2) = P_after(:,2) - min(P_after(:,2));
P_after(:,3) = P_after(:,3) - min(P_after(:,3));

figure(1);clf
pcshow(P_before, 'r'); hold on;
pcshow(P_after, 'b');
xlabel('x'); ylabel('y'); zlabel('z');

%% Create the colop map -- project along the z axis
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
imshow(colormap_before,[]); 
xlabel('x'); ylabel('y');
subplot(1,2,2)
imshow(colormap_after,[]);
xlabel('x'); ylabel('y');

%% Calculate the difference 
diff = colormap_before - colormap_after;
diff(diff < 0.01) = 0;
diff_org = diff;
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

%% Manul measure the pixel error
p_test = [9,8; 26,8; 41,8; 55,8; 9,23; 27,23; 41,23; 58,23; 9,39; 25,41; 40,40; 58,40; 9,55; 25,56; 41,55; 58,55];
label_x = repmat(linspace(8, 56, 4), 1,4);
label_y = repmat(linspace(8, 56, 4), 4,1);
label_y = label_y(:);
% label_mesh = meshgrid(label_x, label_y);
p_label = [label_x',label_y];

diff_testlabel = (p_test - p_label);
error = sum(sqrt(diff_testlabel(:,1).^2 + diff_testlabel(:,2).^2))./length(diff_testlabel);
InchPerPixel = 0.02;    % inch
mTRE_inch = InchPerPixel * error; % inch 
mTRE_mm = mTRE_inch * 2.54;

%% Find the connected region
[L,n] = bwlabel(diff, 8);
centroids = [];
for i = 1 : n
    idx_L = find(L == i);
    [I,J] = ind2sub(size(diff), idx_L); 
    centroids(i,:) = double([uint8(mean(I)), uint8(mean(J))]);
end
diff_result = zeros(length(diff), length(diff));
for i = 1 : length(centroids)
diff_result(centroids(i,1), centroids(i,2)) = 1; 
end
figure(2);clf
imshow(diff_result, []);

label_x = linspace(8, 56, 4);
label_y = linspace(8, 56, 4);
label_mesh = meshgrid(label_x, label_y); 

label_map = zeros(length(diff), length(diff));
count = 1;
for i = 1 : length(label_x)
    for j = 1: length(label_y)
        label_map(centroids(count,1), centroids(count,2)) = 1; 
        label_map(label_x(i), label_y(j)) = 1; 
        count = count + 1; 
    end
end
figure(2);clf
subplot(1,2,1)
imshow(diff_result, []); 
subplot(1,2,2)
imshow(label_map, []); 

% Calculate the pixel values -mTRE
idx_exp = find(diff_result == 1);
[coord_x_exp, coord_y_exp] = ind2sub(size(diff), idx_exp);

idx_label = find(label_map);
[coord_x_label, coord_y_label] = ind2sub(size(diff), idx_label);

error = sum(sqrt((coord_x_exp - coord_x_label).^2 + (coord_y_exp - coord_y_label).^2))./ n;

InchPerPixel = 0.02;
mTRE = InchPerPixel * error;
