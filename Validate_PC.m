% This is the code to check the 3D calibration
% Maybe this is the last chance
clc;clear all;close all;
%% Load the data 
P_MTI_3d = readNPY('C:\Users\gm143\Documents\MATLAB\BTL\Data\exp_3D\P_MIT_target_3d.npy');
P_MTI_2d = readNPY('C:\Users\gm143\Documents\MATLAB\BTL\Data\exp_3D\P_MTI_Calibration_interpolation_2D.npy');
P_MTI_2d = double(P_MTI_2d); 
MTI_tex = readNPY('C:\Users\gm143\Documents\MATLAB\BTL\Data\exp_3D\MTI_textures_calibrate.npy');
MTI_vtx = readNPY('C:\Users\gm143\Documents\MATLAB\BTL\Data\exp_3D\MTI_vertices_calibrate.npy');

MTI_tex_target = readNPY('C:\Users\gm143\Documents\MATLAB\BTL\Data\exp_3D\MTI_textures.npy');
MTI_vtx_target = readNPY('C:\Users\gm143\Documents\MATLAB\BTL\Data\exp_3D\MTI_vertices.npy');

imgpath = strcat('C:\Users\gm143\Documents\MATLAB\BTL\Data\exp_3D\Calibration_2D\', 'P_1.png');
imgtarget = 'C:\Users\gm143\Documents\MATLAB\BTL\Data\exp_3D\target.png';
img_calibrate = imread(imgpath);
img_target = imread(imgtarget);
%% VTX and TEX Preprocessing -- calibrate
% Remove the 0 index
idx_remove = (MTI_vtx(:,1) == 0) & (MTI_vtx(:,2) == 0) & (MTI_vtx(:,2) == 0);
MTI_vtx = MTI_vtx(~idx_remove,:);
MTI_tex = MTI_tex(~idx_remove,:);
% Texture to pixel
MTI_tex(:,1) = int16(MTI_tex(:,1) * 640);
MTI_tex(:,2) = int16(MTI_tex(:,2) * 480);
% Remove the minus index
idx_remove_1 = MTI_tex(:,1) >= 640 | MTI_tex(:,1) <= 0 | MTI_tex(:,2) >= 480 | MTI_tex(:,2) <= 0;
MTI_vtx = MTI_vtx(~idx_remove_1,:);
MTI_tex = MTI_tex(~idx_remove_1,:);
% Color vector 
colorvec = [];
for i = 1:length(MTI_tex(:,1))
    idx_width = MTI_tex(i,1);
    idx_height = MTI_tex(i,2);
    rgb = reshape(img_calibrate(idx_height, idx_width, :), [1,3]);
    colorvec(i,1) = rgb(1); colorvec(i,2) = rgb(2); colorvec(i,3) = rgb(3);
end
colorvec = uint8(colorvec);

%% VTX and TEX Preprocessing -- target
% Remove the 0 index
idx_remove = (MTI_vtx_target(:,1) == 0) & (MTI_vtx_target(:,2) == 0) & (MTI_vtx_target(:,2) == 0);
MTI_vtx_target = MTI_vtx_target(~idx_remove,:);
MTI_tex_target = MTI_tex_target(~idx_remove,:);
% Texture to pixel
MTI_tex_target(:,1) = int16(MTI_tex_target(:,1) * 640);
MTI_tex_target(:,2) = int16(MTI_tex_target(:,2) * 480);
% Remove the minus index
idx_remove_1 = MTI_tex_target(:,1) >= 640 | MTI_tex_target(:,1) <= 0 | MTI_tex_target(:,2) >= 480 | MTI_tex_target(:,2) <= 0;
MTI_vtx_target = MTI_vtx_target(~idx_remove_1,:);
MTI_tex_target = MTI_tex_target(~idx_remove_1,:);
% Color vector 
colorvec_target = [];
for i = 1:length(MTI_tex_target(:,1))
    idx_width = MTI_tex_target(i,1);
    idx_height = MTI_tex_target(i,2);
    rgb = reshape(img_target(idx_height, idx_width, :), [1,3]);
    colorvec_target(i,1) = rgb(1); colorvec_target(i,2) = rgb(2); colorvec_target(i,3) = rgb(3);
end
colorvec_target = uint8(colorvec_target);

%% Stereo 3D
CENTROID = [];
r = 3;
for i = 1:length(P_MTI_2d)
i 
p1 = P_MTI_2d(i,:) + [-r, -r];
p2 = P_MTI_2d(i,:) + [r, r];

x1 = p1(1); y1 = p1(2);
x2 = p2(1); y2 = p2(2);

% 3D colorized point cloud
idx_use = MTI_tex(:,1) >= x1 & MTI_tex(:,1) <= x2 & MTI_tex(:,2) >= y1 & MTI_tex(:,2) <= y2;
tex_use = MTI_tex(idx_use,:);
vtx_use = MTI_vtx(idx_use,:);
color_use = colorvec(idx_use, :);
CENTROID(i,:) = mean(vtx_use);
end
figure(1);clf;pcshow(CENTROID, 'r');
xlabel('x'); ylabel('y'); zlabel('z');
%% 3d to 3d 
idx_miss = [];
P_MTI_3d = P_MTI_3d * 2.54; % cm 
P_STEREO_3d = CENTROID;               
P_STEREO_3d = P_STEREO_3d .* 100;   % cm
P_STEREO_3d(idx_miss,:) = [];
P_MTI_3d(idx_miss, :) = [];
P_MTI_3d = P_MTI_3d(1:length(P_STEREO_3d), :);
figure(1);clf;pcshow(P_STEREO_3d, 'r'); hold on; pcshow(P_MTI_3d, 'b');

[R_final, t_final, P_MTI_tform, P_Stereo_tform] = STEREO_MTI_MAP(P_MTI_3d, P_STEREO_3d);
figure(1);clf; pcshow(P_MTI_tform, 'r'); hold on; pcshow(P_Stereo_tform, 'b');

%% Measure error
error_mean = (P_MTI_tform - P_Stereo_tform) 
% measure in mm
max_x = max(error_mean(:,1)) * 10
max_y = max(error_mean(:,2)) * 10
max_z = max(error_mean(:,3)) * 10

%% Target regions and point cloud
img_name = 'C:\Users\gm143\Documents\MATLAB\BTL\Data\exp_3D\target.png';
img_seg = imread(img_name);
img_seg_grey = rgb2gray(img_seg);
figure(1);clf
[img_seg_crop, rect] = imcrop(img_seg);

[pixel_test, image_dilated, ] = imgColorSeg(img_seg_crop, 2);

% save the crop image
pixel_test_use = pixel_test;
pixel_test_use(:,1) = pixel_test_use(:,1) + rect(2);
pixel_test_use(:,2) = pixel_test_use(:,2) + rect(1);

% check the results for the segmented track 
figure(1);clf;imshow(img_seg_grey); hold on;
plot(pixel_test_use(:,2), pixel_test_use(:,1), 'rx');

%%
% img_name = 'C:\Users\gm143\Documents\MATLAB\BTL\Data\exp_3D\target.png';
% img_seg = imread(img_name);
% [pixel_test, image_dilated, []] = imgColorSeg(img_seg, 10);
% pixel_test_use = pixel_test;
%%
% img_name = 'C:\Users\gm143\Documents\MATLAB\BTL\Data\exp_3D\target.png';
% [pixel_test, image_dilated] = process_tattoo_image(img_name);
% pixel_test_use = pixel_test;
%% convert to color pc
TARGET = [];
colorvec_target_1 = [];
r = 2;
for i = 1:length(pixel_test_use)
i 
p1 = pixel_test_use(i,:) + [-r, -r];
p2 = pixel_test_use(i,:) + [r, r];

x1 = p1(1); y1 = p1(2);
x2 = p2(1); y2 = p2(2);

% 3D colorized point cloud
idx_use = MTI_tex_target(:,1) >= y1 & MTI_tex_target(:,1) <= y2 & MTI_tex_target(:,2) >= x1 & MTI_tex_target(:,2) <= x2;
tex_use = MTI_tex_target(idx_use,:);
vtx_use = MTI_vtx_target(idx_use,:);
color_use = colorvec_target(idx_use, :);
TARGET(i,:) = mean(vtx_use);
colorvec_target_1(i,:) = color_use(1,:);
end

%% denoise the target point cloud
idx_noise = TARGET(:,1) > 0.04;
TARGET(idx_noise,:) = [];
colorvec_target_1(idx_noise,:) = [];
figure(1);clf;pcshow(TARGET, uint8(colorvec_target_1));
xlabel('x'); ylabel('y'); zlabel('z');

%% color target MTI to stereo  
TARGET(isnan(TARGET(:,1)),:) = [];
TARGET = TARGET * 100; % cm
TARGET_mti = TARGET * R_final + t_final;

TARGET_mti = TARGET_mti ./ 2.54;
P_MTI_tform = P_MTI_tform ./ 2.54;

figure(1);clf;
% pcshow(TARGET, 'r');
xlabel('x'); ylabel('y'); zlabel('z'); 
hold on; 
pcshow(TARGET_mti, 'b'); hold on ;
pcshow(P_MTI_tform, 'g'); hold on ;

%% save the data
csvwrite('xyz_points.csv', TARGET_mti);
