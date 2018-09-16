%% calculate the pixel per point
foldername = 'C:\Users\gm143\Documents\MATLAB\BTL\Data\CalibratePixel\exp_2\before\';
MTI_tex_target = readNPY(strcat(foldername, 'MTI_textures.npy'));
MTI_vtx_target = readNPY(strcat(foldername, 'MTI_vertices.npy'));
imgtarget = strcat(foldername, 'target.png');
img_target = imread(imgtarget);

%% Target regions and point cloud
img_seg = imread(imgtarget);
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
idx_noise = TARGET(:,1) > 0.04 |  TARGET(:,3) > 0.18;
TARGET(idx_noise,:) = [];
colorvec_target_1(idx_noise,:) = [];
figure(1);clf;pcshow(TARGET, uint8(colorvec_target_1));
xlabel('x'); ylabel('y'); zlabel('z');

%% 
% R_final = [-0.9995    0.0077    0.0309
%            -0.0166   -0.9542   -0.2989
%             0.0272   -0.2992    0.9538]; 
% t_final = [-2.7846, 5.6894, -0.2564];

R_final = [-0.9997   -0.0129    0.0218
            0.0060   -0.9570   -0.2899
            0.0246   -0.2897    0.9568]; 
t_final = [-2.6552    5.5230   -0.3400];

%% color target MTI to stereo  
TARGET(isnan(TARGET(:,1)),:) = [];
TARGET = TARGET * 100; % cm
TARGET_mti = TARGET * R_final + t_final;

TARGET_mti = TARGET_mti ./ 2.54;
% P_MTI_tform = P_MTI_tform ./ 2.54;

figure(1);clf;
% pcshow(TARGET, 'r');
xlabel('x'); ylabel('y'); zlabel('z'); 
hold on; 
pcshow(TARGET_mti, 'b'); hold on ;
% pcshow(P_MTI_tform, 'g'); hold on ;

%% save the data
csvwrite('xyz_points.csv', TARGET_mti);

