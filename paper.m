%% 
A = P_STEREO_3d; 
B = P_MTI_3d;

% A = A2; 
% B = B;

R = orth(rand(3,3)); % random rotation matrix

if det(R) < 0
    V(:,3) = V(:,3) * -1;
    R = V*U';
end

t = rand(3,1); % random translation

n = 10; % number of points
% A = rand(n,3);
% B = R*A' + repmat(t, 1, n);
% B = B';

[ret_R, ret_t] = rigid_transform_3D(A, B);
n = length(P_STEREO_3d);
A2 = (ret_R*A') + repmat(ret_t, 1, n)
A2 = A2'

% Find the error
err = A2 - B;
err = err .* err;
err = sum(err(:));
rmse = sqrt(err/n); % cm

disp(sprintf("RMSE: %f", rmse));
disp("If RMSE is near zero, the function is correct!");

figure(1);clf 
pcshow(A2, 'r'); hold on;
pcshow(B, 'b'); 

%%
[R_final, t_final, P_MTI_tform, P_Stereo_tform] = STEREO_MTI_MAP(B, A2);
figure(1);clf; pcshow(P_MTI_tform, 'r'); hold on; pcshow(P_Stereo_tform, 'b');

n = length(P_MTI_tform);
err = P_MTI_tform - P_Stereo_tform;
err = err .* err;
err = sum(err(:));
rmse = sqrt(err/n);

%% 


% img_before = imread('before_1280.png');
img_after = imread('after.png'); 
% figure(1);clf; imshow(img_before); 

img_seg = img_after;
img_seg_grey = rgb2gray(img_seg);
figure(1);clf
[img_seg_crop, rect] = imcrop(img_seg);

[pixel_test, image_dilated] = imgColorSeg(img_seg_crop, 3);

% save the crop image
pixel_test_use = pixel_test;
pixel_test_use(:,1) = pixel_test_use(:,1) + rect(2);
pixel_test_use(:,2) = pixel_test_use(:,2) + rect(1);

% check the results for the segmented track 
figure(2);clf;imshow(img_seg_grey); hold on;
plot(pixel_test_use(:,2), pixel_test_use(:,1), 'rx');

%% 
figure(1);clf 
subplot(1,2,1) 
pcshow(MTI_vtx_target, colorvec_target);
subplot(1,2,2)
pcshow(TARGET, 'r');
xlabel('x'); ylabel('y'); zlabel('z');
%% 
% raw image
img = imread('test.png');
img_roi = imcrop(img);
% grey scale image
img_grey = rgb2gray(img_roi);
% threshold the centroid grayscale color region
idx_centroid = find(img_grey>230);
[I,J] = ind2sub(size(img_grey), idx_centroid);
i_cen = mean(I); j_cen = mean(J);

% hsv img 
img_hsv = rgb2hsv(img_roi); 
% threshold the hsv color space


figure(1);clf
subplot(2,2,1)
imshow(img_roi);
subplot(2,2,2)
imshow(img_grey);
subplot(2,2,3)
imshow(img_grey); hold on;
plot(j_cen, i_cen, 'rx', 'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','b'); hold on;
pos = [(j_cen-10) (i_cen-10) 20 20]; 
rectangle('Position',pos,'Curvature',[1 1], 'EdgeColor','r', 'LineWidth',3);

%% Texture point cloud
figure(1);clf 
pcshow(MTI_vtx_target, colorvec_target); 

figure(1);clf
idx_use = MTI_tex_target(:,1) > 235 &  MTI_tex_target(:,1) < 400 & MTI_tex_target(:,2) > 159 & MTI_tex_target(:,2) < 235;

MTI_vtx_target_use = MTI_vtx_target(idx_use,:); 
colorvec_target_use = colorvec_target(idx_use,:); 

figure(1);clf 
pcshow(MTI_vtx_target_use, colorvec_target_use); 


