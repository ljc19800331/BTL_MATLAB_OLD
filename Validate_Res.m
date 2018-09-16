%% Post processing
clc;clear;close all;
CELL = {};
data = [];
%%
for N = 1:10
%% Load the image data
% N = 19; 
foldername = strcat('C:\Users\gm143\Documents\MATLAB\BTL\Data\res2\', num2str(N), '\');
name_before = strcat(foldername, 'before_1.png'); 
name_after = strcat(foldername, 'after_1.png'); 
img_before = imread(name_before); 
img_after = imread(name_after);
pep = 0.28;  % 640 
% pep = 48.6 ./ 260;  % 1280

% Method 1
% img_after = imread('after.png'); 
% figure(1);clf; imshow(img_before); 

% img_seg = img_after;
% img_seg_grey = rgb2gray(img_seg);
% figure(1);clf

[img_crop_before, rect] = imcrop(img_before);
img_crop_after = imcrop(img_after,rect);
% [img_crop_after, rect] = imcrop(img_after);

[pixel_before, img_before_seg, img_map_before] = imgColorSeg(img_crop_before, 2);
[pixel_after, img_after_seg, img_map_after] = imgColorSeg(img_crop_after, 3);

imwrite(img_crop_before, 'img_crop_before.jpg');
imwrite(img_crop_after, 'img_crop_after.jpg');
% imwrite(img_map_before, 'img_map_before.jpg'); 
% imwrite(img_map_after, 'img_map_after.jpg'); 

% save the crop image
% pixel_before_use = pixel_test;
% pixel_test_use(:,1) = pixel_test_use(:,1) + rect(2);
% pixel_test_use(:,2) = pixel_test_use(:,2) + rect(1);
% 
% % check the results for the segmented track 
% figure(2);clf;imshow(img_seg_grey); hold on;
% plot(pixel_test_use(:,2), pixel_test_use(:,1), 'rx');


%% Method 2
% [xy_coordinates, img_before_seg] = imgColorSeg(img_before); 
% [xy_coordinates, img_after_seg] = imgColorSeg(img_after); 

% [xy_coordinates, img_before_seg] = process_tattoo_image(name_before); 
% [xy_coordinates, img_after_seg] = process_tattoo_image(name_after); 

%% Filter out the boundary
mask_before = boundarymask(img_before_seg);
mask_after = boundarymask(img_after_seg);

idx_before = find(mask_before == 1);
idx_after = find(mask_after == 1);

[xv_before, yv_before] = ind2sub(size(mask_before), idx_before);
[xv_after, yv_after] = ind2sub(size(mask_after), idx_after);

[in, on] = inpolygon(xv_after, yv_after, xv_before, yv_before);

figure(1);clf
plot(yv_before, xv_before, 'gx'); % polygon
axis equal
hold on
plot(yv_after(in), xv_after(in),'r+') % points inside
plot(yv_after(~in), xv_after(~in),'bo') % points outside
xlabel('x'); ylabel('y');
grid on;
hold off
xlabel('x'); ylabel('y');

saveas(gcf, 'img_result.jpg'); 

% figure(1);clf
% plot(xv_before, yv_before, 'r+'); hold on; 
% plot(xv_after, yv_after, 'b+');
% 
% figure(1);clf
% subplot(1,2,1)
% imshow(mask_before); 
% subplot(1,2,2)
% imshow(mask_after);

% Calculate the point based distances
% positive error 
idx_positive = in;      % in this region
N_positive = length(find(idx_positive == 1)); 
% single point error 
P_before = [xv_before, yv_before];
P_positive = [xv_after(in), yv_after(in)];
E_positive = [];

for i = 1:length(P_positive)
    p_target = P_positive(i,:);
    [idx_target, value_target] = knnsearch(P_before, p_target, 'k', 1, 'distance', 'euclidean');
    E_positive(i) = value_target;
end

% mean positive error
mSRE_positive_pixel = sqrt(sum(E_positive.^2) ./ length(E_positive));  % pixel
mSRE_positive_mm = pep * mSRE_positive_pixel;                          % mm

% max positive error
max_positive_pixel = max(E_positive);
max_positive_mm = pep * max_positive_pixel;

% negative error
idx_negative = ~in;     % not in this region
N_negative = length(find(idx_negative == 1)); 
P_negative = [xv_after(~in), yv_after(~in)];
E_negative = [];

for i = 1:length(P_negative)
    p_target = P_negative(i,:);
    [idx_target, value_target] = knnsearch(P_before, p_target, 'k', 1, 'distance', 'euclidean');
    E_negative(i) = value_target;
end

% mean negative error
mSRE_negative_pixel = sqrt(sum(E_negative.^2) ./ length(E_negative));  % pixel
mSRE_negative_mm = pep * mSRE_negative_pixel;         % mm

% max negative error
max_negative_pixel = max(E_negative);
max_negative_mm = pep * max_negative_pixel;

data.mSRE_positive_mm = mSRE_positive_mm;
data.mSRE_negative_mm = mSRE_negative_mm;
data.max_positive_mm = max_positive_mm; 
data.max_negative_mm = max_negative_mm;
data.N_positive = N_positive;
data.N_negative = N_negative;
data.N_before = length(idx_before); 
data.N_after = length(idx_after);
CELL{N} = data;
%% 
end 
