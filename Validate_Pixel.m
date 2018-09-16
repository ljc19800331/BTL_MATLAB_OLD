%% Post processing
clc;clear;close all;
%% Load the data
foldername = 'C:\Users\gm143\Documents\MATLAB\BTL\Data\CalibratePixel\exp_1\';

folderbefore = strcat(foldername, 'before\');
MTI_tex_before = readNPY(strcat(folderbefore, 'MTI_textures.npy'));
MTI_vtx_before = readNPY(strcat(folderbefore, 'MTI_vertices.npy'));
imgbefore = strcat(folderbefore, 'target.png');
img_before = imread(imgbefore);

folderafter = strcat(foldername, 'after\');
MTI_tex_after = readNPY(strcat(folderafter, 'MTI_textures.npy'));
MTI_vtx_after = readNPY(strcat(folderafter, 'MTI_vertices.npy'));
imgafter = strcat(folderafter, 'target.png');
img_after = imread(imgafter);

%% 
pixel_before;
pixel_after;

%% Method 1
[img_crop_before, rect] = imcrop(img_before);
img_crop_after = imcrop(img_after,rect);

[pixel_before, img_before_seg, img_map_before] = imgColorSeg(img_crop_before, 2);
[pixel_after, img_after_seg, img_map_after] = imgColorSeg(img_crop_after, 3);

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

pixel_label = [xv_before + rect(2), yv_before + rect(1)];
pixel_in = [xv_after(in) + rect(2), yv_after(in) + rect(1)];
pixel_out = [xv_after(~in) + rect(2), yv_after(~in) + rect(1)];
figure(1);clf; imshow(img_before); hold on; 
plot(pixel_in(:,2), pixel_in(:,1), 'rx'); hold on;
plot(pixel_out(:,2), pixel_out(:,1), 'bx'); hold on;
plot(pixel_label(:,2), pixel_label(:,1), 'gx');

%% Find the corresponding point cloud
pixel_3; 
pixel_4;
pixel_5;
figure(1);clf;pcshow(target_label, 'r'); hold on; 
pcshow(target_in, 'b'); hold on; 
pcshow(target_out, 'g'); hold on; 
xlabel('x'); ylabel('y'); zlabel('z'); 

%% Calculate the distances
E_in = []; 
E_out = [];
% in
for i = 1:length(target_in)
    p_in = target_in(i,:);
    [idx, value] = knnsearch(target_label, p_in, 'k', 1, 'distance', 'euclidean');
    E_in(i) = value;
end

E_in = E_in .* 1000; % mm
% out
for i = 1:length(target_out)
    p_out = target_out(i,:);
    [idx, value] = knnsearch(target_label, p_out, 'k', 1, 'distance', 'euclidean');
    E_out(i) = value;
end

E_out = E_out .* 1000;

%% STAT
max_in = max(E_in);
max_out = max(E_out); 
rmse_in = sqrt(sum(E_in .^ E_in)./length(E_in));
rmse_out = sqrt(sum(E_out .^ E_out)./length(E_out));
N_in = length(E_in);
N_out = length(E_out);
