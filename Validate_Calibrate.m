% ref: https://igl.ethz.ch/projects/ARAP/svd_rot.pdf
% ref: http://nghiaho.com/uploads/code/rigid_transform_3D.m
% Calibration test
clc;clear all;close all;

%% Load the data 
CELL = {};
n_set = 10;

for N = 1:n_set
    
%     if N == 2
%         continue;
%     end
  
data = [];
foldername = strcat('C:\Users\gm143\Documents\MATLAB\BTL\Data\CalibrateNPNP','\pos_', num2str(N), '\');
data.P_MTI_3d = readNPY(strcat(foldername, 'P_MIT_target_3d.npy'));
data.P_MTI_2d = double(readNPY(strcat(foldername, 'P_MTI_Calibration_interpolation_2D.npy')));
data.MTI_tex = readNPY(strcat(foldername, 'MTI_textures.npy'));
data.MTI_vtx = readNPY(strcat(foldername, 'MTI_vertices.npy'));
data.imgpath = strcat(foldername, 'Calibration_2D\', 'P_1.png');
data.img_calibrate = imread(data.imgpath);
CELL{N} = data; 
end 

%% Stereo images
for N = 1 : n_set 
N
MTI_vtx = CELL{N}.MTI_vtx;
MTI_tex = CELL{N}.MTI_tex;
img_calibrate = CELL{N}.img_calibrate;
P_MTI_2d = CELL{N}.P_MTI_2d;

%% vtx and tex
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

%% Stereo 3D
CENTROID = [];
r = 2;
for i = 1:length(P_MTI_2d)

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

CELL{N}.P_STEREO_3d = CENTROID;

end

%% MTI points 
% i = 2
% n_set = 3;
MTI = [];
for i = 1:n_set
%     if i == 2 % | i == 5 | i == 6
%         continue;
%     end
    mti_3d = CELL{i}.P_MTI_3d;
    MTI = [MTI;mti_3d];
end
%% Sterero points 
STEREO = [];
for i = 1:n_set
%     if i == 2 % | i == 5 | i == 6
%         continue;
%     end
    stereo_3d = CELL{i}.P_STEREO_3d;
    STEREO = [STEREO;stereo_3d];
end
%% MAP
P_MTI_3d = MTI;
P_STEREO_3d = STEREO; 

%% 
P_MTI_3d = P_MTI_3d * 2.54; % cm 
% P_STEREO_3d = CENTROID;               
P_STEREO_3d = P_STEREO_3d .* 100;   % cm
% P_STEREO_3d(idx_miss,:) = [];
% P_MTI_3d(idx_miss, :) = [];
% P_MTI_3d = P_MTI_3d(1:length(P_STEREO_3d), :);
figure(1);clf;pcshow(P_STEREO_3d, 'r'); hold on; pcshow(P_MTI_3d, 'b');

[R_final, t_final, P_MTI_tform, P_Stereo_tform] = STEREO_MTI_MAP(P_MTI_3d, P_STEREO_3d);

% [R_final, t_final, P_MTI_tform, P_Stereo_tform] = MapLS(P_MTI_3d, P_STEREO_3d);

figure(1);clf; pcshow(P_MTI_tform, 'r'); hold on; pcshow(P_Stereo_tform, 'b');

%% RANSAC for the calibration





%% Test with the given initial estimation 


%% Test two point sets
% A = P_MTI_3d; 
% B = P_STEREO_3d; 
% [R, t] = rigid_transform_3D(A, B); 
