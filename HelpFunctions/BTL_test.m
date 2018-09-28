% test case 
%% test the sec_pre class
figure(1);clf
elp_time = [];
test = sec_pre();

% 1.txt2pc
tic;
subplot(2,3,1);
% txt_x = '/home/maguangshen/Desktop/matlab_bc/MATLAB/buffer/Data_nonrot_5_100/Scan_x';
% txt_y = '/home/maguangshen/Desktop/matlab_bc/MATLAB/buffer/Data_nonrot_5_100/Scan_y';
% txt_z = '/home/maguangshen/Desktop/matlab_bc/MATLAB/buffer/Data_nonrot_5_100/Scan_z';
txt_x = 'Scan_x'; txt_y = 'Scan_y'; txt_z = 'Scan_z';
pc_test = test.txt2pc(txt_x, txt_y, txt_z);
toc;
elp_time(1) = toc; 
pcshow(pc_test); grid on; title('Test txt2pc section'); 

% 2.stl2pc
tic;
subplot(2,3,2);
stl_name = 'brain_slice.stl';
[pc_test, data_xyz] = test.stl2pc(stl_name);
toc; 
elp_time(2) = toc; 
pcshow(pc_test); grid on; title('Test stl2pc section'); 

% 3.pcLenChange
subplot(2,3,3); 
tic;
flag_type = 'brain';
flag_Len = 'inch';
pc_test = test.pcLenChange(pc_test, flag_type, flag_Len);
toc; elp_time(3) = toc; 
pcshow(pc_test); grid on; title('Test pcLenChange section'); 

% 4.BrainInit
subplot(2,3,4);
tic;
pc_test = test.BrainInit(pc_test);
toc; elp_time(4) = toc;
pcshow(pc_test); grid on; title('Test BrainInit section'); 

% 5.Brain2Img
subplot(2,3,5);
tic;
n_x = 669; n_y = 910;
img_brain = test.Brain2Img(pc_test, n_x, n_y);
elp_time(5) = toc;
image(img_brain, 'CDataMapping', 'scaled'); axis equal;
grid on; title('Test Brain2Img section'); 

% 6.Scan2Img 
subplot(2,3,6);
tic;
pixel_range = 100;
pc_scan = test.txt2pc(txt_x, txt_y, txt_z); 
img_scan = test.Scan2Img(pc_scan, pixel_range);
toc; elp_time(6) = toc;
image(img_scan, 'CDataMapping', 'scaled'); axis equal;
grid on; title('Test Scan2Img section'); 

pc_brain = pc_test;

% 7. Test ImgGen operator
% flag_dof = '6D';
% flag_box = 2; 
% flag_iter = 5;
% Range_mat = [0, 0, -45, 0, 0, 0];
% test.ImgGen6D(pc_brain.Location, flag_dof, flag_box, flag_iter, Range_mat);

%% Test sec_map 
% Test the Coarse Registration
clear Result;
map = sec_map();
flag_iter = 15;
tic
Range_mat = [-45, 0, -15, 30, -45, 0];
[Result] = map.CoarseRgs(img_scan, Range_mat, flag_iter);
toc;

% Test the Fine Registration
flag_box = 2;
[R_output, T_output, pc_obj_2, pc_two, Agl_mat, C_max] = map.FineRgs(pc_brain, pc_scan, Result, Range_mat, flag_box, flag_iter); 
figure(1);clf
pcshow(pc_brain);hold on; 
pc_transformed = pc_obj_2 * R_output + T_output;
pcshow(pc_transformed, 'r'); 

% Viz the convergence region
[id_x, id_y, id_z] = size(C_max);
[X,Y,Z] = ndgrid(1:size(C_max,1), 1:size(C_max,2), 1:size(C_max,3));
pointsize = 30;
scatter3(X(:), Y(:), Z(:), pointsize, C_max(:), 'MarkerFaceColor',[0 1 1]);
rotate3d on; hold on;
[idx] = find(C_max > 0.92);
scatter3(X(idx), Y(idx), Z(idx), pointsize, C_max(idx), 'MarkerFaceColor',[0.5 0.5 0.5]);
          
%% Test the sec_thesis file

% Decode the Cmax
load exp_3_3D_case5.mat;
thesis = sec_thesis();
[CMAX] = thesis.CMAX_MAT(Res_Cmax); 

pre = sec_pre();
stl_name = 'brain_slice.stl';
[pc_brain, data_xyz] = pre.stl2pc(stl_name);
flag_type = 'brain';
flag_Len = 'inch';
pc_brain = pre.pcLenChange(pc_brain, flag_type, flag_Len);
pc_brain = pre.BrainInit(pc_brain);

ERR = [];
ERR_ICP = [];
ERR_ICP2 = [];
flag_viz = 0;
for flag_p = 1:5
for flag_angle = 1:5
%     flag_p = 2; % point index
%     flag_angle = 1;
    Range_mat = [-10, 50, -50, 50, -50, 50];
    % Non ICP
    [err, err_icp, err_icp2] = thesis.Err_6D(flag_p, flag_angle, flag_viz, CMAX, Range_mat, pc_brain);
    % ICP error
    ERR(flag_p, flag_angle) = err; 
    ERR_ICP(flag_p, flag_angle) = err_icp;
    ERR_ICP2(flag_p, flag_angle) = err_icp2;
end
end

% Calculate the RMSE (root mean square error)
RMSE = sqrt(sum(sum(ERR .^ 2))./ size(ERR(:),1)); % The RMSE is 0.4 mm
RMSE_ICP = sqrt(sum(sum(ERR_ICP .^ 2))./ size(ERR_ICP(:),1)); % The RMSE is 0.4 mm
RMSE_ICP2 = sqrt(sum(sum(ERR_ICP2 .^ 2))./ size(ERR_ICP2(:),1)); % The RMSE is 0.4 mm

% test the Err_3D
[RMSE_NOHPR, RMSE_HPR, AVG_RMSE_NOHPR, AVG_RMSE_HPR, RATE_NOHPR, RATE_HPR, rate_nohpr_case, rate_hpr_case] = thesis.Err_3D();
% test sec_viz 4D
flag_iter = 5;

flag_p = 4; flag_angle = 5;
thesis.VIZ_4D(CMAX, Range_mat, flag_iter, flag_p, flag_angle);
% test sec_viz 5D
flag_p = 4; flag_angle = 5;
thesis.VIZ_5D(CMAX, Range_mat, flag_iter, flag_p, flag_angle);
% test sec_viz 6D 
flag_p = 4; flag_angle = 5;
thesis.VIZ_6D(Res_Cmax, flag_p, flag_angle);

%% Modify the STL file
stlscaled(); 

%% test sec_vessel
% problem definition: circle, pink color

%% Save to the pcl or pcd file
test = sec_pre();
stl_name = 'brain_slice.stl';
[pc_test, data_xyz] = test.stl2pc(stl_name);
pc_test = test.pcLenChange(pc_test, flag_type, flag_Len);
pc_test = test.BrainInit(pc_test);
pc_xyz = pc_test.Location;
pcwrite(pc_test,'brain_slice.pcd','Encoding','ascii');
pcwrite(pc_test,'brain_slice.ply','Encoding','ascii');

%% sec_vessel
test = sec_vessel();
Lout = test.m1();

% segment the regions
Res = [];
figure(1);

for idx_test = 1: test.NumColors
    
subplot(3,4, idx_test);
img_obj = imread(test.imgpath);
img_gray = rgb2gray(img_obj);
img_gray(Lout==idx_test) = 255;
img_gray(Lout~=idx_test) = 0;

% Suppose the number of phantom is either 1 or 2 
% This still needs to be hard code

img_use = img_gray;
img_bin = imbinarize(img_use);
L = bwlabel(img_bin);
imshow(img_use);  

end 

%% Vessel detection -- Gabor filter
img = imread('img_vessel_3.jpg');
img_gray = rgb2gray(img);
imshow(img_gray);
imwrite(img_gray, 'img_vessel_gray.jpg');

