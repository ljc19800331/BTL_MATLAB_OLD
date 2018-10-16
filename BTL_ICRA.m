% The class for ICRA 2019
classdef BTL_ICRA
    properties
        folder = 'C:\Users\gm143\Documents\MATLAB\BTL\Data\exp_3D\exp4\';
    end
    methods
        function [pRMSE, nRMSE, pMAE, nMAE] = Err_2D(obj)  
            
            % load the images
            folder_before = strcat(obj.folder,'before.png');
            folder_after = strcat(obj.folder, 'after.png');
            img_before = imread(folder_before);
            img_after = imread(folder_after);
            
            % ROI
            [img_crop_before, rect] = imcrop(img_before);
            img_crop_after = imcrop(img_after,rect);
            
            % Color segmentation
            [pixel_before, img_before_seg, img_map_before] = imgColorSeg(img_crop_before, 2);
            [pixel_after, img_after_seg, img_map_after] = imgColorSeg(img_crop_after, 3);
            
            % Filter out the boundary
            mask_before = boundarymask(img_before_seg);
            mask_after = boundarymask(img_after_seg);
            idx_before = find(mask_before == 1);
            idx_after = find(mask_after == 1);
            [xv_before, yv_before] = ind2sub(size(mask_before), idx_before);
            [xv_after, yv_after] = ind2sub(size(mask_after), idx_after);
            [in, on] = inpolygon(xv_after, yv_after, xv_before, yv_before);
           
            % positive error
            idx_positive = in; 
            P_before = [xv_before, yv_before];
            P_positive = [xv_after(in), yv_after(in)];
            E_positive = [];
            for i = 1:length(P_positive)
                p_target = P_positive(i,:);
                [idx_target, value_target] = knnsearch(P_before, p_target, 'k', 1, 'distance', 'euclidean');
                E_positive(i) = value_target;
            end
            
            % negative error
            idx_negative = ~in;     
            P_negative = [xv_after(~in), yv_after(~in)];
            E_negative = [];
            for i = 1:length(P_negative)
                p_target = P_negative(i,:);
                [idx_target, value_target] = knnsearch(P_before, p_target, 'k', 1, 'distance', 'euclidean');
                E_negative(i) = value_target;
            end

            % MRSE and MAE
            pRMSE = sqrt(sum(E_positive.^2) ./ length(E_positive));
            pMAE =  max(E_positive);
            nRMSE = sqrt(sum(E_negative.^2) ./ length(E_negative)); 
            nMAE = max(E_negative);
            
        end
        function [data] = Err_pixel(obj)
            
            % Load the data 
            folderbefore = strcat(obj.folder, 'before\');
            folderafter = strcat(obj.folder, 'after\');
            img_before = imread(strcat(folderbefore, 'target.png')); 
            img_after = imread(strcat(folderafter, 'target.png'));
            
            % Read the tex and vtx files
            MTI_tex_before = readNPY(strcat(folderbefore, 'MTI_textures.npy'));
            MTI_vtx_before = readNPY(strcat(folderbefore, 'MTI_vertices.npy'));
            MTI_tex_after = readNPY(strcat(folderafter, 'MTI_textures.npy'));
            MTI_vtx_after = readNPY(strcat(folderafter, 'MTI_vertices.npy'));
            
            % clear the pre-data
            [MTI_vtx_before, MTI_tex_before, colorvec_before] = obj.VtxTex_Clean(MTI_vtx_before, MTI_tex_before, img_before); 
            [MTI_vtx_after, MTI_tex_after, colorvec_after] = obj.VtxTex_Clean(MTI_vtx_after, MTI_tex_after, img_after);
        
            % crop and seg the images
            [img_crop_before, rect] = imcrop(img_before);
            img_crop_after = imcrop(img_after,rect);
            [pixel_before, img_before_seg, img_map_before] = imgColorSeg(img_crop_before, 2);
            [pixel_after, img_after_seg, img_map_after] = imgColorSeg(img_crop_after, 3);
            
            % Filter out the boundary
            mask_before = boundarymask(img_before_seg);
            mask_after = boundarymask(img_after_seg);
            idx_before = find(mask_before == 1);
            idx_after = find(mask_after == 1);
            [xv_before, yv_before] = ind2sub(size(mask_before), idx_before);
            [xv_after, yv_after] = ind2sub(size(mask_after), idx_after);
            [in, on] = inpolygon(xv_after, yv_after, xv_before, yv_before);
            
            % Recalculate the index in the color image before cropping
            pixel_label = [xv_before + rect(2), yv_before + rect(1)];
            pixel_in = [xv_after(in) + rect(2), yv_after(in) + rect(1)];
            pixel_out = [xv_after(~in) + rect(2), yv_after(~in) + rect(1)];
            
            % Calculate the colorized voxel
            [target_label, colorvec_label] = obj.Pixel2Voxel(pixel_label, colorvec_before, ...
                                                        MTI_tex_before, MTI_vtx_before);
            [target_in, colorvec_in] = obj.Pixel2Voxel(pixel_in, colorvec_before, ...
                                                        MTI_tex_before, MTI_vtx_before);
            [target_out, colorvec_out] = obj.Pixel2Voxel(pixel_out, colorvec_before, ...
                                                        MTI_tex_before, MTI_vtx_before);
            
            % Calculate the distances
            E_in = []; 
            E_out = [];
            
            % The positive pixel values 
            for i = 1:length(target_in)
                p_in = target_in(i,:);
                [idx, value] = knnsearch(target_label, p_in, 'k', 1, 'distance', 'euclidean');
                E_in(i) = value;
            end                                        
            E_in = E_in .* 1000; % mm       
            
            % The negative pixel values 
            size_out = size(target_out);
            for i = 1:length(target_out)
                if ((size_out(1) == 1) & i == 2) | ((size_out(1) == 2) & i == 3)
                    break;
                end
                p_out = target_out(i,:);
                [idx, value] = knnsearch(target_label, p_out, 'k', 1, 'distance', 'euclidean');
                E_out(i) = value;
            end
            E_out = E_out .* 1000;
            
            % Save all the data in a struct
            data = []; 
            idx_on = find(E_in == 0);
            E_in(idx_on) = []; 
            data.max_in = max(E_in);
            data.max_out = max(E_out); 
            data.rmse_in = sqrt(sum(E_in .^ 2)./length(E_in));
            data.rmse_out = sqrt(sum(E_out .^ 2)./length(E_out));
            data.N_in = length(E_in);
            data.N_out = length(E_out);
                                                    
        end 
        function [CELL] = NPNP(obj, n_set)
            % ref: https://igl.ethz.ch/projects/ARAP/svd_rot.pdf
            % ref: http://nghiaho.com/uploads/code/rigid_transform_3D.m 
            
            CELL = {};
            % n_set = 5;
            folder = 'C:\Users\gm143\Documents\MATLAB\BTL\Data\CalibrateNPNP_2\pos_'; 
            
            % Load the data
            for N = 1:n_set
                data = [];
                foldername = strcat(folder, num2str(N), '\');
                data.P_MTI_3d = readNPY(strcat(foldername, 'P_MIT_target_3d.npy'));
                data.P_MTI_3d = data.P_MTI_3d * 2.54;
                data.P_MTI_2d = double(readNPY(strcat(foldername, 'P_MTI_Calibration_interpolation_2D.npy')));
                data.MTI_tex = readNPY(strcat(foldername, 'MTI_textures.npy'));
                data.MTI_vtx = readNPY(strcat(foldername, 'MTI_vertices.npy'));
                data.imgpath = strcat(foldername, 'Calibration_2D\', 'P_1.png');
                data.img_calibrate = imread(data.imgpath);
                CELL{N} = data; 
            end 
            
            % Read the laser spot from stereo images
            for N = 1 : n_set 
                N
                MTI_vtx = CELL{N}.MTI_vtx;
                MTI_tex = CELL{N}.MTI_tex;
                img_calibrate = CELL{N}.img_calibrate;
                P_MTI_2d = CELL{N}.P_MTI_2d;
                [h, w] = size(rgb2gray(img_calibrate));
                
                % vtx and tex
                % Remove the 0 index
                idx_remove = (MTI_vtx(:,1) == 0) & (MTI_vtx(:,2) == 0) & (MTI_vtx(:,2) == 0);
                MTI_vtx = MTI_vtx(~idx_remove,:);
                MTI_tex = MTI_tex(~idx_remove,:);
                % Texture to pixel
                MTI_tex(:,1) = int16(MTI_tex(:,1) * w);
                MTI_tex(:,2) = int16(MTI_tex(:,2) * h);
                % Remove the minus index
                idx_remove_1 = MTI_tex(:,1) >= w | MTI_tex(:,1) <= 0 | MTI_tex(:,2) >= h | MTI_tex(:,2) <= 0;
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
            
                % laser spots for stereo 3D
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
                CENTROID = CENTROID * 100;
                CELL{N}.P_STEREO_3d = CENTROID;
                
            end 
        end
        function [R_best, t_best] = RANSAC(obj, CELL, n_set, n_load)
            
            % Load the 2D and 3D data
            
            P_MTI_3d = [];
            for i = 1:n_set
                mti_3d = CELL{i}.P_MTI_3d;
                P_MTI_3d = [P_MTI_3d;mti_3d];
            end
            % Sterero points 
            P_STEREO_3d = [];
            for i = 1:n_set
                stereo_3d = CELL{i}.P_STEREO_3d;
                P_STEREO_3d = [P_STEREO_3d;stereo_3d];
            end
            % Load the data
            A = P_STEREO_3d(1:n_load,:);    % moving
            B = P_MTI_3d(1:n_load,:);       % fixed

            % RANSAC parameters
            number = length(A);                 % Total number of points
            num = 50;                           % minimum number of points
            bestInNum = 0;                      % Best fitting line with largest number of inliers
            bestParameter1=0;                   % parameters for best fitting line
            bestParameter2=0;                   % Pending
            inlierRatio = 0.9;                  % inlierratio is 50% 
            threshDist = 0.188;                 % threshold distance is 0.05 -- 1.57 mm

            for iter = 1: n_load
                
            iter
            idx = randperm(number, num); 
            % Apply the transformation
            moving = A(idx,:); fixed = B(idx,:); 
            [R,t] = rigid_transform_3D(moving, fixed); 

            % Apply the transform to all other data
            n = length(A); 
            A2 = (R * A') + repmat(t, 1, n); 
            A2 = A2'; 
            err = A2 - B;
            err = err .* err;
            err = sum(err(:));
            rmse = sqrt(err/n)

            % inliers and outliers
            dis = A2 - B; 
            dis = sqrt(dis(:,1).^2 + dis(:,2).^2 + dis(:,3).^2);
            idx_in = dis < threshDist; 
            idx_out = dis >= threshDist;
            A2_in = A2(idx_in,:); 
            B_in = B(idx_in,:); 
            A2_out = A2(idx_out,:);
            B_out = B(idx_out,:);

            % continue threshold
            dis_new = A2_in - B_in;  
            dis_new = sqrt(dis_new(:,1).^2 + dis_new(:,2).^2 + dis_new(:,3).^2);
            err_new = sqrt(sum(dis_new(:).^2)./length(dis_new)); 
            N_inlier = length(dis_new);

            % condition
            if err_new < threshDist & N_inlier > inlierRatio*number
                count_best = iter;
                R_best = R
                t_best = t
                dis_best = A2_in - B_in;  
                dis_best = sqrt(dis_best(:,1).^2 + dis_best(:,2).^2 + dis_best(:,3).^2);
                err_best = sqrt(sum(dis_best(:).^2)./length(dis_best))
                break;
            end
            end
            
            % final err
            n_final = length(A);
            A_final = (R_best * A') + repmat(t_best, 1, n_final);
            A_final = A_final';
            figure(1);clf; 
            pcshow(A_final, 'r'); hold on; pcshow(B, 'b');

            err_final = A_final - B;
            err_final = err_final .* err_final;
            err_final = sum(err_final(:));
            rmse_final = sqrt(err_final/n)
            
        end
        function [TARGET_mti] = GetTargetPc(obj, R_final, t_final)
           
            % Get the target point cloud
            foldername = 'C:\Users\gm143\Documents\MATLAB\BTL\Data\exp_Phantom\pos_8\before\';
            MTI_tex_target = readNPY(strcat(foldername, 'MTI_textures.npy'));
            MTI_vtx_target = readNPY(strcat(foldername, 'MTI_vertices.npy'));
            imgtarget = strcat(foldername, 'target.png');
            img_target = imread(imgtarget);
            [h, w] = size(rgb2gray(img_target));
            
            % VTX and TEX Preprocessing -- target
            % Remove the 0 index
            idx_remove = (MTI_vtx_target(:,1) == 0) & (MTI_vtx_target(:,2) == 0) & (MTI_vtx_target(:,2) == 0);
            MTI_vtx_target = MTI_vtx_target(~idx_remove,:);
            MTI_tex_target = MTI_tex_target(~idx_remove,:);
            % Texture to pixel
            MTI_tex_target(:,1) = int16(MTI_tex_target(:,1) * w);
            MTI_tex_target(:,2) = int16(MTI_tex_target(:,2) * h);
            % Remove the minus index
            idx_remove_1 = MTI_tex_target(:,1) >= w | MTI_tex_target(:,1) <= 0 | MTI_tex_target(:,2) >= h | MTI_tex_target(:,2) <= 0;
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
            
            % Target regions and point cloud
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
            
            % Convert to colorized point cloud
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

            if (length(color_use) == 0)
                colorvec_target_1(i,:) = [0,0,0];
            else
                colorvec_target_1(i,:) = color_use(1,:);
            end
            end
            
            % denoise the target point cloud
            idx_noise = TARGET(:,1) > 0.04 |  TARGET(:,3) > 0.18;
            TARGET(idx_noise,:) = [];
            colorvec_target_1(idx_noise,:) = [];
            figure(1);clf;pcshow(TARGET, uint8(colorvec_target_1), 'MarkerSize', 50);
            figure(1);clf;pcshow(TARGET, 'r', 'MarkerSize', 50);
            xlabel('x'); ylabel('y'); zlabel('z');
            
            % Transformation information
%             R_final = [
%                         -0.9996    0.0059    0.0272
%                         -0.0137   -0.9547   -0.2973
%                          0.0243   -0.2976    0.9544 ]; 
%             t_final = [-2.7033  5.6484  -0.3263];
            
            % color target MTI to stereo  
            TARGET(isnan(TARGET(:,1)),:) = [];
            TARGET = TARGET * 100;              % cm
            n_final = length(TARGET);
            TARGET_mti = (R_final * TARGET') + repmat(t_final', 1, n_final);
            TARGET_mti = TARGET_mti';
            TARGET_mti = TARGET_mti ./ 2.54;

            figure(1);clf;
            % pcshow(TARGET, 'r');
            xlabel('x'); ylabel('y'); zlabel('z'); 
            hold on; 
            pcshow(TARGET_mti, 'b'); hold on ;
            
            csvwrite('xyz_points.csv', TARGET_mti);
            
        end
        function [MTI_vtx, MTI_tex, colorvec] = VtxTex_Clean(obj, MTI_vtx, MTI_tex, img)
           
            % Get the size of the color image
            [h, w] = size(rgb2gray(img));
            
            % Process and vertex and texture data
            idx_remove = (MTI_vtx(:,1) == 0) & (MTI_vtx(:,2) == 0) & (MTI_vtx(:,2) == 0);
            MTI_vtx = MTI_vtx(~idx_remove,:);
            MTI_tex = MTI_tex(~idx_remove,:);
            
            % Texture to pixel
            MTI_tex(:,1) = int16(MTI_tex(:,1) * w);
            MTI_tex(:,2) = int16(MTI_tex(:,2) * h);
            
            % Remove the minus index
            idx_remove_1 = MTI_tex(:,1) >= 640 | MTI_tex(:,1) <= 0 | MTI_tex(:,2) >= 480 | MTI_tex(:,2) <= 0;
            MTI_vtx = MTI_vtx(~idx_remove_1,:);
            MTI_tex = MTI_tex(~idx_remove_1,:);
            
            % Color vector 
            colorvec = [];
            for i = 1:length(MTI_tex(:,1))
                idx_width = MTI_tex(i,1);
                idx_height = MTI_tex(i,2);
                rgb = reshape(img(idx_height, idx_width, :), [1,3]);
                colorvec(i,1) = rgb(1); colorvec(i,2) = rgb(2); colorvec(i,3) = rgb(3);
            end
            colorvec = uint8(colorvec);            
        end
        function [target, colorvec_target] = Pixel2Voxel(obj, pixel_target, colorvec, MTI_tex, MTI_vtx)
            
            % colorvec is the input 
            % Initialization
            target = [];
            colorvec_target = [];
            r = 2;
            
            % For the 
            for i = 1:length(pixel_target)
                
                p1 = pixel_target(i,:) + [-r, -r];
                p2 = pixel_target(i,:) + [r, r];

                x1 = p1(1); y1 = p1(2);
                x2 = p2(1); y2 = p2(2);

                % 3D colorized point cloud
                idx_use = MTI_tex(:,1) >= y1 & MTI_tex(:,1) <= y2 & ...
                          MTI_tex(:,2) >= x1 & MTI_tex(:,2) <= x2;

                tex_use = MTI_tex(idx_use,:);
                vtx_use = MTI_vtx(idx_use,:);

                if length(vtx_use) == 0
                    continue;
                end

                color_use = colorvec(idx_use, :);
                target(i,:) = mean(vtx_use);
                colorvec_target(i,:) = color_use(1,:);
                
            end

        end
        function [] = Result(obj)
            
            load('planar.mat');
            load('nonplanar.mat');
            planar_pRMSE = []; 
            planar_nRMSE = []; 
            planar_pMAE = []; 
            planar_nMAE = []; 
            planar_in = [];
            planar_out = []; 

            nonplanar_pRMSE = []; 
            nonplanar_nRMSE = []; 
            nonplanar_pMAE = []; 
            nonplanar_nMAE = []; 
            nonplanar_in = [];
            nonplanar_out = []; 

            for i = 1 : length(CELL_planar)

                planar_pRMSE(i) = CELL_planar{i}.rmse_in;

                planar_in(i) = CELL_planar{i}.N_in;

                if isnan(CELL_planar{i}.rmse_out)
                    planar_nRMSE(i) = 0;
                    planar_nMAE(i) = 0;
                    planar_out(i) = 0;
                else
                    planar_nRMSE(i) = CELL_planar{i}.rmse_out;
                    planar_nMAE(i) = CELL_planar{i}.max_out;
                    planar_out(i) = CELL_planar{i}.N_out;
                end

                planar_pMAE(i) = CELL_planar{i}.max_in;

                nonplanar_pRMSE(i) = CELL_nonplanar{i}.rmse_in;
                nonplanar_nRMSE(i) = CELL_nonplanar{i}.rmse_out; 
                nonplanar_pMAE(i) = CELL_nonplanar{i}.max_in;
                nonplanar_nMAE(i) = CELL_nonplanar{i}.max_out;
                nonplanar_in(i) = CELL_nonplanar{i}.N_in;
                nonplanar_out(i) = CELL_nonplanar{i}.N_out;
            end 

            % plot the figure
            figure(1);clf
            plot(planar_pMAE, planar_nMAE, 'r*', 'MarkerSize', 10); 
            hold on;
            plot(nonplanar_pMAE, nonplanar_nMAE, 'bx', 'MarkerSize', 10);
            hold on;
            grid on; 
            xlim([0,2.5]); ylim([0,2.5]);
            ax = gca;
            ax.XGrid = 'on';
            ax.GridColor = [0, 0, 0];
            ax.XColor = [0 0 0]
            xlabel('pMAE (mm)', 'Fontsize', 30);
            ylabel('nMAE (mm)', 'Fontsize', 30);
            legend({'Planar surface', 'Cylindrical surface'}, 'FontSize',30);
            grid on; 
            xlim([0,2.5]); ylim([0,2.5]);
            axis([0 2.5 0 2.5]);
            axis manual; 

            figure(2);clf
            plot(planar_pRMSE, planar_nRMSE, 'r*', 'MarkerSize', 10); 
            hold on;
            plot(nonplanar_pRMSE, nonplanar_nRMSE, 'b*', 'MarkerSize', 10);
            grid on; 
            xlim([0,2.5]); ylim([0,2.5]);
            ax = gca;
            ax.GridColor = [0.1, 0.1, 0.1];
            xlabel(' pRMSE (mm)', 'Fontsize', 30);
            ylabel(' nRMSE (mm)', 'Fontsize', 30);
            legend({'Planar surface', 'Cylindrical surface'}, 'FontSize',30);
            grid on; 
            xlim([0,2.5]); 
            ylim([0,2.5]);
            axis([0 2.5 0 2.5]);
            axis manual; 

            % show the bar histogram 
            figure(3);clf
            bar(1:length(planar_in), planar_in ./ (planar_in + planar_out));
            hold on; 
            bar(1:length(planar_out), planar_out ./ (planar_in + planar_out));
            grid on;
            xlabel('Planar test case index (1~30)', 'Fontsize', 30);
            ylabel('Percentage of voxels', 'Fontsize', 30);
            legend({'positive error', 'negative error'}, 'FontSize',30);

            figure(4);clf 
            bar(1:length(nonplanar_in), nonplanar_in ./ (nonplanar_in + nonplanar_out));
            hold on; 
            bar(1:length(nonplanar_out), nonplanar_out ./ (nonplanar_in + nonplanar_out));
            grid on;
            xlabel('Nonplanar test case index (1~30)', 'Fontsize', 30);
            ylabel('Percentage of voxels', 'Fontsize', 30);
            legend({'positive error', 'negative error'}, 'FontSize',30);

        end
    end
end