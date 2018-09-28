classdef sec_thesis
    properties 
        model_x = 3.346;
        model_y = 4.547;
        model_z = 1.571;
        scanbox = [0.4, 0.5, 0.6];
        scandis = [0.005, 0.01, 0.02];
        X_Label = [ 0.95  0.95  0.95  0.95  0.95  2.35  2.35  2.35  2.35  2.35];
        Y_Label = [ 1.45  1.7   1.95  2.2   2.45  1.45  1.7   1.95  2.2   2.45];
    end
    methods
        function [RMSE_NOHPR, RMSE_HPR, AVG_RMSE_NOHPR, AVG_RMSE_HPR, RATE_NOHPR, RATE_HPR, rate_nohpr_case, rate_hpr_case] = Err_3D(obj)
            % This function is used for RMSE for perpendicular case 
            RMSE_NOHPR = [];
            RMSE_HPR = [];
            sum_rmse_nohpr = 0;
            sum_rmse_hpr = 0;
            rate_nohpr = 0;
            rate_hpr = 0;
            rate_nohpr_case = [];
            rate_hpr_case = [];
            for flag_case = 1:9 
       
            % Load the data -- no HPR
            load exp_Natural.mat;
            DATA = Data_errorg{flag_case}; % inch
            idx_use = find(DATA<0.1);
            DATA_nohpr = DATA(idx_use) * 25.4; 
            % Measure the RMSE for 9 cases -- output: 9 RMSE
            RMSE_nohpr = sqrt(sum(DATA_nohpr .^ 2) ./ length(DATA_nohpr));  % mm
            sum_rmse_nohpr = sum_rmse_nohpr + sum(DATA_nohpr .^ 2);
            rate_nohpr = rate_nohpr + length(DATA_nohpr);
            rate_nohpr_case(flag_case) = length(DATA_nohpr)./10;
            % Load the data -- HPR
            load exp_hprNoSeg.mat;
            DATA = Data_errorg{flag_case}; % inch
            idx_use = find(DATA<0.1);
            DATA_hpr = DATA(idx_use) * 25.4; 
            % Measure the RMSE for 9 cases -- output: 9 RMSE
            RMSE_hpr = sqrt(sum(DATA_hpr .^ 2) ./ length(DATA_hpr));  % mm
            sum_rmse_hpr = sum_rmse_hpr + sum(DATA_hpr .^ 2);
            rate_hpr = rate_hpr + length(DATA_hpr);
            RMSE_NOHPR(flag_case) = RMSE_nohpr; 
            RMSE_HPR(flag_case) = RMSE_hpr;
            rate_hpr_case(flag_case) = length(DATA_hpr)./10;
            end
            RATE_NOHPR = rate_nohpr ./ 90; 
            RATE_HPR = rate_hpr ./ 90;
            AVG_RMSE_NOHPR = sqrt(sum_rmse_nohpr./rate_nohpr);
            AVG_RMSE_HPR = sqrt(sum_rmse_hpr./rate_hpr); 

        end
        function [data] = Err_4D(obj)
           % For the case of 4DOF
           % Just show the demo
           a = 1;
        end
        function [data] = Err_5D(obj, CMAX) 
           % For the case of 5DOF
           % Just show the demo of the data with 5 DOF
           CMAX 
        end
        function [err, err_icp, err_icp2] = Err_6D(obj, flag_p, flag_angle, flag_viz, CMAX, Range_mat, pc_brain)
            % For the case of 6DOF
            % Find the max coefficients 
            % Get the brain model (viz)
            % Find the ground truth -- 3D points 
            % Find the registration result -- 3D points
            % 1-5 points, 5,15,25,35,45 for the angles -- 25 cases

            % Decode and get the coefficients Cmax 
            % CMAX = [p, agl, x, z, y];
         
            % Setting parameter 
%             flag_p = 2; % point index
%             flag_angle = 5; % angle index = 25 angle for y and 0 for x , z
            
            CMAX_xyz = CMAX(flag_p, flag_angle, :, :, :);
            CMAX_xyz = reshape(CMAX_xyz, [size(CMAX_xyz, 3), size(CMAX_xyz, 4), size(CMAX_xyz, 5)]);
            
            % max coeff
            % Range_mat = [-10, 50, -50, 50, -50, 50];
            flag_iter = 5;
            Xrot = Range_mat(1):flag_iter:Range_mat(2);
            Yrot = Range_mat(3):flag_iter:Range_mat(4);
            Zrot = Range_mat(5):flag_iter:Range_mat(6);
            
            ANGLE = [5,15,25,35,45];             
            Agl = ANGLE(flag_angle);
            test = CMAX_xyz;
            
            [max_C, i_C] = max((test(:)));
            [idx_x, idx_z, idx_y] = ind2sub([size(test,1) size(test,2), size(test,3)], i_C);
            cmax = test(idx_x, idx_z, idx_y); 
            
            % specific angle 
            Agl_x_rgs = Xrot(idx_x)
            Agl_y_rgs = -Yrot(idx_y)
            Agl_z_rgs = Zrot(idx_z)
            
            % Viz the scan and brain image and rgs result
                % img_brain
            imgpath = strcat('img_6D/', num2str(Agl_x_rgs), '_', num2str(Agl_y_rgs), '_', num2str(Agl_z_rgs), '.png');
            img_brain = imread(imgpath);
                % img_scan
            flag_dis = 1; 
            txt_x = strcat('p_rot/', 'P_rot_', num2str(flag_p), '/', num2str(Agl), '/', num2str(obj.scandis(flag_dis) * 1000), '/', 'Scan_x');
            txt_y = strcat('p_rot/', 'P_rot_', num2str(flag_p), '/', num2str(Agl), '/', num2str(obj.scandis(flag_dis) * 1000), '/', 'Scan_y');
            txt_z = strcat('p_rot/', 'P_rot_', num2str(flag_p), '/', num2str(Agl), '/', num2str(obj.scandis(flag_dis) * 1000), '/', 'Scan_z');
            
            pre = sec_pre();
            pc_test = pre.txt2pc(txt_x, txt_y, txt_z);
            pixel_range = 100;
            pc_scan = pre.txt2pc(txt_x, txt_y, txt_z); 
            img_scan = pre.Scan2Img(pc_scan, pixel_range);
%             image(img_scan, 'CDataMapping', 'scaled'); axis equal;
%             grid on; title('Test Scan2Img section'); 
            
            Img_1 = img_brain; Img_2 = img_scan; 
            Img_1(isnan(Img_1)) = 0;
            c = normxcorr2( Img_2, Img_1 ); 
            [max_c, imax] = max(abs(c(:)));   
            [ypeak, xpeak] = ind2sub(size(c),imax(1));
            x_center = xpeak - size(Img_2,2)/2;
            y_center = ypeak - size(Img_2,1)/2;
            
            if flag_viz == 1
            figure(1);
            subplot(2,3,1)
            imshow(img_brain); hold on;
            plot(x_center, y_center, 'rx');
            end
            img_brain_rgs = img_brain; 
            
            axis_forward = 1.198;                %inch
            axis_right = 1.872;         %inch
            axis_left = 3.346 - 1.872;  %inch
            
            pc_new = pc_brain.Location;
            pc_new(:,3) = pc_new(:,3) -axis_forward;
            pc_new(:,1) = pc_new(:,1) - axis_left;
            
            theta_x = pi./180 * Agl_x_rgs; % 180
            R_x = [1 0 0 0; ...
                  0 cos(theta_x) sin(theta_x) 0; ...
                  0 -sin(theta_x) cos(theta_x) 0; ...
                  0 0 0 1];
            theta_y = pi./180 * Agl_y_rgs; % 180
            R_y = [cos(theta_y) 0 -sin(theta_y) 0; ...
                  0 1 0 0; ...
                  sin(theta_y) 0 cos(theta_y) 0; ...
                  0 0 0 1];
            theta_z = pi./180 * Agl_z_rgs;  % 180
            R_z = [cos(theta_z) sin(theta_z) 0 0; ...
                 -sin(theta_z) cos(theta_z) 0 0; ...
                  0 0 1 0; ...
                  0 0 0 1];
            pc_brain_rgs = pc_new * R_x(1:3,1:3) * R_y(1:3,1:3) * R_z(1:3,1:3); 
            Range = 0.025;
            
            [img_l, img_w] = size(img_brain);
            
            move_x = x_center ./ img_w * (max(pc_brain_rgs(:,1) - min(pc_brain_rgs(:,1)))) + min(pc_brain_rgs(:,1)); % The distance between origin and the target x
            move_y = obj.model_y - y_center ./ img_l * obj.model_y; % The origin to the target y
            
            idx = pc_brain_rgs(:,1) > (move_x-Range) & ...
                  pc_brain_rgs(:,1) < (move_x+Range) & ...
                  pc_brain_rgs(:,2) > (move_y-Range) & ...
                  pc_brain_rgs(:,2) < (move_y+Range);
            pc_rgs = pc_brain_rgs(idx,:); 
            idx = pc_rgs(:,3) > mean(pc_rgs(:,3));
            pc_rgs = pc_rgs(idx,:);
            if flag_viz == 1
            figure(1);
            subplot(2,3,3)
            pcshow(pc_brain_rgs);hold on; pcshow(pc_rgs,'r');
            end
            pc_register = [mean(pc_rgs(:,1)), mean(pc_rgs(:,2)), mean(pc_rgs(:,3))];
            
            % Find the ground truth from the original model -- 3D 
            Agl_x_org = 0;
            Agl_y_org = -ANGLE(flag_angle);
            Agl_z_org = 0;
         
            theta_x = pi./180 * Agl_x_org; % 180
            R_x = [1 0 0 0; ...
                  0 cos(theta_x) sin(theta_x) 0; ...
                  0 -sin(theta_x) cos(theta_x) 0; ...
                  0 0 0 1];
            theta_y = pi./180 * Agl_y_org; % 180
            R_y = [cos(theta_y) 0 -sin(theta_y) 0; ...
                  0 1 0 0; ...
                  sin(theta_y) 0 cos(theta_y) 0; ...
                  0 0 0 1];
            theta_z = pi./180 * Agl_z_org;  % 180
            R_z = [cos(theta_z) sin(theta_z) 0 0; ...
                 -sin(theta_z) cos(theta_z) 0 0; ...
                  0 0 1 0; ...
                  0 0 0 1];
            pc_brain_org = pc_new * R_x(1:3,1:3) * R_y(1:3,1:3) * R_z(1:3,1:3); 
            
            % Find the closest point 
            Range = 0.025;
            move_x = obj.X_Label(flag_p+5); % The distance between origin and the target x
            move_y = obj.model_y - obj.Y_Label(flag_p+5); % The origin to the target y
            
            idx = pc_brain_org(:,1) > (-axis_left + move_x-Range) & ...
                  pc_brain_org(:,1) < (-axis_left + move_x+Range) & ...
                  pc_brain_org(:,2) > (move_y-Range) & ...
                  pc_brain_org(:,2) < (move_y+Range);
            pc_gt_org = pc_brain_org(idx,:); 
            
            idx = pc_gt_org(:,3) > mean(pc_gt_org(:,3));
            pc_gt_org = pc_gt_org(idx,:);
            if flag_viz == 1
            figure(1);
            subplot(2,3,2)
            pcshow(pc_brain_org);hold on; pcshow(pc_gt_org,'r');
            end
            pc_groundtruth = [mean(pc_gt_org(:,1)), mean(pc_gt_org(:,2)), mean(pc_gt_org(:,3))];
            
            err = sqrt(sum((pc_groundtruth - pc_register).^2)   ) .* 25.4; % mm
            
            % Perform the ICP registration and find the 3D points -- 3D 
                % Downsample
                pc_rgs_down = pcdownsample(pointCloud(pc_brain_rgs),'gridAverage', 0.025);
                pc_rgs_down = pc_rgs_down.Location;
                
                % HPR operator
                Campos_x =  min(pc_rgs_down(:,1)) + (  max(pc_rgs_down(:,1)) -  min(pc_rgs_down(:,1))  )/2;
                Campos_y =  min(pc_rgs_down(:,2)) + (  max(pc_rgs_down(:,2)) -  min(pc_rgs_down(:,2))  )/2;
                Campos_z = 3 * max(pc_rgs_down(:,3)); 
                param = 3.5;
                C = [Campos_x, Campos_y, Campos_z];
                p = pc_rgs_down;
                visiblePtInds = HPR(p,C,param);
                pc_hpr = pc_rgs_down(visiblePtInds,:);
                
                % Predict the position
                    % The center coordinates
                    range = 0.025; 
                    x_crop = x_center ./ img_w * (max(pc_brain_rgs(:,1) - min(pc_brain_rgs(:,1)))) + min(pc_brain_rgs(:,1));
                    y_crop = obj.model_y - y_center ./ img_l * obj.model_y;
                    % 
                    idx = pc_hpr(:,1) > (-axis_left + move_x-Range) & ...
                          pc_hpr(:,1) < (-axis_left + move_x+Range) & ...
                          pc_hpr(:,2) > (move_y-Range) & ...
                          pc_hpr(:,2) < (move_y+Range);
                    pc_hpr_center = pc_hpr(idx,:); 
                    if flag_viz == 1
                    figure(1);
                    subplot(2,3,4)
                    pcshow(pc_hpr); hold on; pcshow(pc_hpr_center, 'r');
                    end
                    pc_register_hpr = mean(pc_hpr_center);
                % Predict the scan position
                    pc_scan = pcdownsample(pc_scan,'gridAverage', 0.025);
                    pc_scan_use = pc_scan.Location;
                    pc_scan_use = sec_GSseg( pc_scan_use);
                    theta_y = -pi;               % 180 degree
                    R_y = [cos(theta_y) 0 -sin(theta_y) 0; ...
                          0 1 0 0; ...
                          sin(theta_y) 0 cos(theta_y) 0; ...
                          0 0 0 1];
                    pc_scan_use = pc_scan_use * R_y(1:3,1:3);
%                     pc_icp2 = pc_icp2 * R_y(1:3,1:3);
%                     pc_icp2(1) = pc_icp2(1) - min(pc_scan_use(:,1));
%                     pc_icp2(2) = pc_icp2(2) - min(pc_scan_use(:,2));
                    pc_scan_use(:,1) = pc_scan_use(:,1) - min(pc_scan_use(:,1));
                    pc_scan_use(:,2) = pc_scan_use(:,2) - min(pc_scan_use(:,2));
                    pc_scan_use(:,3) = pc_scan_use(:,3) - min(pc_scan_use(:,3));
                    
                    Range = 0.025;
                    x_center_icp2 = min(pc_scan_use(:,1)) + (max(pc_scan_use(:,1)) - min(pc_scan_use(:,1)))./2;
                    y_center_icp2 = min(pc_scan_use(:,2)) + (max(pc_scan_use(:,2)) - min(pc_scan_use(:,2)))./2;                    
                    idx_icp2 =  pc_scan_use(:,1) > (x_center_icp2-Range) & ...
                                pc_scan_use(:,1) < (x_center_icp2+Range) & ...
                                pc_scan_use(:,2) > (y_center_icp2-Range) & ...
                                pc_scan_use(:,2) < (y_center_icp2+Range);
                    pc_icp2 = pc_scan_use(idx_icp2, :);
                    idx = pc_icp2(:,3) > mean(pc_icp2(:,3));
                    pc_icp2 = mean(pc_icp2(idx,:));
                    
                    center_obj_1 = [max(pc_scan_use(:,1)) - min(pc_scan_use(:,1)), ...
                                    max(pc_scan_use(:,2)) - min(pc_scan_use(:,2)), ...
                                    max(pc_scan_use(:,3)) - min(pc_scan_use(:,3))]/2 + ...
                                    [min(pc_scan_use(:,1)), min(pc_scan_use(:,2)), min(pc_scan_use(:,3))];
                    center_obj_2 = pc_register_hpr; 
                    T_cen = repmat(center_obj_2 - center_obj_1, length(pc_scan_use), 1);
                    %T_cen_icp2 = repmat(center_obj_2 - center_obj_1, length(pc_icp2), 1);
                    pc_firstpos = pc_scan_use + T_cen;
                    pc_icp2 = pc_icp2 + center_obj_2 - center_obj_1;
                    
                    if flag_viz == 1
                    figure(1); subplot(2,3,5)
                    pcshow(pc_hpr); hold on; pcshow(pc_firstpos, 'r');hold on; pcshow(pc_icp2,'g');
                    end
                    % ICP 
                    % optimize the result based on ICP method
                    [tform, movingReg] = pcregrigid(pointCloud(pc_firstpos), pointCloud(pc_hpr), 'MaxIterations', 100, 'Tolerance', [0.0001, 0.00009]);
                    R_rot = tform.T(1:3,1:3);
                    T_rot = repmat(tform.T(4,1:3), length(pc_firstpos),1);
                    pc_icp2 = pc_icp2 * R_rot + tform.T(4,1:3);
                    pc_secondpos = pc_firstpos * R_rot + T_rot;
                    % viz the result
                  
                    % Find the new register center
                    range = 0.025; 
                    x_center_icp = min(pc_secondpos(:,1)) + (max(pc_secondpos(:,1)) - min(pc_secondpos(:,1)))./2;
                    y_center_icp = min(pc_secondpos(:,2)) + (max(pc_secondpos(:,2)) - min(pc_secondpos(:,2)))./2;
                    % 
                    idx = pc_secondpos(:,1) > (x_center_icp-Range) & ...
                          pc_secondpos(:,1) < (x_center_icp+Range) & ...
                          pc_secondpos(:,2) > (y_center_icp-Range) & ...
                          pc_secondpos(:,2) < (y_center_icp+Range);
                    pc_register_icp = pc_secondpos(idx,:); 
                    idx = pc_register_icp(:,3) > mean(pc_register_icp(:,3));
                    pc_register_icp = mean(pc_register_icp(idx,:));
                    if flag_viz == 1
                    figure(1);subplot(2,3,6);
                    pcshow(pc_hpr); hold on; pcshow(pc_secondpos,'r');
                    hold on; pcshow(pc_register_icp, 'b'); 
                    hold on; pcshow(pc_icp2, 'g'); 
                    hold on; plot3(pc_groundtruth(1), pc_groundtruth(2), pc_groundtruth(3), 'kx');
                    end
                    
                    if flag_viz == 1
                       figure(2);clf
                       pcshow(pc_brain_org); hold on; 
                       pcshow(pc_brain_rgs, 'r'); hold on;
                       pcshow(pc_groundtruth,'k');
%                        hold on; pcshow(pc_register, 'r'); 
%                        hold on; pcshow(pc_register_icp, 'b'); 
%                        hold on; pcshow(pc_icp2, 'g');
                    end
                    
                err_icp = sqrt(sum( (pc_groundtruth - pc_register_icp).^2 )) .* 25.4;    % mm
                err_icp2 = sqrt(sum( (pc_groundtruth - pc_icp2).^2  )) .* 25.4;  % mm
        end
        function [CMAX] = CMAX_MAT(obj, Res_Cmax)
            CMAX = [];
            for idx_p = 1: length(Res_Cmax)
                V_p = Res_Cmax{idx_p};
                for idx_agl = 1: length(V_p)
                    V_agl = V_p{idx_agl};
                    [L_x, L_z, L_y] = size(V_agl);
                    for idx_x = 1: L_x
                        for idx_z = 1: L_z
                            for idx_y = 1: L_y
                                data_cmax = V_agl(idx_x, idx_z, idx_y); 
                                CMAX(idx_p, idx_agl, idx_x, idx_z, idx_y) = data_cmax;
                            end
                        end
                    end
                end
            end 
        end
        function [] = VIZ_4D(obj, CMAX, Range_mat, flag_iter, flag_p, flag_angle)
            % Set up the value in the other two DOF
            Xrot = Range_mat(1):flag_iter:Range_mat(2);
            Yrot = Range_mat(3):flag_iter:Range_mat(4);
            Zrot = Range_mat(5):flag_iter:Range_mat(6);
            Agl_x = 0; idx_x = find(Xrot == 0); 
            Agl_z = 0; idx_z = find(Zrot == 0);
            % Find the value with Agl_x and Agl_z
            CMAX_xyz = CMAX(flag_p, flag_angle, :, :, :);
            CMAX_xyz = reshape(CMAX_xyz, [size(CMAX_xyz, 3), size(CMAX_xyz, 4), size(CMAX_xyz, 5)]);
            DATA_4D = CMAX_xyz(idx_x, idx_z, :);
            DATA_4D = DATA_4D(:);
            figure;clf;
            plot(DATA_4D, '-*');
            grid on; 
            xlabel('X angle', 'Fontsize', 20);
            ylabel('Max coefficients', 'Fontsize', 20);
        end
        function [] = VIZ_5D(obj, CMAX, Range_mat, flag_iter, flag_p, flag_angle)
            Xrot = Range_mat(1):flag_iter:Range_mat(2);
            Yrot = Range_mat(3):flag_iter:Range_mat(4);
            Zrot = Range_mat(5):flag_iter:Range_mat(6);
            Agl_x = 0; idx_x = find(Xrot == 0); 
            CMAX_xyz = CMAX(flag_p, flag_angle, :, :, :);
            CMAX_xyz = reshape(CMAX_xyz, [size(CMAX_xyz, 3), size(CMAX_xyz, 4), size(CMAX_xyz, 5)]);
            DATA_5D = CMAX_xyz(idx_x, :, :);
            DATA_5D = reshape(DATA_5D, [size(DATA_5D, 2), size(DATA_5D, 3)]);
            figure;clf;
            surf(DATA_5D), shading flat; rotate3d on;
            grid on;
            xlabel('X angle', 'Fontsize', 20);
            ylabel('Y angle', 'Fontsize', 20);
            zlabel('Max coefficients', 'Fontsize', 20);
        end
        function [] = VIZ_6D(obj, Res_Cmax, flag_p, flag_angle)
            figure;clf
            res_1 = Res_Cmax{flag_p}{flag_angle};
            [N_x, N_z, N_y] = size(res_1); 
            for flag_x = 1: N_x
                for flag_y = 1: N_y
                    for flag_z = 1: N_z
                        pos = [flag_x, flag_y, flag_z];
                        res = res_1(flag_x, flag_y, flag_z);
                        if res > 0.92 
                            plot3(pos(1), pos(2), pos(3), 'rx'); hold on; 
                        else
                            plot3(pos(1), pos(2), pos(3), 'bx'); hold on; 
                        end 
                    end 
                end
            end 
            grid on;
            xlabel('X angle', 'Fontsize', 20); 
            ylabel('Y angle', 'Fontsize', 20);
            zlabel('Z angle', 'Fontsize', 20);
            max_coff = max(max(max(res_1))); 
            rotate3d on;
        end 
    end
end