% This is the sec_map class 
% This class includes the following function
% CoarseRgs: Coarse registration
% FineRgs: Fine registration
% VizRgs: Visualize the registration result

classdef sec_map 
    properties
        model_x = 3.346;                % The lengh of model measure in inch
        model_y = 4.547;
        model_z = 1.571;
        scanbox = [0.4, 0.5, 0.6];      % The size of scanbox
        scandis = [0.005, 0.01, 0.02];  % The resolution of scanbox
        X_Label = [ 0.95  0.95  0.95  0.95  0.95  2.35  2.35  2.35  2.35  2.35];    % The label of the position measured in inch
        Y_Label = [ 1.45  1.7   1.95  2.2   2.45  1.45  1.7   1.95  2.2   2.45];
        model_x_img = 669;
        model_y_img = 910;
    end
    methods
        function [Result] = CoarseRgs(obj, img_scan, Range_mat, flag_iter)
            
            % Range_mat is the degree range of the search angles
            % [xmin, xmax, ymin, ymax, zmin, zmax];
            % Range_mat = [0, 0, -45, 0, 0, 0];
            % flag_iter = 5;
            flag_box = 2;
            
            X_rot = Range_mat(1):flag_iter:Range_mat(2);
            Y_rot = Range_mat(3):flag_iter:Range_mat(4);
            Z_rot = Range_mat(5):flag_iter:Range_mat(6);
            
            % input: The brain model
            % output: The center of the laser in the brain image
            C_max = {};     % The max coefficients
            P_center = {};  % The coordinates corresponding to the max coefficients
            count = 0;
            count_x = 0;
            
          %  for flag_xrot = Range_mat(1):flag_iter:Range_mat(2)
          Result = {};
          parfor count_x = 1: length(X_rot)
              flag_xrot = X_rot(count_x);
              [Res] = Loop(flag_xrot, Y_rot, Z_rot, img_scan);
              Result{count_x} = Res;
          end
        end
        function [R_output, T_output, pc_obj_2, pc_two, Agl_mat, C_max] = FineRgs(obj, pc_brain, pc_scan, Result, Range_mat, flag_box, flag_iter)

            X_rot = Range_mat(1):flag_iter:Range_mat(2);
            Y_rot = Range_mat(3):flag_iter:Range_mat(4);
            Z_rot = Range_mat(5):flag_iter:Range_mat(6);
            
            % Decode the data
            [a,b] = size(Result);
            C_max = [];
            for i = 1: length(Result)   % x
                res_1 = Result{i};
                [a,b] = size(res_1);
                for j = 1:a     % y
                    for k = 1:b % z 
                        C_max(i,j,k) = res_1{j,k}(3); 
                    end
                end
            end
            
            % Here: C_max is a 3D array 
            % D1: xrot, D2, yrot, D3, zrot
            for k=1:size(C_max,3)
               a=C_max(:,:,k);
               [maxv,idx]=max(a(:));
               [ii,jj]=ind2sub([size(C_max,1) size(C_max,2)],idx);
               out{k}=[maxv ii jj];
            end
            max_cx = [];
            for k = 1: length(out)
                max_cx(k) = out{k}(1);
            end
            [max_v, idx_x] = max(max_cx);
            idx_y = out{idx_x}(3); idx_z = out{idx_x}(2); 
            % read xrot, yrot, zrot
            xrot = X_rot(idx_x);
            yrot = Y_rot(idx_y);
            zrot = Z_rot(idx_z);
            Agl_mat = [xrot, yrot, zrot];
            p_center = [Result{idx_x}{idx_y,idx_z}(1), Result{idx_x}{idx_y,idx_z}(2)]; 
      
            % Get the target brain model(with angle) 
            theta_x = pi./180 * xrot;
            R_x = [1 0 0 0; ...
                   0 cos(theta_x) sin(theta_x) 0; ...
                   0 -sin(theta_x) cos(theta_x) 0; ...
                   0 0 0 1];
            theta_y = pi./180 * yrot;
            R_y = [cos(theta_y) 0 -sin(theta_y) 0; ...
                   0 1 0 0; ...
                   sin(theta_y) 0 cos(theta_y) 0; ...
                   0 0 0 1]; 
            theta_z = pi./180 * zrot;
            R_z = [cos(theta_z) sin(theta_z) 0 0; ...
                   -sin(theta_z) cos(theta_z) 0 0; ...
                   0 0 1 0; ...
                   0 0 0 1];
            
            % HPR -- The HPR operator
            pc_use = pc_brain.Location * R_x(1:3,1:3) * R_y(1:3,1:3) * R_z(1:3,1:3);
            
            pc_use(:,1) = pc_use(:,1) - min(pc_use(:,1));
            pc_use(:,2) = pc_use(:,2) - min(pc_use(:,2));
            pc_use(:,3) = pc_use(:,3) - min(pc_use(:,3));
            
            Campos_x =  min(pc_use(:,1)) + (max(pc_use(:,1) -  min(pc_use(:,1))))/2;
            Campos_y =  min(pc_use(:,2)) + (max(pc_use(:,2) -  min(pc_use(:,2))))/2;
            Campos_z = 3 * max(pc_use(:,3)); 
            
            param = 3.5;
            C = [Campos_x, Campos_y, Campos_z];
            p = pc_use;
            visiblePtInds = HPR(p,C,param);
            pc_hpr = pc_use(visiblePtInds,:);  
               
            % Corresponding image size
            imgpath = strcat( 'img_6D/', num2str(xrot), '_', num2str(yrot), '_', num2str(zrot), '.png'); 
            img_brain = imread(imgpath);
            [img_l, img_w] = size(img_brain);
            
            % The center coordinates
            x_crop = p_center(1) / img_w * max(pc_use(:,1)); 
            y_crop = (obj.model_y - p_center(2) / img_l * max(pc_use(:,2))  );
            
            % crop the 2 size point cloud region
            p_center = [x_crop, y_crop];
            range_two = obj.scanbox(flag_box)*1.25;
            range_one = obj.scanbox(flag_box)*0.75;
            idx_two = pc_hpr (:,1) > (p_center(1)-range_two) & pc_hpr (:,1) < (p_center(1) + range_two) ...
            & pc_hpr (:,2) > (p_center(2)-range_two) & pc_hpr (:,2) < (p_center(2)+range_two);
            idx_one = pc_hpr (:,1) > (p_center(1)-range_one) & pc_hpr (:,1) < (p_center(1) + range_one) ...
            & pc_hpr (:,2) > (p_center(2)-range_one) & pc_hpr (:,2) < (p_center(2)+range_one);
            pc_two = pc_hpr(idx_two,:);
            pc_one = pc_hpr(idx_one,:);   
            
            % So far, we get pc_one, pc_two, pc_use(pc_brain), pc_scan
            % ICP for fine alignment
            
            % Segment the region pointcloud
            pc_obj_1 = sec_GSseg(pc_two);
            pc_obj_2 = sec_GSseg(pc_scan.Location);
            
            % Change the view of the scan point cloud
            theta_y = -pi;               % 180 degree
            R_y = [cos(theta_y) 0 -sin(theta_y) 0; ...
                  0 1 0 0; ...
                  sin(theta_y) 0 cos(theta_y) 0; ...
                  0 0 0 1];
            pc_obj_2 = pc_obj_2 * R_y(1:3,1:3);
            pc_scan_output = pc_obj_2;
            
            % normalized again 
            pc_obj_2(:,1) = pc_obj_2(:,1) - min(pc_obj_2(:,1));
            pc_obj_2(:,2) = pc_obj_2(:,2) - min(pc_obj_2(:,2));
            pc_obj_2(:,3) = pc_obj_2(:,3) - min(pc_obj_2(:,3));
            
            % Downsample the point cloud
            pc_obj_2 = pcdownsample(pointCloud(pc_obj_2),'gridAverage', 0.025);   % Based on the property of the brain model
            pc_obj_2 = pc_obj_2.Location;
            
            % Centroid(align) the pointcloud
            center_obj_1 = [max(pc_one(:,1)) - min(pc_one(:,1)), ...
                            max(pc_one(:,2)) - min(pc_one(:,2)), ...
                            max(pc_one(:,3)) - min(pc_one(:,3))]/2 + ...
                            [min(pc_one(:,1)), min(pc_one(:,2)), min(pc_one(:,3))];
            
            center_obj_2 = [max(pc_obj_2(:,1)) - min(pc_obj_2(:,1)), ...
                            max(pc_obj_2(:,2)) - min(pc_obj_2(:,2)), ...
                            max(pc_obj_2(:,3)) - min(pc_obj_2(:,3))]/2 + ...
                            [min(pc_obj_2(:,1)), min(pc_obj_2(:,2)), min(pc_obj_2(:,3))];
            
            T_cen = repmat(center_obj_1 - center_obj_2, length(pc_obj_2), 1);
            pc_firstpos = pc_obj_2 + T_cen;
            
            % optimize the result based on ICP method
            [tform, movingReg] = pcregrigid(pointCloud(pc_firstpos), pointCloud(pc_obj_1), 'MaxIterations', 100, 'Tolerance', [0.0001, 0.00009]);
            R_rot = tform.T(1:3,1:3);
            T_rot = repmat(tform.T(4,1:3), length(pc_obj_2),1);
                        
            % The final position of the two point cloud
            pc_finalpos = (pc_obj_2 + T_cen) * R_rot + T_rot;
            R_output = R_rot;
            T_output = T_cen * R_rot + T_rot;
            
        end
        function [X,Y,Z] = VizRgs(obj, C_max)
            [id_x, id_y, id_z] = size(C_max);
            [X,Y,Z] = ndgrid(1:size(C_max,1), 1:size(C_max,2), 1:size(C_max,3));
            pointsize = 30;
            figure;clf;
            scatter3(X(:), Y(:), Z(:), pointsize, C_max(:), 'MarkerFaceColor',[0 1 1]);
            rotate3d on; hold on;
            [idx] = find(C_max > 0.92);
            scatter3(X(idx), Y(idx), Z(idx), pointsize, C_max(idx), 'MarkerFaceColor',[0.5 0.5 0.5]);
        end 
    end
end
