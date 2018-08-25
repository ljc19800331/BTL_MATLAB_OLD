% class: sec_pre
% This class is designed for preprocessing the data
% input: Different data format
% output: Pointcloud data format
%% 
classdef sec_pre
    properties
        model_x = 3.346;                % The lengh of model measure in inch
        model_y = 4.547;
        model_z = 1.571;
        scanbox = [0.4, 0.5, 0.6];     % The size of scanbox
        scandis = [0.005, 0.01, 0.02]; % The resolution of scanbox
        X_Label = [ 0.95  0.95  0.95  0.95  0.95  2.35  2.35  2.35  2.35  2.35];    % The label of the position measured in inch
        Y_Label = [ 1.45  1.7   1.95  2.2   2.45  1.45  1.7   1.95  2.2   2.45];
        model_x_img = 669;
        model_y_img = 910;
    end
    methods
        function data = txt2pc(obj, txt_x, txt_y, txt_z)
            
            % Read the data
            scan_x = load(txt_x);
            scan_y = load(txt_y);
            scan_z = load(txt_z);
           
            % Assume txt_x, txt_y and txt_z have same dimensions
            pc = zeros(length(scan_x), 3);
            
            % Make zero 
            scan_x = scan_x - min(scan_x);
            scan_y = scan_y - min(scan_y);
            scan_z = scan_z - min(scan_z);
            
            % Convert data to point cloud 
            pc(:,1) = scan_x;
            pc(:,2) = scan_y;
            pc(:,3) = scan_z;
            data = pointCloud(pc);
            
        end   
        function [pc, data_xyz] = stl2pc(obj, stl_name)
            
            % pc: pointcloud
            % data_xyz: the numpy data of the pointcloud
            data_stl = stlread(stl_name);
            data_use = unique(data_stl.vertices(:,1:3),'rows');
            data_xyz = data_use; 
            data_xyz(:,1) = data_xyz(:,1) - min(data_xyz(:,1));
            data_xyz(:,2) = data_xyz(:,2) - min(data_xyz(:,2));
            data_xyz(:,3) = data_xyz(:,3) - min(data_xyz(:,3));
            pc = pointCloud(data_xyz);
            
        end
        function data = ply2pc(obj, ply_name)
            ptCloud = pcread(ply_name);
            pt_xyz = zeros(length( ptCloud.Location(:,1) ), 3);
            pt_xyz(:,1) = ptCloud.Location(:,1) - min(ptCloud.Location(:,1));
            pt_xyz(:,2) = ptCloud.Location(:,2) - min(ptCloud.Location(:,2));
            pt_xyz(:,3) = ptCloud.Location(:,3) - min(ptCloud.Location(:,3));
            data = pointCloud(pt_xyz);
        end
        function data = pcd2pc(obj, pcd_name)
            ptCloud = pcread(pcd_name);
            pt_xyz = zeros(length( ptCloud.Location(:,1) ), 3);
            pt_xyz(:,1) = ptCloud.Location(:,1) - min(ptCloud.Location(:,1));
            pt_xyz(:,2) = ptCloud.Location(:,2) - min(ptCloud.Location(:,2));
            pt_xyz(:,3) = ptCloud.Location(:,3) - min(ptCloud.Location(:,3));
            data = pointCloud(pt_xyz);
        end
        function data = pcLenChange(obj, pc, flag_type, flag_Len)
            % Convert pc to pc(with len unit change)
            % Load the data
            pc_xyz = pc.Location;
            x_range = max(pc_xyz(:,1)) - min(pc_xyz(:,1));
            y_range = max(pc_xyz(:,2)) - min(pc_xyz(:,2));
            z_range = max(pc_xyz(:,3)) - min(pc_xyz(:,3));
            
            % Conver to the "inch" or "cm" unit
            if flag_type == "brain"
                if flag_Len == "cm"
                   x_len = obj.model_x * 2.54;
                   y_len = obj.model_y * 2.54;
                   z_len = obj.model_z * 2.54;
                   c_x = x_len./x_range;        % cm/N
                   c_y = y_len./y_range;
                   c_z = z_len./z_range;
                   pc_xyz(:,1) = pc_xyz(:,1).* c_x;
                   pc_xyz(:,2) = pc_xyz(:,2).* c_y;
                   pc_xyz(:,3) = pc_xyz(:,3).* c_z;
                end

                if flag_Len == "inch"
                   x_len = obj.model_x;
                   y_len = obj.model_y;
                   z_len = obj.model_z;
                   c_x = x_len./x_range;         % cm/N
                   c_y = y_len./y_range;
                   c_z = z_len./z_range;
                   pc_xyz(:,1) = pc_xyz(:,1).* c_x;
                   pc_xyz(:,2) = pc_xyz(:,2).* c_y;
                   pc_xyz(:,3) = pc_xyz(:,3).* c_z;
                end
            end
            
            if flag_type == "scan"
                if flag_Len == "cm"
                   c = 2.54; % cm/inch
                   pc_xyz(:,1) = pc_xyz(:,1).* c;
                   pc_xyz(:,2) = pc_xyz(:,2).* c;
                   pc_xyz(:,3) = pc_xyz(:,3).* c;
                end

                if flag_Len == "inch"
                   c = 1; % inch/inch
                   pc_xyz(:,1) = pc_xyz(:,1).* c;
                   pc_xyz(:,2) = pc_xyz(:,2).* c;
                   pc_xyz(:,3) = pc_xyz(:,3).* c;
                end
            end
            
            data = pointCloud(pc_xyz);
        end
        function data = BrainInit(obj, pc_brain)
            
            % Move the brain model to the expected view
            % In theory, it is not necessary
            theta_x = pi; % 180
            R_x = [1 0 0 0; ...
                  0 cos(theta_x) sin(theta_x) 0; ...
                  0 -sin(theta_x) cos(theta_x) 0; ...
                  0 0 0 1];
            theta_y = pi; % 180
            R_y = [cos(theta_y) 0 -sin(theta_y) 0; ...
                  0 1 0 0; ...
                  sin(theta_y) 0 cos(theta_y) 0; ...
                  0 0 0 1];
            theta_z = 0;  % 180
            R_z = [cos(theta_z) sin(theta_z) 0 0; ...
                 -sin(theta_z) cos(theta_z) 0 0; ...
                  0 0 1 0; ...
                  0 0 0 1];
            tform = affine3d(R_x*R_y*R_z);
            pc_brain_1 = pctransform(pc_brain, tform);
            
            % Normalize the data
            pc_brain_2 = [];
            pc_brain_2(:,1) = pc_brain_1.Location(:,1) - min(pc_brain_1.Location(:,1));
            pc_brain_2(:,2) = pc_brain_1.Location(:,2) - min(pc_brain_1.Location(:,2));
            pc_brain_2(:,3) = pc_brain_1.Location(:,3) - min(pc_brain_1.Location(:,3));  
            data = pointCloud(pc_brain_2);
            
        end   
        function data = Brain2Img(obj, pc_brain, n_x, n_y)
            
            % Load the data
            x_obj = pc_brain.Location(:,1);
            y_obj = pc_brain.Location(:,2);
            v_obj = pc_brain.Location(:,3);
            
            % n_x: Number of pixels in x axis for brain img
            % n_y: Number of pixels in y axis for scan img
            [xq,yq] = meshgrid(min(x_obj):(max(x_obj)-min(x_obj))/n_x:max(x_obj), ...
                       min(y_obj):max((y_obj)-min(y_obj))/n_y:max(y_obj));
            vq = griddata(x_obj, y_obj, v_obj, xq, yq, 'natural');
            
            % flp the image for other applications 
            data = -flipud(vq);
            data(isnan(data)) = 0; 
            
        end
        function data = Scan2Img(obj, pc_scan, pixel_range)
            
            mat_xyz = []; 
            % mat_xyz: The matrix to save all the x,y,z data
            x = pc_scan.Location(:,1); mat_xyz(:,1) = x;
            y = pc_scan.Location(:,2); mat_xyz(:,2) = y;
            z = pc_scan.Location(:,3); mat_xyz(:,3) = z;
            
            %% Consider x as the signal and find the peaks
            
            t = 1:length(x); 
            % Upper peak
            u_middle = (max(x)-min(x))/2 + min(x);
            [pks_u, locs_u] = findpeaks(x,'Threshold',0,'MinPeakDistance',pixel_range,'MinPeakHeight',u_middle);
            
            % Lower peak
            l_middle = (max(-x)-min(-x))/2 + min(-x);
            [pks_l, locs_l] = findpeaks(-x, 'Threshold',0,'MinPeakDistance',pixel_range,'MinPeakHeight',l_middle);

            % Check if the first value is correct
            if locs_l(1) < locs_u(1)
               locs_l(1) = []; 
            end
            
            % Save the data into a matrix
            n_length = max( length(locs_l),length(locs_u));
            peak_locs = zeros( n_length, 2 ); 
            peak_locs(1:length(locs_u),1) = locs_u; peak_locs(1:length(locs_l),2) = locs_l;  

            % The last number is the length of the data
            peak_locs(max( length(locs_l),length(locs_u)),2) = length(x);

            %Get the peak vector
            peak_vec = zeros( n_length *2, 1);
            peak_vec( 1 : length(locs_u)) = locs_u;
            peak_vec( length(locs_u)+1 : length(locs_u)+length(locs_l)) = locs_l;
            if peak_vec( n_length*2,1) == 0 
                peak_vec( n_length*2,1) = length(x);
            end
            
            peak_vec = sort(peak_vec);
            
            % Crop the data and save in an image
            img_mat = zeros(n_length,n_length);
            n1 = 1;     %The first number, by default
            flag = 1;
    
            for i = 1: length(peak_vec)
                n2 = peak_vec(i);
                vec_z = z(n1:n2);
                % different flag value represents different direction for the raster
                % scan lines
                if mod(flag,2) == 1
                if length(vec_z) < 2*n_length
                    img_mat( (2*n_length+1)-i ,(2*n_length + 1) - (1:length(vec_z)) ) = vec_z(1:length(vec_z));
                    img_mat( (2*n_length+1)-i ,(2*n_length + 1) - (length(vec_z):length(vec_z)+1)) = 0;
                else
                img_mat( (2*n_length+1)-i , (2*n_length + 1) - (1:2*n_length) ) = vec_z(1:2*n_length);
                end
                n1 = peak_vec(i);
                end

                if mod(flag,2) == 0
                if length(vec_z) < 2*n_length
                     img_mat( (2*n_length+1)-i , 1:length(vec_z) ) = vec_z(1:length(vec_z));
                else
                img_mat( (2*n_length+1)-i , 1:2*n_length ) = vec_z(1:2*n_length);
                end
                n1 = peak_vec(i);
                end
                flag = flag + 1;
            end
            
            data = img_mat;
            
        end
        function [] = ImgGen6D(obj, pc_brain, flag_dof, flag_box, flag_iter, Range_mat)
            
            % "D" here means the degree of freedom
            % "flag_dof" means the degree of freedom
            % pc_brain: The original brain model after normalization and
            % Zeros 
            % flag_box = 2;
            % flag_iter = 30; 
            % Range_mat = [0, 0, -45, 0, 0, 0];
            
            if flag_dof == "6D"
                for flag_dis = 1:3 
                for flag_xrot = Range_mat(1):flag_iter:Range_mat(2)
                for flag_yrot = Range_mat(3):flag_iter:Range_mat(4)
                for flag_zrot = Range_mat(5):flag_iter:Range_mat(6)

                % check the file if it exists
                file_check = strcat('img_6D', '/', num2str(flag_xrot), ...
                                    '_', num2str(flag_yrot), ...
                                    '_', num2str(flag_zrot), '.png');
                if exist(file_check, 'file') == 2
                    continue;
                end

                % determine the rotation angle 
                theta_x = pi./180 * flag_xrot;
                R_x = [1 0 0 0; ...
                      0 cos(theta_x) sin(theta_x) 0; ...
                      0 -sin(theta_x) cos(theta_x) 0; ...
                      0 0 0 1];
                theta_y = pi./180 * flag_yrot;
                R_y = [cos(theta_y) 0 -sin(theta_y) 0; ...
                      0 1 0 0; ...
                      sin(theta_y) 0 cos(theta_y) 0; ...
                      0 0 0 1]; 
                theta_z = pi./180 * flag_zrot;
                R_z = [cos(theta_z) sin(theta_z) 0 0; ...
                     -sin(theta_z) cos(theta_z) 0 0; ...
                      0 0 1 0; ...
                      0 0 0 1];

                % perform the mapping
                pc_use = pc_brain * R_x(1:3,1:3) * R_y(1:3,1:3) * R_z(1:3,1:3);
                pc_use(:,1) = pc_use(:,1) - min(pc_use(:,1));
                pc_use(:,2) = pc_use(:,2) - min(pc_use(:,2));
                pc_use(:,3) = pc_use(:,3) - min(pc_use(:,3));

                % hpr operator
                % design the parameter as 3 higher then the
                % value
                Campos_x =  min(pc_use(:,1)) + (max(pc_use(:,1) -  min(pc_use(:,1))))/2;
                Campos_y =  min(pc_use(:,2)) + (max(pc_use(:,2) -  min(pc_use(:,2))))/2;
                Campos_z = 3 * max(pc_use(:,3)); 
                param = 3.5;
                C = [Campos_x, Campos_y, Campos_z];
                p = pc_use;
                visiblePtInds = HPR(p,C,param);
                pc_hpr = pc_use(visiblePtInds,:);

                % generate the image
                pixel_scan = obj.scanbox(flag_box) / obj.scandis(flag_dis);
                x_box = obj.scanbox(flag_box); y_box = obj.scanbox(flag_box); % the lenght of scan region
                x_brain = max(pc_use(:,1)); y_brain = max(pc_use(:,2));
                n_x = round(pixel_scan * x_brain / x_box); n_y = round(pixel_scan * y_brain / y_box);

                x_obj = pc_hpr(:,1);
                y_obj = pc_hpr(:,2);
                v_obj = pc_hpr(:,3);

                [xq,yq] = meshgrid(min(x_obj):(max(x_obj)-min(x_obj))/n_x:max(x_obj), ...
                               min(y_obj):max((y_obj)-min(y_obj))/n_y:max(y_obj));
                vq = griddata(x_obj,y_obj,v_obj,xq,yq, 'natural');
                img_hpr = -flipud(vq);
                img_hpr(isnan(img_hpr)) = 0;  

                savepath = strcat('img_6D', '/', num2str(flag_xrot), '_', num2str(flag_yrot), '_', num2str(flag_zrot), '.png');
                imwrite(mat2gray(img_hpr), savepath);

                end
                end
                end
                end
            end
        end
        function [] = BackCulling(obj, stlname)
            
            % F and V are both the vertices and faces
            [F,V] = stlread(stlname);
            N = normals(V,F);
            BC = barycenter(V,F);
            
            % Define the backface culling processing
            back_facing = sum(N.*bsxfun(@minus,BC,campos),2)<=0;

            t = tsurf(F,V,'EdgeColor','none','FaceLighting','phong');
            view(2);
            axis equal;
            camproj('persp');
            t.FaceVertexCData = 1*(sum(N.*bsxfun(@minus,BC,campos),2)<=0);
            t.FaceVertexCData(sum(N.*bsxfun(@minus,BC,campos),2)>0) = nan;
            test = 1*(sum(N.*bsxfun(@minus,BC,campos),2)<=0);
            rotate3d on;
            apply_ambient_occlusion();
            
        end
    end    
end
    