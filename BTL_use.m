%% BTL_ICRA project
% 1. ANN_CAlibration -- finished
% 2. MapStereoMTI -- finished
% 3. MapStereoMTI_Interpolation -- pending
% 4. TestDesign: Design the test data
% 5. FastPc
% 6. ErrorMap
% 7. sec_GSsegSurf 

classdef BTL_ICRA
    
    properties
        MTI2D_NPY_name = 'P_MTI_Calibration_interpolation_2D.npy';
        MTI3D_NPY_name = 'P_MTI_Calibration_interpolation_3D.npy';
        STEREO3D_NPY_name = 'P_STEREO_Calibration_3D.npy';
        Template_img_name = 'Aug27.jpg';
        Texture_name = 'MTI_textures.npy';
        Vertex_name = 'MTI_vertices.npy'
        img_name = 'Aug27.jpg';
        idx_miss = [1816, 1936, 1937, 1965, 1966, 1967, 1968, 2083, 2084, 2085, 2086, 2087, 2115, 2116, 2117, 2118, 2119, 2233, 2234, 2235, 2236, 2237, 2266, 2267, 2268, 2269, 2384, 2385, 2386, 2417];
    end
    methods 
        function [net, P_MTI_2d, P_MTI_3d] = Calibration_ANN(obj)
            
            % Read the data
            fprintf("Read the data from Python \n");
            img_MTI = imread(obj.Template_img_name);
            P_MTI_2d = readNPY(obj.MTI2D_NPY_name);
            P_MTI_3d = readNPY(obj.MTI3D_NPY_name);
            
            % preprocessing
            fprintf("The 2D and 3D array should be the same \n");
            P_MTI_2d = double(P_MTI_2d);
            P_MTI_3d(obj.idx_miss,:) = [];
            P_MTI_3d = P_MTI_3d(1: length(P_MTI_2d),:);
            
            % Train the model
            fprintf("Begin to train the model \n");
            net = fitnet(4); % NUMBER OF NEURONS IN HIDDEN LAYER
            [net tr y e] = train(net, P_MTI_2d', P_MTI_3d');
           
        end 
        function [R_final, t_final, P_MTI_tform, P_Stereo_tform] = Calibration_3D23D(obj)
            
            % Read the data
            P_MTI_3d = readNPY('P_MIT_Calibration_3D.npy');
            P_STEREO_3d = readNPY('P_STEREO_Calibration_3D.npy');

            % Preprocessing the data
            P_MTI_3d = P_MTI_3d .* 2.54;        % inch to cm
            P_STEREO_3d = P_STEREO_3d .* 100;   % mm to cm
            idx_remove = P_STEREO_3d(:,1) == 0 & P_STEREO_3d(:,2) == 0 & P_STEREO_3d(:,3) == 0; 
            P_STEREO_3d = P_STEREO_3d(~idx_remove, :);
            
            % Region segmentation
            flag_seg = input('Do you want to calculate the centroids? (1/0)');
            if(flag_seg == 1)
                [data_seg, Cell_regions] = sec_GSsegSurf(P_STEREO_3d);
                P_STEREO_3d_Seg = [];
                for i = 1: length(Cell_regions)
                    data_3d = Cell_regions{i};
                    centroid = mean(data_3d);
                    P_STEREO_3d_Seg(i,:) = centroid;
                end 
            end
           
            % Perform the mapping
            [R_final, t_final, P_MTI_tform, P_Stereo_tform] = STEREO_MTI_MAP(P_MTI_3d, P_STEREO_3d_Seg);

        end
        function [data_use, Cell_2d, Cell_3d, P_MTI_2d, P_MTI_3d] = Calibration_Interpolation(obj, X_laser, Y_laser)
            
            % Read the data
            P_MTI_2d = readNPY(obj.MTI2D_NPY_name);
            P_MTI_3d = readNPY(obj.MTI3D_NPY_name);
            
            % Preprocessing the data
            fprintf("The 2D and 3D array should be the same \n");
            P_MTI_2d = double(P_MTI_2d);
            P_MTI_3d(obj.idx_miss,:) = [];
            P_MTI_3d = P_MTI_3d(1: length(P_MTI_2d),:);
            
            % interpolation process
            data_use = [];             
            Cell_2d = {};          
            Cell_3d = {};
            
            N_near = 20;     % The number of nearby pixels for KNN
            
            for i = 1 : length(X_laser) 
                i
                % find the ith pixel
                pixel = [X_laser(i), Y_laser(i)];
                % The N nearby pixels
                [idx_2d, value_2d] = knnsearch(P_MTI_2d, pixel, 'k', N_near, 'distance', 'euclidean');
                % Nearby pixel coordinates -- 2D
                pixel_2d = P_MTI_2d(idx_2d,:);
                Cell_2d{i} = pixel_2d;
                % Nearby point coordinates -- 3D
                idx_3d = idx_2d;
                point_3d = P_MTI_3d(idx_3d,:);
                Cell_3d{i} = point_3d;
                % The mean of the nearby 3D points
                point_final = mean(point_3d);
                % Save the data
                data_use(i,:) = point_final;
            end
            
        end
        function [pixel_test] = TestDesign(obj, flag_method)
            
            if(flag_method == 1)
                fprintf("Rectangular region based on pixels \n");
                p1 = [243, 196];
                p2 = [382, 288];
                X_laser = linspace(p1(1), p2(1), (p2(1) - p1(1) + 1));
                Y_laser = linspace(p1(2), p2(2), (p2(2) - p1(2) + 1));
                x_length = length(X_laser); 
                y_length = length(Y_laser);
                Y_laser = repmat(Y_laser, [x_length, 1] );
                Y_laser = Y_laser(:)';
                X_laser = repmat(X_laser, [1, y_length] );
                pixel_test = [X_laser',Y_laser'];
            end
            
            if(flag_method == 2)
                fprintf("Random region based on manual drawing operation \n");
                [pixel_test, image_dilated] = process_tattoo_image(obj.img_name);
            end
        end
        function [] = FastPc(obj, p1, p2)

            % Read the data
            MTI_vtx = readNPY(obj.Texture_name);
            MTI_tex = readNPY(obj.Vertex_name);
            img = imread(obj.Template_img_name);
            x1 = p1(1); y1 = p1(2);
            x2 = p2(1); y2 = p2(2);
            
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
            for i = 1:length(MTI_tex(:,1))
                idx_width = MTI_tex(i,1);
                idx_height = MTI_tex(i,2);
                rgb = reshape(img(idx_height, idx_width, :), [1,3]);
                colorvec(i,1) = rgb(1); colorvec(i,2) = rgb(2); colorvec(i,3) = rgb(3);
            end

            colorvec = uint8(colorvec);

            % 3D colorized point cloud
            idx_use = MTI_tex(:,1) >= x1 & MTI_tex(:,1) <= x2 & MTI_tex(:,2) >= y1 & MTI_tex(:,2) <= y2;
            MTI_tex = MTI_tex(idx_use,:);
            MTI_vtx = MTI_vtx(idx_use,:);
            colorvec = colorvec(idx_use, :);
            
        end
    end
end 



