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
        function [] = Err_pixel(obj)
            
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
            [MTI_vtx_before, MTI_tex_before, colorvec_before] = VtxTex_Clean(MTI_vtx_before, MTI_tex_before, img_before); 
            [MTI_vtx_after, MTI_tex_after, colorvec_after] = VtxTex_Clean(MTI_vtx_after, MTI_tex_after, img_after);
        
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
            
            
        end 
        
        function [MTI_vtx, MTI_tex, colorvec] = VtxTex_Clean(MTI_vtx, MTI_tex, img)
           
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
        function [] = Pixel2Voxel()
            
        end
    end
end






