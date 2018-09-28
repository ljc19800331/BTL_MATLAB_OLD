% Problem definition
% Black object in the while background (well defined background)
% The objects are all connected regions
% Label the object from background

classdef sec_vessel
    properties
        imgpath = 'test_img1.png'; 
        NumColors = 2;
    end 
    methods
        function [Lout] = m1(obj)
            
            % Load the images
            img = imread(obj.imgpath);
            imshow(img); 
            
            % based on LAB algorithm
            img_lab = rgb2lab(img);
            [L, N] = superpixels(img_lab, 20000, 'isInputLab', true);
            BW = boundarymask(L);
            
            % convert the sp to cell array
            pixelIdxList = label2idx(L);  
            
            % Median color for each superpixel -- color feature
            meanColor = zeros(N,3);
            [m,n] = size(L);
            for  i = 1:N
                meanColor(i,1) = mean(img_lab(pixelIdxList{i}));
                meanColor(i,2) = mean(img_lab(pixelIdxList{i}+m*n));
                meanColor(i,3) = mean(img_lab(pixelIdxList{i}+2*m*n));
            end
            
            % Cluster the color feature using kmeans function
            numColors = obj.NumColors;
            [idx,cmap] = kmeans(meanColor,numColors,'replicates',2);
            cmap = lab2rgb(cmap);
            Lout = zeros(size(img_lab,1),size(img_lab,2));
            for i = 1:N
                Lout(pixelIdxList{i}) = idx(i);
            end
            
        end
    end
end
        




