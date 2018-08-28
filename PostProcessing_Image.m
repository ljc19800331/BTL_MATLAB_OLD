% Post processing
%% Load the image
img_before = imread('C:\Users\gm143\Documents\MATLAB\BTL\exp_3\before.png');
img_ANN = imread('C:\Users\gm143\Documents\MATLAB\BTL\exp_3\after_firstcut_interpolation.png');
img_interpolation = imread('C:\Users\gm143\Documents\MATLAB\BTL\exp_3\after_secondcut_interpolation.png'); 

%% Greyscale processing
img_before_grey = rgb2gray(img_before);
img_ANN_grey = rgb2gray(img_ANN);
img_interpolation_grey = rgb2gray(img_interpolation);

figure(1);clf;
subplot(1,3,1)
imshow(img_before);
subplot(1,3,2)
imshow(img_ANN); 
subplot(1,3,3) 
imshow(img_interpolation); 

%% Find the edge region
[pixel_before, image_before] = process_tattoo_image('C:\Users\gm143\Documents\MATLAB\BTL\exp_3\before.png');
[pixel_first, image_first] = process_tattoo_image('C:\Users\gm143\Documents\MATLAB\BTL\exp_3\after_firstcut_interpolation.png');
[pixel_second, image_second] = process_tattoo_image('C:\Users\gm143\Documents\MATLAB\BTL\exp_3\after_secondcut_interpolation.png');

%% Adjust the offset 
% pixel_ANN(:,1) = pixel_ANN(:,1) + 4;
% pixel_Interpolation(:,1) = pixel_Interpolation(:,1) + 4;

%% Viz the final results 
[h, w, d] = size(image_before);
img_greytest = zeros(h ,w);
img_before_greytest = zeros(h ,w);
img_first_greytest = zeros(h ,w);
img_second_greytest = zeros(h ,w);

for i = 1 : length(pixel_before)
    img_greytest(pixel_before(i,1), pixel_before(i,2)) = 1;
end

for i = 1 : length(pixel_first)
    img_first_greytest(pixel_first(i,1), pixel_first(i,2)) = 1;
end

for i = 1 : length(pixel_second)
    img_second_greytest(pixel_second(i,1), pixel_second(i,2)) = 1;
end

figure(1);clf;
subplot(1,3,1)
imshow(image_before); 
subplot(1,3,2) 
imshow(img_first_greytest); 
subplot(1,3,3)
imshow(img_second_greytest); 

%% Measure the error
% define the area first
N_region_1 = length(find(image_before==1));
N_region_2 = length(find(image_first==1));
N_region_3 = length(find(image_second==1));

% determine the mismatch pixels
miss_12 = find((img_greytest == img_first_greytest)==0);
Npixel_12 = length(find((img_greytest == img_first_greytest)==0));

miss_13 = find((img_greytest == img_second_greytest)==0);
Npixel_13 = length(find((img_greytest == img_second_greytest)==0));

% determine the mismatched pixels
image_before(miss_12) = 255;
img_first_greytest(miss_12) = 255;

image_before(miss_13) = 255;
img_second_greytest(miss_13) = 255;

% Viz the missing pixels
figure(2);clf;
subplot(1,3,1)
imshow(image_before, []); 
subplot(1,3,2) 
imshow(img_first_greytest, []); 
subplot(1,3,3)
imshow(img_second_greytest, []); 

