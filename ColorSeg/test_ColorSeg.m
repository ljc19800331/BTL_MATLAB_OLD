
figure(1);

for flag_img = 1:5

flag_img
  
% load the images
subplot(2,3,flag_img); 
test = sec_ColorSeg();
test.imgpath = strcat('IMG_378',num2str(flag_img-1),'.JPG');

% Define the number of color levels
test.NumColors = 10;
Lout = test.method_1(); % the output label matrix of the whole image

% segment the regions
Res = [];

for idx_test = 1: test.NumColors
    
img_obj = imread(test.imgpath);
img_gray = rgb2gray(img_obj);
img_gray(Lout==idx_test) = 255;
img_gray(Lout~=idx_test) = 0;

% Suppose the number of phantom is either 1 or 2 
% This still needs to be hard code

img_use = img_gray;
img_bin = imbinarize(img_use);
L = bwlabel(img_bin);

% if it is only 1 or 2 region, the 3 region should be empty
[r, c] = find(L==3);
if length(r) == 0
    Res(idx_test) = 1;
else
    Res(idx_test) = 0; 
end
end 

% If two similar results happen, choose the larger size region
if length(find(Res==1)) > 1
    idx_1 = find(Res == 1 );
    % compare the size
    Res_1 = [];
    for i =  1: length(idx_1)
        mat = (Lout == idx_1(i));
        Res_1(i) = length(find(mat==1));
    end
    [value_use, idx_use] = max(Res_1);
    idx_out = idx_1(idx_use);
else
    idx_out = find(Res == 1 );
end

% Based on the output index, choose the binary region
    img_gray = rgb2gray(img_obj);
    img_gray(Lout==idx_out) = 255;
    img_gray(Lout~=idx_out) = 0;
    img_use = img_gray;
    imshow(img_use);
    
end
