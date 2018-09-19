function [xy_coordinates, img_output, Lout] = imgColorSeg(img_template, NumColors)

img_lab = rgb2lab(img_template);
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

numColors = NumColors;
[idx,cmap] = kmeans(meanColor,numColors,'replicates',2);
cmap = lab2rgb(cmap);

Lout = zeros(size(img_lab,1),size(img_lab,2));

for i = 1:N
    Lout(pixelIdxList{i}) = idx(i);
end

figure(1);clf 
imshow(Lout,[]);

%% Segment the color region
N_target = input('Which index is the target?');

idx_target = find(Lout == N_target);

mask_out = zeros(size(img_lab,1),size(img_lab,2));
mask_out(idx_target) = 1;

figure(1);clf 
imshow(mask_out);

imwrite(mask_out, 'mask_before.png');

filterImage = mask_out;

disp('Clean image.');
h1 = figure; title('Clean the image...')
while(ishandle(h1))
    imshow(filterImage);
    try
        disp('Select region to delete...');
        rect = getrect(h1);
    catch
        disp('caught');
    end

    rect = floor(rect);
    xmin = max(rect(2),1);
    xmax = min(rect(2) + rect(4),size(filterImage,1));
    ymin = max(rect(1),1);
    ymax = min(rect(1) + rect(3),size(filterImage,2));
    filterImage(xmin:xmax,ymin:ymax) = 0;
end

filter_new = filterImage;
h2 = figure; imshow(filter_new); title('Filter New');
img_output = filter_new;
[M1, M2] = ind2sub(size(filter_new), find(filter_new == 1));
xy_coordinates = [M1, M2];

