% PNP recheck problem
data = load('worldToImageCorrespondences.mat');
cameraParams1 = data.cameraParams;

cameraParams2 = cameraParameters('IntrinsicMatrix', ...
    [929.769, 0, 0;0, 929.769, 0; 628.21, 357.027, 1]);

%% 
imagePoints = readNPY('P_MTI_Calibration_interpolation_2D_before.npy');
imagePoints = double(imagePoints);

worldPoints = readNPY('P_MTI_Calibration_interpolation_3D_before.npy');

% change from inch to mm
worldPoints = worldPoints(1:length(imagePoints),:);
worldPoints = worldPoints % .* 25.4; % inch to mm

[worldOrientation, worldLocation] = estimateWorldCameraPose(imagePoints, worldPoints, cameraParams2);

%% Reproject image points to 3D points
reproj_imagePoints = worldToImage(cameraParams2, worldOrientation, worldLocation, worldPoints);

figure(1);clf
imshow('template.png');
hold on;
plot(imagePoints(:,1), imagePoints(:,2), 'rx');
hold on;
plot(reproj_imagePoints(:,1), reproj_imagePoints(:,2), 'bx');

%% Check the results
figure(1);clf 
pcshow(worldPoints,'VerticalAxis','Y','VerticalAxisDir','down', ...
     'MarkerSize',30);
hold on
plotCamera('Size',10,'Orientation',worldOrientation,'Location',...
     worldLocation);
hold off

%% 
[worldOrientation,worldLocation] = estimateWorldCameraPose(...
     data.imagePoints,data.worldPoints,data.cameraParams);

pcshow(data.worldPoints,'VerticalAxis','Y','VerticalAxisDir','down', ...
     'MarkerSize',30);
hold on
plotCamera('Size',10,'Orientation',worldOrientation,'Location',...
     worldLocation);
hold off

% Camera paramter 

% image points

% World coordinates

% Check 3D to 2D 

% Reproject 3D to 2D to check the results
