%% RANSAC for calibration
load('data_ransac.mat');

%% Load the data 
A = P_STEREO_3d;    % moving
B = P_MTI_3d;       % fixed

%% Plot the data points
% figure(1); clf; pcshow(A, 'r'); hold on; pcshow(B, 'b');
number = length(A);                 % Total number of points
num = 50;                           % minimum number of points
bestInNum = 0;                      % Best fitting line with largest number of inliers
bestParameter1=0;                   % parameters for best fitting line
bestParameter2=0;                   % Pending
inlierRatio = 0.33;                 % inlierratio is 50% 
threshDist = 0.05;                  % threshold distance is 0.05

for iter = 1: 10000

iter
idx = randperm(number, num); 

% Apply the transformation
moving = A(idx,:); fixed = B(idx,:); 
[R,t] = rigid_transform_3D(moving, fixed); 

% Apply the transform to all other data
ransac_1; 

% inliers and outliers
dis = A2 - B; 
dis = sqrt(dis(:,1).^2 + dis(:,2).^2 + dis(:,3).^2);
idx_in = dis < 0.05; 
idx_out = dis >= 0.05;
A2_in = A2(idx_in,:); 
B_in = B(idx_in,:); 
A2_out = A2(idx_out,:);
B_out = B(idx_out,:);

% continue threshold
dis_new = A2_in - B_in;  
dis_new = sqrt(dis_new(:,1).^2 + dis_new(:,2).^2 + dis_new(:,3).^2);
err_new = sqrt(sum(dis_new(:).^2)./length(dis_new)); 
N_inlier = length(dis_new);
% condition
if err_new < 0.05 & N_inlier > inlierRatio*number
    count_best = iter;
    R_best = R;
    t_best = t
    dis_best = A2_in - B_in;  
    dis_best = sqrt(dis_best(:,1).^2 + dis_best(:,2).^2 + dis_best(:,3).^2);
    err_best = sqrt(sum(dis_best(:).^2)./length(dis_best)); 
    break;
end

end

%% final err
n_final = length(A);
A_final = (R_best * A') + repmat(t_best, 1, n_final);
A_final = A_final';
figure(1);clf; 
pcshow(A_final, 'r'); hold on; pcshow(B, 'b');

err_final = A_final - B;
err_final = err_final .* err_final;
err_final = sum(err_final(:));
rmse_final = sqrt(err_final/n)

