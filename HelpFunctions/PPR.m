% This script is used for invisible points removal
% Projection hidden point removal

test = sec_pre();
stl_name = 'brain_slice.stl';
[pc_test, data_xyz] = test.stl2pc(stl_name);
flag_type = 'brain';
flag_Len = 'inch';
pc_test = test.pcLenChange(pc_test, flag_type, flag_Len);
pc_test = test.BrainInit(pc_test);
pcshow(pc_test); grid on; title('Test stl2pc section'); 

% project the points to a 2D Plane 
pc_xyz = pc_test.Location;  
X = pc_xyz(:,1);
Y = pc_xyz(:,2);

% Choose the target region
% x:0-1, y:0.5-1.5
idx = pc_xyz(:,1) > 0 & pc_xyz(:,1) < 1 & pc_xyz(:,2) > 0.5 & pc_xyz(:,2) < 1.5;
pc_tar = pc_xyz(idx,:);
ColorVec_tar = pc_tar(:,3);
% Find the basic point (closet to the plane)
[value_basic, idx_basic] = max(pc_tar(:,3));
figure(1);subplot(1,2,1); pcshow(pc_tar);

pc_tar(:,3) = 0;
subplot(1,2,2);
pcshow(pc_tar, ColorVec_tar); hold on; 
plot3(pc_tar(idx_basic,1), pc_tar(idx_basic,2), pc_tar(idx_basic,3), 'rx');

% Project the 3D point cloud to the plane
ColorVec = pc_xyz(:,3);
%C = zeros( length(X), 3);
pc_xyz(:,3) = 0; 
pcshow(pc_xyz, ColorVec);

% Show the example code
ColorVec_exp1 = pc_exp1(:,3);
pc_exp1(:,3) = 0;
pcshow(pc_exp1, ColorVec_exp1);

