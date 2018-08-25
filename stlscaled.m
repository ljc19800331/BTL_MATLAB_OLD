function [] = stlscaled()

test = sec_pre();

stl_name = 'brain_tumor.stl';
data_stl = stlread(stl_name);
V = data_stl.vertices; 
F = data_stl.faces;

% change the size for each dimension
x_range = max(V(:,1)) - min(V(:,1));
y_range = max(V(:,2)) - min(V(:,2));
z_range = max(V(:,3)) - min(V(:,3));

% Scale for each parameter
x_len = test.model_x;
y_len = test.model_y;
z_len = test.model_z;
c_x = x_len./x_range;         % cm/N
c_y = y_len./y_range;
c_z = z_len./z_range;

V(:,1) = V(:,1).* c_x;
V(:,2) = V(:,2).* c_y;
V(:,3) = V(:,3).* c_z;

% Rotate the model
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

V = V * R_x(1:3,1:3) * R_y(1:3,1:3) * R_z(1:3,1:3);
   
V(:,1) = V(:,1) - min(V(:,1));
V(:,2) = V(:,2) - min(V(:,2));
V(:,3) = V(:,3) - min(V(:,3));

pcshow(V);

writename = input('input the saved file name '); 
stlwrite(writename, F, V);

