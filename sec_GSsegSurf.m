function [data, Cell] = sec_GSsegSurf(pt_use)

% pt_use = pc_STEREO.Location;

figure(1);clf;

%pt_use = data_icp;
x_0 = pt_use;
y_0 = pt_use(1,:);
[n_0, d_0] = knnsearch(x_0, y_0, 'k', 10, 'distance', 'euclidean');

x = pt_use;         % The main data set
y = pt_use(n_0,:);  % The closest n_0 points
result = [n_0(:)];  % The current result of the cloest points index 
Cell = {};          % The Cell -- final results 
flag_group = 0;     
data_viz = [y];
count = 0;

for i = 1: 1000

    %count = count + 1

    [n, d] = knnsearch(x, y, 'k', 10, 'distance', 'euclidean');

    idx_out = find(d(:) > 0.1);    % index that needs to be removed out 
    
    % index to tbe kept for each iteration
    n_iter = n(:);
    n_iter(idx_out) = [];
    n = n_iter;

    n_idx = unique(n(:));

    % storage the results
    result = [result; n_idx];

    y = pt_use(n_idx, :, :);
    data_viz = [data_viz;y];

    % Determine the next x and y
    idx_remove = n_idx;
    pt_use(idx_remove,:) = [];

    x = pt_use;

    % viz the result
    figure(1);hold on; 
    xlabel('x');ylabel('y');zlabel('z');
    pcshow(y,'r');

    % stop
    if length(y) == 0 & length(x) == 0
        flag_group = flag_group + 1;
        Cell{flag_group} = data_viz;
        break;
    end

    if length(y) == 0 & length(x) ~= 0
        flag_group = flag_group + 1;
        Cell{flag_group} = data_viz;
        data_viz = [];
        y = pt_use(1,:);
        data_viz = [data_viz;y];
    end

end
%toc;
% Viz the result
L = [];
% figure(1);clf
for i = 1: length(Cell) 
    res_viz = Cell{i};  
%     figure(1);pcshow(res_viz,rand(1,3));
%     hold on;
    L = [L;length(res_viz)];
end

% remove the other point cloud
[L_max, idx_max] = max(L);   
if length(Cell) == 0
    data_use = x_0;
else
    data_use = Cell{idx_max};
end
data = data_use;
