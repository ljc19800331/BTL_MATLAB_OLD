target_label = [];
colorvec_target_label = [];
r = 2;
for i = 1:length(pixel_label)
i 
p1 = pixel_label(i,:) + [-r, -r];
p2 = pixel_label(i,:) + [r, r];

x1 = p1(1); y1 = p1(2);
x2 = p2(1); y2 = p2(2);

% 3D colorized point cloud
idx_use = MTI_tex_before(:,1) >= y1 & MTI_tex_before(:,1) <= y2 & MTI_tex_before(:,2) >= x1 & MTI_tex_before(:,2) <= x2;
tex_use = MTI_tex_before(idx_use,:);
vtx_use = MTI_vtx_before(idx_use,:);
color_use = colorvec_before(idx_use, :);
target_label(i,:) = mean(vtx_use);
colorvec_target_label_1(i,:) = color_use(1,:);
end