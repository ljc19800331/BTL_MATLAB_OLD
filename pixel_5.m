target_out = [];
colorvec_target_out = [];
r = 2;
for i = 1:length(pixel_out)
i 
p1 = pixel_out(i,:) + [-r, -r];
p2 = pixel_out(i,:) + [r, r];

x1 = p1(1); y1 = p1(2);
x2 = p2(1); y2 = p2(2);

% 3D colorized point cloud
idx_use = MTI_tex_after(:,1) >= y1 & MTI_tex_after(:,1) <= y2 & MTI_tex_after(:,2) >= x1 & MTI_tex_after(:,2) <= x2;
tex_use = MTI_tex_after(idx_use,:);
vtx_use = MTI_vtx_after(idx_use,:);
color_use = colorvec_after(idx_use, :);
target_out(i,:) = mean(vtx_use);
colorvec_target_out_1(i,:) = color_use(1,:);
end