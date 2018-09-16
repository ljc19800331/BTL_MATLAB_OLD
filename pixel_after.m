%% VTX and TEX Preprocessing -- after
% Remove the 0 index
idx_remove = (MTI_vtx_after(:,1) == 0) & (MTI_vtx_after(:,2) == 0) & (MTI_vtx_after(:,2) == 0);
MTI_vtx_after = MTI_vtx_after(~idx_remove,:);
MTI_tex_after = MTI_tex_after(~idx_remove,:);
% Texture to pixel
MTI_tex_after(:,1) = int16(MTI_tex_after(:,1) * 640);
MTI_tex_after(:,2) = int16(MTI_tex_after(:,2) * 480);
% Remove the minus index
idx_remove_1 = MTI_tex_after(:,1) >= 640 | MTI_tex_after(:,1) <= 0 | MTI_tex_after(:,2) >= 480 | MTI_tex_after(:,2) <= 0;
MTI_vtx_after = MTI_vtx_after(~idx_remove_1,:);
MTI_tex_after = MTI_tex_after(~idx_remove_1,:);
% Color vector 
colorvec_after = [];
for i = 1:length(MTI_tex_after(:,1))
    idx_width = MTI_tex_after(i,1);
    idx_height = MTI_tex_after(i,2);
    rgb = reshape(img_after(idx_height, idx_width, :), [1,3]);
    colorvec_after(i,1) = rgb(1); colorvec_after(i,2) = rgb(2); colorvec_after(i,3) = rgb(3);
end
colorvec_after = uint8(colorvec_after);