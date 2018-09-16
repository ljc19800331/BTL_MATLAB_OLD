%% VTX and TEX Preprocessing -- before
% Remove the 0 index
idx_remove = (MTI_vtx_before(:,1) == 0) & (MTI_vtx_before(:,2) == 0) & (MTI_vtx_before(:,2) == 0);
MTI_vtx_before = MTI_vtx_before(~idx_remove,:);
MTI_tex_before = MTI_tex_before(~idx_remove,:);
% Texture to pixel
MTI_tex_before(:,1) = int16(MTI_tex_before(:,1) * 640);
MTI_tex_before(:,2) = int16(MTI_tex_before(:,2) * 480);
% Remove the minus index
idx_remove_1 = MTI_tex_before(:,1) >= 640 | MTI_tex_before(:,1) <= 0 | MTI_tex_before(:,2) >= 480 | MTI_tex_before(:,2) <= 0;
MTI_vtx_before = MTI_vtx_before(~idx_remove_1,:);
MTI_tex_before = MTI_tex_before(~idx_remove_1,:);
% Color vector 
colorvec_before = [];
for i = 1:length(MTI_tex_before(:,1))
    idx_width = MTI_tex_before(i,1);
    idx_height = MTI_tex_before(i,2);
    rgb = reshape(img_before(idx_height, idx_width, :), [1,3]);
    colorvec_before(i,1) = rgb(1); colorvec_before(i,2) = rgb(2); colorvec_before(i,3) = rgb(3);
end
colorvec_before = uint8(colorvec_before);