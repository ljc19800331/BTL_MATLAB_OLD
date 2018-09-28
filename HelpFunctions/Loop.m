function [Res] = Loop( flag_xrot, Y_rot, Z_rot, img_scan)

    Res = {}; 
    count = 0;  
    tic
    for count_y = 1: length(Y_rot) 
        flag_yrot = Y_rot(count_y);
    %  for flag_yrot = Range_mat(3):flag_iter:Range_mat(4)
    for count_z = 1: length(Z_rot)   
        flag_zrot = Z_rot(count_z);
    %  for flag_zrot = Range_mat(5):flag_iter:Range_mat(6)
        count = count + 1;

    % Load the brain images

        imgpath = strcat( 'img_6D/', num2str(flag_xrot), '_', ...
                          num2str(flag_yrot), '_', ...
                          num2str(flag_zrot), '.png');
        img_brain = imread(imgpath);

    % Map the two images
        Img_1 = img_brain; Img_2 = img_scan; 
        Img_1(isnan(Img_1)) = 0;
        c = normxcorr2( Img_2, Img_1 ); 
        [max_c, imax] = max(abs(c(:)));   
        [ypeak, xpeak] = ind2sub(size(c),imax(1));
        x_center = xpeak - size(Img_2,2)/2;
        y_center = ypeak - size(Img_2,1)/2;
        data = [x_center, y_center, max_c, flag_xrot, flag_yrot, flag_zrot];
        Res{count_y, count_z} = [data]; 
        
    end
    end
    toc;
end
