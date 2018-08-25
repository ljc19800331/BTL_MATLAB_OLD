% Map Stereo MTI
function [R_final, t_final, P_MTI_tform, P_Stereo_tform] = STEREO_MTI_MAP(P_MTI_3d, P_STEREO_3d_Seg)

    %%  Change the coordinate of the MTI points
    R_subx = [-1, 0, 0; 
          0, 1, 0;
          0, 0, 1];
    P_MTI_3d = P_MTI_3d * R_subx;

    %% change rotation of the moving object
    fixed = P_MTI_3d;
    moving = P_STEREO_3d_Seg;
    R_move = rotx(180); 
    moving = moving * R_move;
    
    %% normalize the data
    fixed_minvec = [min(fixed(:,1)), min(fixed(:,2)), min(fixed(:,3))];
    moving_minvec = [min(moving(:,1)), min(moving(:,2)), min(moving(:,3))];
    fixed = fixed - fixed_minvec;
    moving = moving - moving_minvec; 
    t_normal_MTI = fixed_minvec;
    t_normal_Stereo = -moving_minvec;
    
    %% ICP
    [tform, movingReg] = pcregrigid(pointCloud(moving), pointCloud(fixed), 'Tolerance', [0.0001, 0.00009]);
    R_ICP = tform.T(1:3,1:3);
    t_ICP = tform.T(4,1:3);
    moving_tform = moving * R_ICP + t_ICP;
    
    %% Final transformation
    R_final = R_move * R_ICP;
    t_final = t_normal_Stereo * R_ICP + t_normal_MTI + t_ICP;

    %% Tform results
    % MTI: fixed -- Stereo: moving
    P_Stereo_tform = P_STEREO_3d_Seg * R_final + t_final; 
    P_MTI_tform = P_MTI_3d;
   
end




