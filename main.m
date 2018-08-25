% MAIN
test = BTL_ICRA();
[R_final, t_final, P_MTI_tform, P_Stereo_tform] = test.Calibration_3D23D();
%%
X_laser = [316, 338, 360, 383, 316, 338, 360, 382, 316, 339, 360, 382, 315, 339, 360, 382]; 
Y_laser = [212, 212, 212, 212, 234, 234, 234, 234, 255, 255, 255, 255, 274, 274, 274, 274];
[data_use] = test.Calibration_Interpolation(X_laser, Y_laser);
%%
test.PixelFromImg();