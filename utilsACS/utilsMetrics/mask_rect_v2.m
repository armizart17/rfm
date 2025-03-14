function [rect_mask] = mask_rect_v2(x, z, cx, cz, Lx, Lz)
% function [rect_mask] = mask_rect_v2(x, z, cx, cz, Lx, Lz)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION : Generic function for metric calculation, extract rectangular ROI (V2.0)
% rectangular region 
% INPUTS: 
%         - x  : Dim2 MATLAB (left-right) x_axis vector (only matters inital and final value)
%         - z  : Dim1 MATLAB (up-down) z_axis vector (only matters inital and final value) 
%         - cx : center horizontal position
%         - cz : center vertical position
%         - Lx : length horizontal
%         - Lz : length vertical
% Example Region of Interest (ROI)
%
%               ---------- Lx  ----------
%               |                       |
%              Lz       (cx, cz)        Lz
%               |                       |
%               ---------- Lx  ----------
%
% OUTPUTS: 
%         - rect_mask  : Mask (ROI are '1', others are '0')
% AUTHORs: Edmundo Arom Miranda 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [X_m, Z_m]= meshgrid(x, z);
    x0 = cx - Lx/2;
    z0 = cz - Lz/2;
    
    rect_mask = and( X_m >= x0, X_m <= x0+Lx ) & ...
                and( Z_m >= z0, Z_m <= z0+Lz );
    

end