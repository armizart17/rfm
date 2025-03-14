function [rect_mask, im_out] = mask_rect(x, z, cx1, cx2, cz1, cz2, im_in)
% function [rect_mask, im_out] = mask_rect(x, z, cx1, cx2, cz1, cz2, im_in)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION : rectangular mask based on 
% rectangular region
% INPUTS: 
%         - x   : Dim2 MATLAB (left-right) x_axis vector (only matters inital and final value)
%         - z   : Dim1 MATLAB (up-down) z_axis vector (only matters inital and final value) 
%         - cx1 : left position  
%         - cx2 : right position
%         - cz1 : up position 
%         - cz2 : down position
%         - im_in: input image
% Example Region of Interest (ROI)
%
%           (cx1, cz1)  ---------  (cx2, cz1)  
%               |                       |
%               |                       |
%           (cx1, cz2)  ---------  (cx2, cz2)  
%
% OUTPUTS: 
%         - rect_mask  : Mask (ROI are '1', others are '0')
%         - im_out     : Output image (Only ROI values, others are NaN) 
% AUTHORs: Edmundo Arom Miranda 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [M, N] = size(im_in);
    x_big = linspace(x(1), x(end), N);
    z_big = linspace(z(1), z(end), M);
    [X_m,Z_m]= meshgrid(x_big, z_big);
    
    rect_mask = and( X_m >= cx1, X_m <= cx2 ) & ...
                and( Z_m >= cz1, Z_m <= cz2 );
    
    R_region = im_in .* rect_mask;
    R_region(R_region == 0) = NaN;
    im_out = R_region;

end

