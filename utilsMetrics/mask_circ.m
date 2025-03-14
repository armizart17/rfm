function [circ_mask, im_out] = mask_circ(x, z, cx, cz, radious, im_in)
% function [circ_mask, im_out] = mask_circ(x, z, cx, cz, radious, im_in)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION : circular mask based on center and radious
% 
% INPUTS: 
%         - x   : Dim2 MATLAB (left-right) x_axis vector (only matters inital and final value)
%         - z   : Dim1 MATLAB (up-down) z_axis vector (only matters inital and final value) 
%         - cx  : center position Dim2
%         - cz  : center position Dim1
%         - radious  : radious
%         - im_in: input image
% Example Region of Interest (ROI)
%
%        (---- radious ----  (cx, cz)  ---- radious -----)  
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
    [X_m,Z_m]= meshgrid(x_big,z_big);

    RR_m = sqrt( (X_m - cx).^2 + (Z_m - cz).^2 );
    circ_mask = RR_m <= radious;

    R_region = im_in .* circ_mask;
    R_region(R_region == 0) = NaN;
    im_out = R_region;

end