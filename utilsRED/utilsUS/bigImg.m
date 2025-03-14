function [im_out_big] = bigImg(im_in, ref_out)
%function [im_out_big] = bigImg(im_in, ref_out)

    [M, N] = size(im_in);  % INPUT IMAGE
    [P, Q] = size(ref_out); % REFERENCE OUT

    [X1, Y1] = meshgrid( 1:N, 1:M );
    [X2, Y2] = meshgrid( linspace(1,N,Q), linspace(1,M,P) );   

    im_out_big = interp2(X1, Y1, im_in, X2, Y2);
    
end