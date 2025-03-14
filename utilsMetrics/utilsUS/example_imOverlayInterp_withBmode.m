%% Simple Script example to overlay SMALL colorImg to bmodeFull
% Function does interpolation inside

% INTERP OVERLAY BMODE + COLORIMAGE EXAMPLE

units = 1E3; % for [mm]
% units = 1E2; % for [cm]

bmodeFull       = ;
colorImg        = ;
range_bmode     = [-60 0];
range_img       = [0 1];
transparency    = 0.7;
x_img           = *units;
z_img           = *units;
xFull           = *units;
zFull           = *units;

[X, Z] = meshgrid(xFull, zFull);
roi = and(X >= x_img(1), X <= x_img(end)) & ...
      and(Z >= z_img(1), Z <= z_img(end));

figure, 

[~,hB,hColor] = imOverlayInterp(bmodeFull, colorImg, range_bmode, range_img, ...
                                transparency, x_img, z_img, roi, xFull, zFull);

hold on;
contour(xFull, zFull, roi, 1,'r--', 'LineWidth', 2)
hold off;

xlabel('Lateral'), ylabel('Axial');
hColor.Label.String = 'dB\cdotcm^{-1}\cdotMHz^{-1}';
title('ColorImage+Bmode')

%% FUNCTION
% function [hF,hB,hColor] = imOverlayInterp(B, SWS, climB, clim, alpha, x, z, ROI, xBm, zBm)
% % function [hF,hB,hColor] = imOverlayInterp(B, SWS, climB, clim, alpha, x, z, ROI, xBm, zBm)
% % IMOVERLAY(B,F) displays the image SWS transparently over the image B.
% %   Inputs:
% %       B      - Background image (grayscale)
% %       SWS    - Overlay image
% %       climB  - Display range for background image
% %       clim   - Display range for overlay image
% %       alpha  - Transparency level (0 = transparent, 1 = opaque)
% %       x, z   - Coordinates for the axes SWS
% %       ROI:    Region of interest
% %       xBm, zBm   - Coordinates for the axes SWS
% %
% %   Outputs:
% %       hF     - Handle to the overlay image
% %       hB     - Handle to the background image
% %       hColor - Handle to the colorbar
% 
% 
%     [X, Z] = meshgrid(x, z);
%     [Xq, Zq] = meshgrid(xBm, zBm);
% 
%     imgInterp = interp2(X,Z,SWS,Xq,Zq);
%     emptyRegion = isnan(imgInterp);
%     newRoi = ~emptyRegion & ROI;
%     % [hF,hB,hColor] = imoverlay2(B,imgInterp,climB,clim,alpha,x,z,newRoi,xBm,zBm);
% 
% 
%     B = repmat(mat2gray(double(B),double(climB)),[1,1,3]);
% 
%     hB = imagesc(xBm, zBm, B);
%     axis image on;
%     % xlabel('\bf x [mm]')
%     % ylabel('\bf z [mm]')
%     colormap(gray)
% 
%     hColor = colorbar; colormap("turbo");
%     hold on;
%     factor = 0.01;
%     if isempty(clim) % set 5% more of min and max values
%         clim = [ (1-factor)*min(SWS(:)), (1+factor)*max(SWS(:))] ;
%     end
%     if isequal(clim, [0 0])
%         clim = [-1 1];
%     end
% 
%     % Ensure clim is in ascending order
%     if clim(1) > clim(2)
%         clim = sort(clim); % Sort clim to ensure it is ascending
%     end
% 
%     hF = imagesc(x,z,imgInterp,clim);
% 
%     % If images are different sizes, map the front image to back coordinates
%     set(hF,'XData',get(hB,'XData'),'YData',get(hB,'YData'))
% 
%     % Make the foreground image transparent
%     alphadata = (1-alpha).*(newRoi);
%     set(hF,'AlphaData',alphadata);
%     hold off
%     axis image
% 
% end

