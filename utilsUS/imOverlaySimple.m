function [hF, hB, hColor] = imOverlaySimple(B, SWS, climB, clim, alpha, x, z, xBm, zBm)
% IMOVERLAYSIMPLE overlays the SWS image transparently over the B image.
% This version does not use an ROI argument and applies transparency uniformly.
%
%   Inputs:
%       B      - Background image (grayscale)
%       SWS    - Overlay image
%       climB  - Display range for background image
%       clim   - Display range for overlay image
%       alpha  - Transparency level (0 = transparent, 1 = opaque)
%       x, z   - Coordinates for the axes SWS
%       xBm, zBm   - Coordinates for the axes SWS
%
%   Outputs:
%       hF     - Handle to the overlay image
%       hB     - Handle to the background image
%       hColor - Handle to the colorbar

    % Normalize background image to [0, 1] range using climB
    B = repmat(mat2gray(double(B), double(climB)), [1, 1, 3]);

    % Display the background image
    hB = imagesc(xBm, zBm, B);
    axis image on;
    colormap(gray);
    hold on;

    % Set color limits for the overlay image
    factor = 0.01;
    if isempty(clim) % Default to 5% range expansion
        clim = [(1-factor)*min(SWS(:)), (1+factor)*max(SWS(:))];
    end
    if isequal(clim, [0 0]) % Handle edge cases
        clim = [-1 1];
    end
    % Ensure clim is in ascending order
    if clim(1) > clim(2)
        clim = sort(clim); % Sort clim to ensure it is ascending
    end
    % Display the overlay image
    hF = imagesc(x, z, SWS, clim);
    
    % Apply uniform transparency
    alphadata = (1-alpha) * ones(size(SWS)); % Uniform transparency
    set(hF, 'AlphaData', alphadata);

    % Add a colorbar for the overlay image
    hColor = colorbar;
    colormap(gcf, "turbo");
    hold off;
    axis image;
end


% function [hF, hB, hColor] = imOverlayNoInterp(B, SWS, climB, clim, alpha, x, z, ROI)
% % IMOVERLAYNOINTERP overlays the SWS image transparently over the B image
% % without interpolation, assuming they are already the same size.
% %
% %   Inputs:
% %       B      - Background image (grayscale)
% %       SWS    - Overlay image
% %       climB  - Display range for background image
% %       clim   - Display range for overlay image
% %       alpha  - Transparency level (0 = transparent, 1 = opaque)
% %       x, z   - Coordinates for the axes
% %       ROI    - Region of interest (binary mask)
% %
% %   Outputs:
% %       hF     - Handle to the overlay image
% %       hB     - Handle to the background image
% %       hColor - Handle to the colorbar
% 
%     % Normalize background image to [0, 1] range using climB
%     B = repmat(mat2gray(double(B), double(climB)), [1, 1, 3]);
% 
%     % Display the background image
%     hB = imagesc(x, z, B);
%     axis image on;
%     colormap(gray);
%     hold on;
% 
%     % Set color limits for the overlay image
%     factor = 0.01;
%     if isempty(clim) % Default to 5% range expansion
%         clim = [(1-factor)*min(SWS(:)), (1+factor)*max(SWS(:))];
%     end
%     if isequal(clim, [0 0]) % Handle edge cases
%         clim = [-1 1];
%     end
% 
%     % Display the overlay image
%     hF = imagesc(x, z, SWS, clim);
% 
%     % Apply transparency using the ROI mask
%     alphadata = (1-alpha) .* ROI;
%     set(hF, 'AlphaData', alphadata);
% 
%     % Add a colorbar for the overlay image
%     hColor = colorbar;
%     colormap("turbo");
%     hold off;
%     axis image;
% end
