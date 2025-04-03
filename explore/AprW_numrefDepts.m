numRefsDepths = 45;  % Number of reference depths to use

tic;
for jj = 1:n  % Loop over lateral positions (x_j)
    for ii = 1:m  % Loop over depth positions (z_k)

        % Temporary storage for this location, preallocated with numRefsDepths rows
        y_temp = nan(numRefsDepths, p_ufr);
        x_temp = nan(numRefsDepths, p_ufr);
        
        % Determine the neighborhood of reference depths based on ii and numRefsDepths
        half_window = floor(numRefsDepths/2);
        if ii - half_window < 1
            % If near the top boundary, take the first numRefsDepths indices (or as many as available)
            r_window = 1:min(numRefsDepths, m_r);
        elseif ii + half_window > m_r
            % If near the bottom boundary, take the last numRefsDepths indices
            r_window = max(1, m_r - numRefsDepths + 1):m_r;
        else
            % Otherwise, center the window around ii
            if mod(numRefsDepths, 2) == 0
                % Even number: use half_window indices above and half_window below (adjusted)
                r_window = (ii - half_window):(ii + half_window - 1);
            else
                % Odd number: symmetric around ii
                r_window = (ii - half_window):(ii + half_window);
            end
        end
        
        % Loop over the selected reference depths in r_window
        for rr = 1:length(r_window)
            r = r_window(rr);  % Actual reference depth index

            % Compute the column difference for the current depth and reference depth
            y_col = squeeze(( log(RSp_k_ufr(ii, jj, :)) - log(RSp_r_ufr(r, jj, :)) ) / (4*df_MHz) * Np2dB);
            X = z_ACS_cm(ii) - z_ACS_r_cm(r);
            
            y_temp(rr, :) = y_col;
            x_temp(rr, :) = X;
        end
        
        % Store the temporary results
        x_temp_all{ii, jj} = x_temp;
        y_temp_all{ii, jj} = y_temp;
        
    end
end
t = toc;
fprintf('Loop with neighborhood ref Elapsed time %.2f \n', t);


% NOW ESTIMATION CGS RFM

a_rfm = zeros(m, n); 
for jj = 1:n  % Loop over lateral positions (x_j)
    for ii = 1:m  % Loop over depth positions (z_k)

        x_temp = x_temp_all{ii, jj};
        y_temp = y_temp_all{ii, jj};

        X_vec = x_temp(:);

        y_vec = y_temp(:);
        a_rfm(ii, jj) = -(X_vec' * X_vec) \ (X_vec' * y_vec);
    end
end
%
% RFM
[m_a, s_a, cv_a] = deal(calc2dStats{1}(a_rfm), calc2dStats{2}(a_rfm), calc2dStats{3}(a_rfm));
caxis_acs = [0 1.1];
fontSize = 14;

figure, 
imagesc(x_ACS * 1e3, z_ACS * 1e3, a_rfm, caxis_acs); % Convert to mm
axis("image")
colorbar; colormap("turbo")
xlabel('Lateral [mm]');
ylabel('Depth [mm]');
hb2=colorbar; ylabel(hb2,'dB\cdotcm^{-1}\cdotMHz^{-1}')
% title('Local Attenuation Coefficient');
title(sprintf('NumRefs %d | RFM AC (GT= %.2f)\n%.3f $\\pm$ %.3f,  \\%%CV = %.2f', ...
               numRefsDepths,alpha_sam, m_a, s_a, cv_a), ...
      'Interpreter', 'latex');
set(gca,'fontsize',fontSize)