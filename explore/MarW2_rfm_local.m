% Define reference depth intervals (1-mm spacing)
zr_spacing = 1e-3;  % 1 mm in meters
zr_values = z_ACS(1):zr_spacing:z_ACS(end); % Select depths at 1-mm steps

% Find closest indices for the reference depths
zr_indices = arrayfun(@(zr) find(abs(z_ACS - zr) == min(abs(z_ACS - zr)), 1), zr_values);
num_zr = length(zr_indices);

% Preallocate normalized spectrum array
RSp_norm = zeros(m, n, p, num_zr); 

% Compute Normalized Spectrum using Multiple zr Values
for zr_idx = 1:num_zr
    zr = zr_indices(zr_idx);
    Sp_ref = squeeze(Sp(zr, :, :)); % Extract reference spectrum at zr

    % Normalize spectra for each depth
    for ii = 1:m
        for jj = 1:n
            RSp_norm(ii, jj, :, zr_idx) = Sp(ii, jj, :) ./ Sp_ref(jj, :);
        end
    end
end

% Handle division errors (NaN, Inf)
RSp_norm(isnan(RSp_norm) | isinf(RSp_norm)) = 1; % Replace NaN/Infs with 1

%%
% Initialize matrices for LSM
y_vec = [];
X_mat = [];

for r = 1:num_zr % Loop over reference depths
    for k = 1:m   % Loop over depth positions
        for i = 2:p % Loop over frequency bins (starting from 2 for f_{i-1})
            % Define y = -ln(RSnor)
            y = -log(RSp_norm(k, :, i, r)); % Use all lateral positions

            % Define X = 4 (fi - fi-1) (zk - zr)
            X = 4 * (band(i) - band(i-1)) * (z_ACS(k) - z_ACS(zr_indices(r)));

            % Store values for least squares regression
            y_vec = [y_vec; y(:)];
            X_mat = [X_mat; repmat(X, length(y), 1)]; % Expand for all lateral positions
        end
    end
end

% Solve for 'a' using least squares: a = (X'X)^{-1}X'y
a_hat = (X_mat' * X_mat) \ (X_mat' * y_vec);

fprintf('Estimated Attenuation Coefficient: %.6f\n', a_hat);
