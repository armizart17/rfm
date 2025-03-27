% Preallocate attenuation matrix and storage arrays
a_local_ufr2 = zeros(m, n);  
Np2dB = 8.686;  % Conversion factor from Nepers to dB

% Preallocate cell arrays for storing x_temp and y_temp
x_temp_all = cell(m, n);
y_temp_all = cell(m, n);

df_MHz = min(diff(band_ufr));

tic;
for jj = 1:n  % Loop over lateral positions (x_j)
    for ii = 1:m  % Loop over depth positions (z_k)
        
        % Temporary storage for this location
        y_temp = zeros(p_ufr - 1, m_r);  % (Frequency, Reference depth)
        x_temp = zeros(p_ufr - 1, m_r);  % (Frequency, Reference depth)
        
        for r = 1:m_r  % Loop over reference depths
            for i = 2:p_ufr  % Loop over frequency bins

                % WAY ORIGINAL
                % y = ( log(RSp_k_ufr(ii, jj, i)) - log(RSp_r_ufr(r, jj, i)) );
                % X = -4 * (band_ufr(i) - band_ufr(i-1)) * (z_ACS_cm(ii) - z_ACS_r_cm(r)) / Np2dB;

                % WAY WITH UNITS
                y = ( log(RSp_k_ufr(ii, jj, i)) - log(RSp_r_ufr(r, jj, i)) ) / (4*df_MHz)*Np2dB;
                X = z_ACS_cm(ii) - z_ACS_r_cm(r);
     
                % Store y and X values in the temporary matrices
                y_temp(i-1, r) = y;  % (Frequency, Reference depth)
                x_temp(i-1, r) = X;  % (Frequency, Reference depth)
            end
        end
        
        % Store x_temp and y_temp for later use
        x_temp_all{ii, jj} = x_temp;
        y_temp_all{ii, jj} = y_temp;
        
        % % (Optional) Directly solve for local attenuation if desired
        % if ~isempty(y_temp)
        %     X_vec = x_temp(:);
        %     y_vec = y_temp(:);
        %     a_local_ufr2(ii, jj) = (X_vec' * X_vec) \ (X_vec' * y_vec);
        % end
    end
end
t = toc;
fprintf('Elapsed time: %.2f seconds\n', t);

%%

%% NOW ESTIMATION 

a_rfm = zeros(m, n); % Preallocate local attenuation matrix (depth x lateral)

for jj = 1:n  % Loop over lateral positions (x_j)
    for ii = 1:m  % Loop over depth positions (z_k)

        x_temp = x_temp_all{ii, jj};
        y_temp = y_temp_all{ii, jj};

        X_vec = x_temp(:);
        y_vec = y_temp(:);

        if (ii==fix(m/2) && jj==fix(n/6)) || (ii==fix(m/2) && jj==fix(n/2)) || (ii==fix(m/2) && jj==fix(5*n/6)) % different depths
        % if (jj==fix(n/2) && ii==fix(m/6)) || (jj==fix(n/2) && ii==fix(m/2)) || (jj==fix(n/2) && ii==fix(5*m/6)) % different laterals
            freq1 = freqL+0.5; freq2 = 0.5*(freqL+freqH); freq3 = freqH-0.5;
            idx_f1 = find(band_ufr >= freq1, 1, 'first'); idx_f2 = find(band_ufr >= freq2, 1, 'first'); idx_f3 = find(band_ufr >= freq3, 1, 'first');
            

            [slope_f1, ~, ~, ~] = fit_linear(x_temp(:, idx_f1), y_temp(:, idx_f1), 2); 
            [slope_f2, ~, ~, ~] = fit_linear(x_temp(:, idx_f2), y_temp(:, idx_f2), 2); 
            [slope_f3, ~, ~, ~] = fit_linear(x_temp(:, idx_f3), y_temp(:, idx_f3), 2); 

            figure;
            set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1000, 600]); % [x, y, width, height]          

            title_f1 = sprintf('Freq %.2fMHz (a=%.2f)', band_ufr(idx_f1), slope_f1);
            plot(x_temp(:, idx_f1), y_temp(:, idx_f1), 'r.-', 'DisplayName', title_f1 ); 
            hold on; 
            title_f2 = sprintf('Freq %.2fMHz (a=%.2f)', band_ufr(idx_f2), slope_f2);
            plot(x_temp(:, idx_f2), y_temp(:, idx_f2), 'b.-', 'DisplayName', title_f2 ); 
            title_f3 = sprintf('Freq %.2fMHz (a=%.2f)', band_ufr(idx_f3), slope_f3);
            plot(x_temp(:, idx_f3), y_temp(:, idx_f3), 'k.-', 'DisplayName', title_f3 ); 
            hold off; grid on;
            title(sprintf('Data at ii = %d, jj = %d', ii, jj));
            xlabel('X_t [cm]'); ylabel('RS_{norm} [dB/MHz]');
            ylim([-20 20]);
            legend('Location','best')
            set(gca, 'FontSize', 22)

            figure, 
            set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1400, 600]); % [x, y, width, height] 
            plot(y_vec); grid on;
            xlabel('Freq Groupes'); ylabel('RS_{norm} [dB/MHz]');
            title(sprintf('Data at ii = %d, jj = %d', ii, jj));
            x_max = length(y_vec);  % Maximum x value
            xticks(0:360:x_max);     % Set ticks at every 90 units
            set(gca, 'FontSize', 22)
        end
        %%%%%%%%%%%%%%%%%%%%%% PLOT FPR %%%%%%%%%%%%%%%%%%%%%%


        a_rfm(ii, jj) = -(X_vec' * X_vec) \ (X_vec' * y_vec);

    end
end


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
title(sprintf('RFM Local AC (GT= %.2f)\n%.3f $\\pm$ %.3f,  \\%%CV = %.2f', ...
               alpha_sam, m_a, s_a, cv_a), ...
      'Interpreter', 'latex');
set(gca,'fontsize',fontSize)

%% EXAMPLE ONE CASE @ ii & jj

ii = floor(m/2); jj =  floor(n/2);
% Retrieve the stored data
x_temp_data = x_temp_all{ii, jj};
y_temp_data = y_temp_all{ii, jj};

% Preallocate the result vector
a_vector = zeros(size(y_temp_data, 1), 1);  
snr1 = zeros(size(y_temp_data, 1), 1);  
snr2 = zeros(size(y_temp_data, 1), 1);  

% Perform linear fit for each row
for i = 1:size(y_temp_data, 1)
    % Extract the current row
    X_row = x_temp_data(i, :);
    Y_row = y_temp_data(i, :);    

    %%%%%%%%%%%%%%%% TV Denoising %%%%%%%%%%%%%%%%
    mu = 5;
    tol = 1e-4;

    snr1(i) = mean(Y_row)/std(Y_row);
    [M, N] = size(Y_row(:));
    [y_opt] = IRLS_TV(Y_row(:),speye(M*N),mu,M,N,tol,ones(size(M*N)),ones(M*N,1));
    Y_row = y_opt';
    y_temp_data(i, :) = Y_row;
    snr2(i) = mean(Y_row)/std(Y_row);
    %%%%%%%%%%%%%%%% TV Denoising %%%%%%%%%%%%%%%%

    % Solve the linear regression problem
    a_vector(i) = -(X_row * X_row') \ (X_row * Y_row');
end

% a_vector = - arrayfun(@(i) (x_temp_data(i, :) * x_temp_data(i, :)') \ ...
%     (x_temp_data(i, :) * y_temp_data(i, :)'), 1:size(y_temp_data, 1))';

% Plotting data for a particular frequency
freq1 = 3.5; freq2 = 6; freq3 = 8.5;
idx_f1 = find(band >= freq1, 1, 'first');
idx_f2 = find(band >= freq2, 1, 'first');
idx_f3 = find(band >= freq3, 1, 'first');

figure;
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1000, 600]); % [x, y, width, height]

subplot(2,1,1)
title_f1 = sprintf('Freq %.2fMHz (a=%.2f)', band(idx_f1), a_vector(idx_f1));
plot(x_temp_data(idx_f1, :), y_temp_data(idx_f1, :), 'r.-', 'DisplayName', title_f1 ); 
hold on; grid on;
title_f2 = sprintf('Freq %.2fMHz (a=%.2f)', band(idx_f2), a_vector(idx_f2));
plot(x_temp_data(idx_f2, :), y_temp_data(idx_f2, :), 'b.-', 'DisplayName', title_f2 ); 
title_f3 = sprintf('Freq %.2fMHz (a=%.2f)', band(idx_f3), a_vector(idx_f3));
plot(x_temp_data(idx_f2, :), y_temp_data(idx_f3, :), 'k.-', 'DisplayName', title_f3 ); 
% axis('auto')

title(sprintf('Data at ii = %d, jj = %d', ii, jj));
xlabel('X_t [cm]'); ylabel('RS_{norm} [dB/MHz]');
legend('Location','best')

%
subplot(2,1,2)
plot(band_ufr, a_vector), grid on
yline(mean(a_vector), 'k--')

title(sprintf('a @ ii = %d, jj = %d : %.2f ± %.2f ', ii, jj, mean(a_vector), std(a_vector)));

%% ALL IN MATRICES (ERRORS)

Np2dB = 8.686;  % Conversion factor from Nepers to dB

% Initialize big matrices/vectors
X_big_mat = [];
y_big_vec = [];
index_map = zeros(m, n);  % To store the starting index of each (ii, jj) pair

total_entries = 0;  % To keep track of the size of X_big_mat and y_big_vec

tic;
for jj = 1:n  % Loop over lateral positions (x_j)
    for ii = 1:m  % Loop over depth positions (z_k)
        
        y_vec = []; % Initialize y vector for this location
        X_mat = []; % Initialize X matrix for this location
        
        for r = 1:m_r  % Loop over reference depths
            for i = 2:p_ufr  % Loop over frequency bins
                % Compute y = log(RSnorm) at this depth & lateral position
                y = log(RSp_k_ufr(ii, jj, i)) - log(RSp_r_ufr(r, jj, i));
                
                % Define X = -4 * (fi - fi-1) * (zk - zr)
                X = -4 * (band_ufr(i) - band_ufr(i-1)) * (z_ACS_cm(ii) - z_ACS_r_cm(r)) / Np2dB;

                % Store values for least squares regression
                y_vec = [y_vec; y(:)];
                X_mat = [X_mat; X];
            end
        end
        
        % Save the current index position for this (ii, jj) pair
        index_map(ii, jj) = total_entries + 1;  
        
        % Update total entries count
        total_entries = total_entries + length(y_vec);
        
        % Accumulate X_mat and y_vec into the big matrix and vector
        X_big_mat = [X_big_mat; X_mat];  % Appending rows
        y_big_vec = [y_big_vec; y_vec];  % Appending rows
    end
end
t = toc;
fprintf('Building X_big_mat and y_big_vec Elapsed time: %.2f seconds\n', t);

% Solving the large system using Conjugate Gradient Solver
a_local_ufr_vec = cgs(X_big_mat' * X_big_mat, X_big_mat' * y_big_vec, 1e-6, 100);

figure, 
imagesc(reshape(a_local_ufr_vec, [m, n])), colorbar, colormap("jet");


%% TEST  (CHECK NOTEBOOK)

% Parameters
num_blocks = 5;     % Number of row vectors (number of blocks)
block_size = 3;     % Size of each row vector (number of columns in each block)

% Generate random row vectors
row_vectors = rand(num_blocks, block_size);  % (num_blocks x block_size)

% Preallocate the sparse matrix
X_array_sparse = sparse(num_blocks, num_blocks * block_size);

% Fill in the sparse matrix with row vectors at appropriate positions
for i = 1:num_blocks
    col_start = (i - 1) * block_size + 1;  % Starting column index for each block
    col_end = i * block_size;              % Ending column index for each block
    X_array_sparse(i, col_start:col_end) = row_vectors(i, :);  % Assign row vector to matrix
end

% Display the result (as full matrix for visualization)
disp(full(X_array_sparse));

% Display the size of the result
disp(['Size of X_array_sparse: ', mat2str(size(X_array_sparse))]);


band = linspace(10,20,5);
f = band(:);

X_array_sparse = kron( speye(num_blocks*block_size), f' );

% Display the size of the result
disp(['Size of X_array_sparse: ', mat2str(size(X_array_sparse))]);


figure, imagesc(full(X_array_sparse));