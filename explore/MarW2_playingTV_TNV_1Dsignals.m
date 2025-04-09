
n = linspace(0, 4*pi, 100);
y = 0.5*randn(100,1) + sin(n(:)); % Generate noisy signal

function x_denoised = tv_denoise(y, lambda, iter)
    % Implements Chambolle's algorithm for Total Variation Denoising (1D)
    % 
    % Inputs:
    %   y      - Noisy input signal (1D array)
    %   lambda - Regularization parameter (higher = smoother)
    %   iter   - Number of iterations (default: 50)
    %
    % Output:
    %   x_denoised - Denoised signal
    y = y(:);
    if nargin < 3
        iter = 50; % Default number of iterations
    end
    
    % Initialize variables
    x = y; 
    p = zeros(length(y), 1); % p should be one element smaller than x
    tau = 0.25;

    for k = 1:iter
        % Compute divergence (div_p) and update x
        div_p = [diff(p); 0];  
        x_new = y - lambda * div_p; 
        
        % Update p, making sure dimensions match
        div_x_nw = [diff(x_new); 0];
        p = p + tau * div_x_nw;
        p = p ./ max(1, abs(p)); 
        
        % Update x
        x = x_new;
    end
    
    x_denoised = x;
end

%%

lambda = 0.2;  % Regularization parameter
x_denoised = tv_denoise(y, lambda, 100);

figure;
plot(y, 'r'); hold on; grid on;
plot(x_denoised, 'b', 'LineWidth', 2);
legend('Noisy Signal', 'Denoised Signal');
title('Total Variation Denoising - 1D');


% TV Denoising
mu = 0.2;
tol = 1e-4;
[M, N] = size(y);

[y_opt] = IRLS_TV(y,speye(M*N),mu,M,N,tol,ones(size(M*N)),ones(M*N,1));

figure;
plot(y, 'r', 'DisplayName', 'Original'); hold on; grid on;
plot(x_denoised, 'b--', 'LineWidth', 2, 'DisplayName', 'TV v1');
plot(y_opt, 'k--', 'LineWidth', 2, 'DisplayName', 'TV IRLS');
legend('Location', 'best');
title('Total Variation Denoising - 1D');

%% TNV Denoising 1D signals


clc; clear; close all;

% Generate synthetic 1D parallel signals with different noise levels
M = 5;   % Number of signals
N = 200; % Signal length
x = linspace(0, 2*pi, N);
X_clean = repmat(sin(x), M, 1); % Parallel clean signals

% Define noise levels
noise_levels = 1*[0.05, 0.1, 0.2, 0.3, 0.5]; % Different noise strengths
X_noisy = X_clean; % Initialize

for i = 1:M
    X_noisy(i, :) = X_clean(i, :) + noise_levels(i) * randn(1, N);
end

% Apply TNV denoising
lambda = 1.2;   % Regularization parameter
max_iter = 200; % Maximum iterations
tol = 1e-4;     % Convergence tolerance

X_denoised = total_nuclear_variation_denoising(X_noisy, lambda, max_iter, tol);

% Plot results
figure;
subplot(3,1,1);
plot(x, X_clean', 'k', 'LineWidth', 1.5);
title('Original Signals');

subplot(3,1,2);
plot(x, X_noisy', 'r');
title('Noisy Signals with Different Noise Levels');
hold on;
plot(x, X_denoised', 'b', 'LineWidth', 1.5);

subplot(3,1,3);
plot(x, X_denoised', 'b', 'LineWidth', 1.5);
title('Denoised Signals with TNV');

% Function for TNV Denoising
function X_denoised = total_nuclear_variation_denoising(X_noisy, lambda, max_iter, tol)
    [M, N] = size(X_noisy);
    D = diff(eye(N));  % Finite difference operator
    Dt = D';  % Transpose of D

    X = X_noisy;  % Initialize with noisy input
    step_size = 0.05;  % Fixed step size

    for iter = 1:max_iter
        DX = X * Dt;  % Compute finite differences
        
        % Singular Value Thresholding (SVT)
        [U, S, V] = svd(DX, 'econ');
        S_thresh = max(S - lambda, 0);  % Soft thresholding
        DX_new = U * S_thresh * V';  % Reconstruct denoised differences

        % Update X
        X_new = X_noisy - step_size * (X * Dt - DX_new) * D;

        % Convergence check
        if norm(X_new - X, 'fro') / norm(X, 'fro') < tol
            break;
        end
        X = X_new;
    end

    X_denoised = X;
end


%% Create signals (Mar W4)

% Parameters
num_points      = 100;         % Number of points in each column (x-axis samples)
num_cols        = 50;          % Number of columns (different intercepts)
slope           = 2;                 % Desired slope of the linear function
intercept_range = [-5, 5];  % Range of intercepts (randomly chosen for each column)
noise_level     = 1.5;          % Noise standard deviation

x = linspace(0, 10, num_points)';  % x-axis (column vector)
% x_array = x*ones(1, num_points);
x_array = repmat(x, [1, num_cols]);

% Preallocate the array to store noisy columns
y_array = zeros(num_points, num_cols);

% Generate noisy columns
for i = 1:num_cols
    intercept       = rand * (intercept_range(2) - intercept_range(1)) + intercept_range(1);  % Random intercept
    % disp(intercept);
    y_clean         = slope * x + intercept;  % Generate clean line
    y_noisy         = y_clean + noise_level * randn(size(y_clean));  % Add noise
    y_array(:, i)   = y_noisy;  % Store in the array
end

figure;
plot(x, y_array, 'LineWidth', 1);
title('Generated Linear Functions with Noise');
xlabel('x-axis');
ylabel('y-axis');
grid on;

% CGS SOLUTION ALL
b = y_array(:);
A = x_array(:);
At = A';
slope_cgs = cgs(At*A, At*b);
fprintf('Slope CGS: %f (GT= %.2f)\n', slope_cgs, slope);

% CGS PER ITERATION
vect = @(x) x(:);
slope_cgs_col = zeros(1, num_cols);
for i = 1:num_cols
    b  = vect(y_array(:, i));
    A  = vect(x_array(:, i));
    At = A';
    slope_cgs_col(1, i) = cgs(At*A, At*b);
end
fprintf('Slope CGS mean: %f (GT= %.2f)\n', mean(slope_cgs_col), slope);
figure, 
stem(slope_cgs_col),
yline(mean(slope_cgs_col), 'k--', 'DisplayName', num2str(mean(slope_cgs_col)))
xlabel('Ncol'), ylabel('Slope')
title('CGS per iteration')
legend('Location', 'best')

A_full    = kron(speye(num_cols), x);
At_full   = A_full';

b_full    = y_array(:);

slope_full = cgs(At_full*A_full, At_full*b_full);
slope_full = slope_full';

figure,
stem(slope_full, 'k')
yline(mean(slope_full), 'k--', 'DisplayName', num2str(mean(slope_full)))
xlabel('Ncol'), ylabel('Slope')
title('CGS ensemble')
legend('Location', 'best')

figure, 
stem(slope_full - slope_cgs_col, 'r')
title('\Delta')

%%
test_col = [1; 5 ;7]
mat_test = kron(test_col, speye(3))
mat_test2 = kron(speye(3), test_col)

figure, imagesc(mat_test2)
%%



%% TNV Regularization specs

% Parameters
num_points      = 100;         % Number of points in each column (x-axis samples)
num_cols        = 50;          % Number of columns (different intercepts)
slope           = 2;           % Desired slope of the linear function
intercept_range = [-5, 5];     % Range of intercepts (randomly chosen for each column)
noise_level     = 1.5;         % Noise standard deviation

x = linspace(0, 10, num_points)';  % x-axis (column vector)
% x_array = x*ones(1, num_points);
x_array = repmat(x, [1, num_cols]);

% Preallocate the array to store noisy columns
y_clean = zeros(num_points, num_cols);
y_array = zeros(num_points, num_cols);

% Generate noisy columns
for i = 1:num_cols
    intercept       = rand * (intercept_range(2) - intercept_range(1)) + intercept_range(1);  % Random intercept
    % disp(intercept);
    y_clean(:, i)   = slope * x + intercept;  % Generate clean line
    y_noisy         = y_clean(:, i)  + noise_level * randn(size(y_clean(:, i)));  % Add noise
    y_array(:, i)   = y_noisy;  % Store in the array
end

% CGS SOLUTION ALL
b = y_array(:);
A = x_array(:);
At = A';
slope_cgs = cgs(At*A, At*b);
fprintf('Slope CGS: %f (GT= %.2f)\n', slope_cgs, slope);

fprintf('================================================\n')
%%
%%%%%%%%%%%%%%%% TNV  %%%%%%%%%%%%%%%%
lambda   = 2;  % Regularization parameter
max_iter = 20;
tol = 2e-3;
[y_denoised, cost1, err1, fid1, reg1, iter1] = TNV_regularization(y_array, lambda, max_iter, tol);

b = y_denoised(:);
A = x_array(:);
At = A';
slope_cgs = cgs(At*A, At*b);
fprintf('TNV Denoised Slope CGS: %f (GT= %.2f)\n', slope_cgs, slope);

%%%%%%%%%%%%% TNV PRIMAL DUAL EMZ  %%%%%%%%%%%%%
weightEstimators = ones(num_cols, 1);
lambda = 2;
tau = 0.01;
maxIter = 300;
tol = 0.5e-3;
stableIter = 5;

[y_denoised2, cost2, err2, fid2, reg2] = pdo_den_wtnv_1d(y_array, lambda, tau, maxIter, tol, stableIter, weightEstimators);

b = y_denoised2(:);
A = x_array(:);
At = A';
slope_cgs = cgs(At*A, At*b);
fprintf('EMZ TNV Denoised Slope CGS: %f (GT= %.2f)\n', slope_cgs, slope);

%%

y_clean = y_temp;
x =  x_temp(:, 1);
figure;
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1800, 600]); % [x, y, width, height]   

subplot(141)
plot(x, y_clean(:,1:10:50), 'LineWidth', 1);
title('Clean Linear Functions');
xlabel('x-axis');
ylabel('y-axis');
ylim([-5 25]);
grid on;

subplot(142)
plot(x, y_array(:,1:10:50), 'LineWidth', 1);
title('Generated Linear Functions with Noise');
xlabel('x-axis');
ylabel('y-axis');
ylim([-5 25]);
grid on;

subplot(143)
plot(x, y_denoised(:,1:10:50), 'LineWidth', 1);
title('TNV Denoised linear functions');
xlabel('x-axis');
ylabel('y-axis');
ylim([-5 25]);
grid on;

subplot(144)
plot(x, y_denoised2(:,1:10:50), 'LineWidth', 1);
title('TNVv2 Denoised linear functions');
xlabel('x-axis');
ylabel('y-axis');
ylim([-5 25]);
grid on;


%%

figure, 
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1800, 600]); % [x, y, width, height]   

subplot(141)
plot(cost1, 'DisplayName', 'Cost1')
hold on, grid on
plot(cost2, 'DisplayName', 'Cost2')
hold off;
xlabel('Iter')
title('Cost')
legend('Location', 'northeast')
xlim([0 250])

subplot(142)
plot(fid1, 'DisplayName', 'Fid1')
hold on, grid on
plot(fid2, 'DisplayName', 'Fid2')
hold off;
xlabel('Iter')
title('Fid')
legend('Location', 'northeast')
xlim([0 250])

subplot(143)
plot(reg1, 'DisplayName', 'Reg1')
hold on, grid on
plot(reg2, 'DisplayName', 'Reg2')
hold off;
xlabel('Iter')
title('Reg')
legend('Location', 'northeast')
xlim([0 250])

subplot(144)
plot(err1, 'DisplayName', 'Err1')
hold on, grid on
plot(err2, 'DisplayName', 'Err2')
yline(tol, 'k--')
hold off;
xlabel('Iter')
title('Error')
legend('Location', 'northeast')
xlim([0 250])
ylim([0 0.005])
%%

% function [a_local_tnv, y_denoised] = TNV_regularization(x_temp, y_temp, lambda, max_iter)
%     % Inputs:
%     % - x_temp: The repeated vector of X (size: 47 x 90)
%     % - y_temp: The matrix to be denoised (size: 47 x 90)
%     % - lambda: Regularization parameter (default: 1e-3)
%     % - max_iter: Maximum number of iterations (default: 50)
% 
%     if nargin < 3, lambda = 1e-3; end
%     if nargin < 4, max_iter = 50; end
% 
%     % Flatten inputs
%     x_vec = x_temp(:);       % Flatten X
%     y_vec = y_temp(:);        % Flatten Y
% 
%     % Initialize variables
%     y_denoised = y_temp;  % Start with the original measurements
%     tol = 1e-6;
% 
%     for iter = 1:max_iter
%         % Step 1: Solve Linear System with `cgs()`
%         XTX = (x_vec' * x_vec + lambda * eye(length(x_vec)));
%         XTy = (x_vec' * y_denoised(:));
% 
%         % Using Conjugate Gradient Solver
%         a_local_tnv = cgs(XTX, XTy, 1e-6, 100);  % CGS instead of direct inversion
% 
%         % Step 2: Compute Residuals and Reshape to Matrix Form
%         y_denoised_matrix = reshape(y_vec - x_vec * a_local_tnv, size(y_temp));
% 
%         % Step 3: Apply TNV Denoising to the Matrix
%         y_denoised_matrix = apply_TNV(y_denoised_matrix, lambda);
% 
%         % Step 4: Flatten the result back to a vector
%         y_denoised = y_denoised_matrix(:);
% 
%         % Step 5: Check Convergence
%         if norm(x_vec * a_local_tnv - y_vec) / norm(y_vec) < tol
%             break;
%         end
%     end
% 
%     fprintf('TNV Regularization completed in %d iterations\n', iter);
% end
% 
% 
% function y_tnv = apply_TNV(y_matrix, lambda)
%     % Apply TNV denoising using Singular Value Thresholding (SVT)
%     [U, S, V] = svd(y_matrix, 'econ');
% 
%     % Apply soft-thresholding to the singular values
%     S_thresholded = max(S - lambda, 0);
% 
%     % Reconstruct the matrix
%     y_tnv = U * S_thresholded * V';
% end


%% RFM TNV plot

x =  x_temp(:, 1);
freq_range = round(linspace(1, 90, 5));

figure;
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1800, 600]); % [x, y, width, height]   

tiledlayout(1,3)

nexttile
plot(x, y_temp(:,freq_range), 'LineWidth', 1);
title('FPR');
xlabel('x-axis');
ylabel('y-axis');
% ylim([-5 25]);
grid on;

nexttile
plot(x, y_denoised1(:,freq_range), 'LineWidth', 1);
title('TNV v1');
xlabel('x-axis');
ylabel('y-axis');
% ylim([-5 25]);
grid on;

nexttile
plot(x, y_denoised2(:,freq_range), 'LineWidth', 1);
title('TNV v2 ');
xlabel('x-axis');
ylabel('y-axis');
% ylim([-5 25]);
grid on;

%%

ii = fix(m/2);
jj = fix(n/2);

% ii = fix(m/4);
% ii = fix(m/2);
% ii = fix(3*m/4);

x_temp = x_temp_all{ii, jj};
y_temp = y_temp_all{ii, jj};

% ALL TOGETHER
X_vec = x_temp(:);
y_vec = y_temp(:);
aaa = -(X_vec' * X_vec) \ (X_vec' * y_vec)
aaa = abs( (X_vec' * X_vec) \ (X_vec' * y_vec) )

lambda = 1;
% TNV v1
maxIter1 = 100;
tol1 = 2e-3;

% TNV2 
weights2 = ones(length(band_ufr), 1);
tau2 = 0.01;
maxIter2 = 300;
tol2 = 0.5e-3;
stableIter2 = 10;

%%%%%%%%%%%%%%%%%%%%% TNV DENOISING JOINT 1 %%%%%%%%%%%%%%%%%%%%%

[y_denoised1, cost1, err1, fid1, reg1, iter1] = TNV_regularization(y_temp, ...
    lambda, maxIter1, tol1);
y_temp1 = y_denoised1;
%%%%%%%%%%%%%%%%%%%%% TNV DENOISING JOINT 1 %%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% TNV DENOISING JOINT 2 %%%%%%%%%%%%%%%%%%%%%

[y_denoised2, cost2, err2, fid2, reg2] = pdo_den_wtnv_1d(y_temp, ...
    lambda, tau2, maxIter2, tol2, stableIter2, weights2);

y_temp2 = y_denoised2;
%%%%%%%%%%%%%%%%%%%%% TNV DENOISING JOINT 2 %%%%%%%%%%%%%%%%%%%%%

% =====================  ORIGINAL

A_full    = kron(speye(p_ufr), x_temp(:, 1));
At_full   = A_full';
b_full    = y_temp(:);

slope_full = -cgs(At_full*A_full, At_full*b_full);
slope_full = (slope_full)';

figure,
stem(slope_full, 'k')
yline(mean(slope_full), 'k--', 'DisplayName', [num2str(mean(slope_full)), '\pm', num2str(std(slope_full))  ])
yline(median(slope_full), 'k--', 'DisplayName', ['Med ', num2str(median(slope_full)) ])

xlabel('Ncol'), ylabel('Slope')
title('ORIGINAL')
legend('Location', 'best')

% CONSIDERING INTERCEPT
A1_full    = kron( speye(p_ufr), x_temp(:, 1) );
A2_full    = kron( speye(p_ufr), ones(size( x_temp(:, 1))) );
A_full = horzcat(A1_full, A2_full);

At_full   = A_full';
b_full    = y_temp(:);

x_opt = -cgs(At_full*A_full, At_full*b_full);
slope_full = x_opt(1:end/2);
inter_full = x_opt(end/2+1:end);

figure,
subplot(121)
stem(band_ufr, slope_full, 'k')
yline(mean(slope_full), 'k--', 'DisplayName', [num2str(mean(slope_full)), '\pm', num2str(std(slope_full))  ])
yline(median(slope_full), 'b--', 'DisplayName', ['Med ', num2str(median(slope_full)) ])
xlabel('Ncol'), ylabel('Slope')
title('ORIGINAL SLOPE')
legend('Location', 'best')

subplot(122)
stem(band_ufr, inter_full, 'k')
yline(mean(inter_full), 'k--', 'DisplayName', [num2str(mean(inter_full)), '\pm', num2str(std(inter_full))  ])
xlabel('Ncol'), ylabel('Slope')
title('ORIGINAL INTERCEPT')
legend('Location', 'best')

slope_fit_fx = zeros(1, p_ufr);
inter_fit_fx = zeros(1, p_ufr);
for pp = 1:p_ufr
    [slope_fx, intercept_fx, ~, ~] = fit_linear(x_temp(:, pp), y_temp(:, pp), 2);
    slope_fit_fx(1,pp) =    -(slope_fx) ;
    inter_fit_fx(1,pp) =    intercept_fx ;
end

figure, 
subplot(121)
plot(band_ufr, slope_fit_fx, '.-')
yline(mean(slope_fit_fx), 'k--', 'DisplayName', [num2str(mean(slope_fit_fx)), '\pm', num2str(std(slope_fit_fx))  ])
yline(median(slope_fit_fx), 'k--', 'DisplayName', ['Med ', num2str(median(slope_fit_fx))  ])

title('Slope Fit')
xlabel('Freq [MHz]')
legend('Location', 'best')
grid on;

subplot(122)
plot(band_ufr, inter_fit_fx, '.-')
yline(mean(inter_fit_fx), 'k--', 'DisplayName', [num2str(median(inter_fit_fx)), '\pm', num2str(std(inter_fit_fx))  ])
title('Intercept Fit')
xlabel('Freq [MHz]')
legend('Location', 'best')
grid on;


 % if (ii==fix(m/2) && jj==fix(n/6)) || (ii==fix(m/2) && jj==fix(n/2)) || (ii==fix(m/2) && jj==fix(5*n/6)) % different depths
 %        % if (jj==fix(n/2) && ii==fix(m/6)) || (jj==fix(n/2) && ii==fix(m/2)) || (jj==fix(n/2) && ii==fix(5*m/6)) % different laterals
freq1 = freqL+0.5; freq2 = 0.5*(freqL+freqH); freq3 = freqH-0.5;
idx_f1 = find(band_ufr >= freq1, 1, 'first'); idx_f2 = find(band_ufr >= freq2, 1, 'first'); idx_f3 = find(band_ufr >= freq3, 1, 'first');


[slope_f1, intercept_f1, ~, ~] = fit_linear(x_temp(:, idx_f1), y_temp(:, idx_f1), 2); 
[slope_f2, intercept_f2, ~, ~] = fit_linear(x_temp(:, idx_f2), y_temp(:, idx_f2), 2); 
[slope_f3, intercept_f3, ~, ~] = fit_linear(x_temp(:, idx_f3), y_temp(:, idx_f3), 2); 

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
        % end

% ===================== TNV v1

A_full    = kron(speye(p_ufr), x_temp(:, 1));
At_full   = A_full';
b_full    = y_denoised1(:);

slope_full = -cgs(At_full*A_full, At_full*b_full);
slope_full = slope_full';

figure,
stem(band_ufr , slope_full, 'k')
yline(mean(slope_full), 'k--', 'DisplayName', [num2str(mean(slope_full)), '\pm', num2str(std(slope_full))  ])
xlabel('Ncol'), ylabel('Slope')
title('TNV v1')
legend('Location', 'best')



[slope_f1, intercept_f1, ~, ~] = fit_linear(x_temp(:, idx_f1), y_denoised1(:, idx_f1), 2); 
[slope_f2, intercept_f2, ~, ~] = fit_linear(x_temp(:, idx_f2), y_denoised1(:, idx_f2), 2); 
[slope_f3, intercept_f3, ~, ~] = fit_linear(x_temp(:, idx_f3), y_denoised1(:, idx_f3), 2); 

figure;
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1000, 600]); % [x, y, width, height]          

title_f1 = sprintf('Freq %.2fMHz (a=%.2f)', band_ufr(idx_f1), slope_f1);
plot(x_temp(:, idx_f1), y_denoised1(:, idx_f1), 'r.-', 'DisplayName', title_f1 ); 
hold on; 
title_f2 = sprintf('Freq %.2fMHz (a=%.2f)', band_ufr(idx_f2), slope_f2);
plot(x_temp(:, idx_f2), y_denoised1(:, idx_f2), 'b.-', 'DisplayName', title_f2 ); 
title_f3 = sprintf('Freq %.2fMHz (a=%.2f)', band_ufr(idx_f3), slope_f3);
plot(x_temp(:, idx_f3), y_denoised1(:, idx_f3), 'k.-', 'DisplayName', title_f3 ); 
hold off; grid on;
title(sprintf('TNVv1 Data at ii = %d, jj = %d', ii, jj));
xlabel('X_t [cm]'); ylabel('RS_{norm} [dB/MHz]');
ylim([-20 20]);
legend('Location','best')
set(gca, 'FontSize', 22)

% ===================== TNV v2

A_full    = kron(speye(p_ufr), x_temp(:, 1));
At_full   = A_full';
b_full    = y_denoised2(:);

slope_full = -cgs(At_full*A_full, At_full*b_full);
slope_full = slope_full';

figure,
stem(slope_full, 'k')
yline(mean(slope_full), 'k--', 'DisplayName', [num2str(mean(slope_full)), '\pm', num2str(std(slope_full))  ])
xlabel('Ncol'), ylabel('Slope')
title('TNV v2')
legend('Location', 'best')

[slope_f1, ~, ~, ~] = fit_linear(x_temp(:, idx_f1), y_denoised2(:, idx_f1), 2); 
[slope_f2, ~, ~, ~] = fit_linear(x_temp(:, idx_f2), y_denoised2(:, idx_f2), 2); 
[slope_f3, ~, ~, ~] = fit_linear(x_temp(:, idx_f3), y_denoised2(:, idx_f3), 2); 

figure;
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1000, 600]); % [x, y, width, height]          

title_f1 = sprintf('Freq %.2fMHz (a=%.2f)', band_ufr(idx_f1), slope_f1);
plot(x_temp(:, idx_f1), y_denoised2(:, idx_f1), 'r.-', 'DisplayName', title_f1 ); 
hold on; 
title_f2 = sprintf('Freq %.2fMHz (a=%.2f)', band_ufr(idx_f2), slope_f2);
plot(x_temp(:, idx_f2), y_denoised2(:, idx_f2), 'b.-', 'DisplayName', title_f2 ); 
title_f3 = sprintf('Freq %.2fMHz (a=%.2f)', band_ufr(idx_f3), slope_f3);
plot(x_temp(:, idx_f3), y_denoised2(:, idx_f3), 'k.-', 'DisplayName', title_f3 ); 
hold off; grid on;
title(sprintf('TNVv2 Data at ii = %d, jj = %d', ii, jj));
xlabel('X_t [cm]'); ylabel('RS_{norm} [dB/MHz]');
ylim([-20 20]);
legend('Location','best')
set(gca, 'FontSize', 22)

%% APR W1

%% USE THIS FROM APRIL AND NOW ON
%%%%%%%%%%%%%%%%%%%% FAST WAY %%%%%%%%%%%%%%%%%%%%

% UFR strategy
bw_ufr = [3 9];
freqL = bw_ufr(1); freqH = bw_ufr(2);
range = bandFull >= freqL & bandFull <= freqH;

band_ufr    = bandFull(range);
p_ufr       = length(band_ufr);

RSp_k_ufr   = RSp_k(:,:,range);
RSp_r_ufr   = RSp_r(:,:,range);


% % DENOISING TNV RSP
mu_range  = 10.^linspace(log10(0.1), log10(10), 15);
tau_range = [0.0005 0.001 0.005 0.010 0.025 0.05];

% mu_range  = 10.^linspace(log10(0.1), log10(7), 5);
% tau_range = [0.01 ];

maxIter     = 1000;
stableIter  = 20;
tol         = 1e-3; % tolerance error

% gamma = 1;
% vec_gamma = 1 + (10 - 1) * (1 - linspace(0, 1, p_ufr).^gamma);
% weigthChannels = vec_gamma;
weigthChannels = ones(1, p_ufr);

tic
for tt = 1:length(tau_range)

for uu = 1:length(mu_range)

mu  = mu_range(uu);
tau = tau_range(tt);


[RSp_k_ufr_opt, cost, error, fid, reg] = pdo_den_wtnv(RSp_k_ufr, mu, tau, maxIter, tol, stableIter, weigthChannels);

[RSp_r_ufr_opt, cost, error, fid, reg] = pdo_den_wtnv(RSp_r_ufr, mu, tau, maxIter, tol, stableIter, weigthChannels);


% RFM METHOD

% Convert depth values to cm
z_ACS_cm = z_ACS * 1e2;      % Convert from meters to cm
z_ACS_r_cm = z_ACS_r * 1e2;  % Convert reference depths to cm

% Delta MHz 
df_MHz = band_ufr(2) - band_ufr(1);

% Preallocate cell arrays for storing x_temp and y_temp
x_temp_all = cell(m, n);
y_temp_all = cell(m, n);

% tic;
for jj = 1:n  % Loop over lateral positions (x_j)
    for ii = 1:m  % Loop over depth positions (z_k)

        % Temporary storage for this location
        y_temp = nan(m_r, p_ufr);  % (Reference depth, Frequency) ** m_r
        x_temp = nan(m_r, p_ufr);  % (Reference depth, Frequency) ** m_r

        for r = 1:m_r  % Loop over reference depths
            % if (ii==1 && r==1)
            y_col = squeeze( ( (RSp_k_ufr_opt(ii, jj, :)) - (RSp_r_ufr_opt(r, jj, :)) ) /(4*df_MHz) *Np2dB ); % p_ufr x 1

            X = z_ACS_cm(ii) - z_ACS_r_cm(r);

            y_temp(r, :) = y_col; %**
            x_temp(r, :) = X; % **

        end
        x_temp_all{ii, jj} = x_temp;
        y_temp_all{ii, jj} = y_temp;

    end
end
% t = toc;
% fprintf('Loop way Elapsed time %.2f \n', t);
%%%%%%%%%%%%%%%%%%%% FAST WAY %%%%%%%%%%%%%%%%%%%%

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
    
data_ACS_TNV.a_rfm_mu(:,:,uu) = a_rfm;
end

fulldataACS_TNV{tt} =  data_ACS_TNV;

end

t = toc;
fprintf('Loop way Elapsed time %.2f \n', t);

%%
%%
% FIGURES 

%%%%%%%%%%%%%%%%%%%%%%%%% TNV %%%%%%%%%%%%%%%%%%%%%%%%

caxis_acs = [0 1.1];

font = 30;
for tt = 1 : length(tau_range)

    figure (tt+40)
    set(tt+40,'units','normalized','outerposition',[0 0 1 1]);
    set(gca,'FontSize',30);
    
    for uu = 1 : length(mu_range)

        subplot (2,3,uu)
        a_rfm = fulldataACS_TNV{tt}.a_rfm_mu(:,:,uu);

        [m_a, s_a, cv_a] = deal(calc2dStats{1}(a_rfm), calc2dStats{2}(a_rfm), calc2dStats{3}(a_rfm));

        imagesc(a_rfm, caxis_acs), grid on 
        set(gca,'FontSize',12);
        % title(['\mu = ',num2str(mu_range(uu))], 'FontSize', font-15);
        title(sprintf('TNV-RFM ($\\mu$ = %.2f, GT = %.2f)\n%.3f $\\pm$ %.3f, CV = %.2f\\%%', ...
               mu_range(uu), alpha_sam, m_a, s_a, cv_a), ...
        'Interpreter', 'latex');

        % xlabel('Lateral [cm]', 'FontSize', 14) 
        % ylabel('Axial [cm]', 'FontSize', 14)
        colormap jet, 
        h1 = colorbar;
        % ylabel(h1,'dB.cm^{-1}.MHz^{-1}','FontSize', 14);
        
    end
    sgtitle(['\bf TNV RFM \tau = ', ...
        num2str(tau_range(tt)) ], ...
       'FontSize', font ); 
end
