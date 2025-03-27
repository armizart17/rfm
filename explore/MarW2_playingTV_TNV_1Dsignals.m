
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
legend('Location', 'best')

%% TNV Regularization specs

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
lambda   = 1;  % Regularization parameter
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
title('EMZ TNV Denoised linear functions');
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
