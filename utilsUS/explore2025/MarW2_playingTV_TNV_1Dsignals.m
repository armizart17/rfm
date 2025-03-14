
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
M = 1;   % Number of signals
N = 200; % Signal length
x = linspace(0, 2*pi, N);
X_clean = repmat(sin(x), M, 1); % Parallel clean signals

% Define noise levels
noise_levels = 10*[0.05, 0.1, 0.2, 0.3, 0.5]; % Different noise strengths
X_noisy = X_clean; % Initialize

for i = 1:M
    X_noisy(i, :) = X_clean(i, :) + noise_levels(i) * randn(1, N);
end

% Apply TNV denoising
lambda = 1.2   % Regularization parameter
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
