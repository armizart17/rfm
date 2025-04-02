function [y_denoised, cost, err, fid, reg, iter] = TNV_regularization(y_array, lambda, max_iter, tol)
%TNV_REGULARIZATION Applies Total Nuclear Variation (TNV) denoising
%   y_denoised = tnv_regularization(y_array, lambda, max_iter)
%
%   Inputs:
%       - y_array:     Matrix of noisy signals (num_points, num_cols)
%       - lambda:      Regularization parameter (higher = more denoising)
%       - max_iter:    Maximum number of iterations (default: 100)
%
%   Outputs:
%       - y_denoised:  Denoised signal matrix (same size as y_array)

    % if nargin < 3
    %     max_iter = 100;  % Default maximum iterations
    % end

    % Initialize the denoised signal matrix
    y_denoised = y_array;
    
    % Convergence threshold
    % tol = 1e-6;  

    err = zeros(max_iter, 1);
    fid = err;
    reg = err;
    cost = err;
    cost(1) = 100;
    
    % Iteratively apply Singular Value Thresholding (SVT)
    for iter = 1:max_iter
        % Perform SVD
        [U, S, V] = svd(y_denoised, 'econ');
        
        % Apply soft-thresholding to the singular values
        S_thresholded = max(S - lambda, 0);
        
        % Reconstruct the denoised matrix
        y_new = U * S_thresholded * V';

        % fid(iter) = 0.5*norm(y_new - y_denoised, 'fro');
        fid(iter) = 0.5*sum((y_new(:) - y_array(:)).^2);
        reg(iter) = lambda*TV_nuclear(y_denoised);
        cost(iter+1) = fid(iter) + reg(iter);

        % v1
        % e = norm(y_new - y_denoised, 'fro') / norm(y_denoised, 'fro');
        
        % v2
        e = abs(cost(iter+1) - cost(iter)) / cost(iter+1);

        err(iter) = e;
        % Check for convergence
        if e < tol
            break;
        end
        
        % figure(18), 
        % plot(y_new(:,1:10:50), 'LineWidth', 1);
        % title(['TNV v1 Iter: ', num2str(iter)]);
        % xlabel('x-axis');
        % ylabel('y-axis');
        % ylim([-5 25]);
        % grid on;
        % 
        % pause(0.01)
        
        % Update the denoised matrix
        y_denoised = y_new;
        
    end
    
    fprintf('TNV Regularization completed in %d iterations\n', iter);
end

function [u] = opD(v)
    % Computes finite differences along rows
    [H, C] = size(v);
    % Compute differences along the rows (like taking gradient along depth)
    u1 = [diff(v, 1, 1); zeros(1, C)];  % Difference along rows (H direction)
    % Combine results along 3rd dimension
    u = cat(3, u1);
    % u = cat(3, [diff(v, 1, 1); zeros(1, C)]);  % Difference along the rows
end

function val = TV_nuclear(v)
    % Computes finite differences along rows
    [H, C] = size(v);
    % Compute differences along the rows (like taking gradient along depth)
    u1 = [diff(v, 1, 1); zeros(1, C)];  % Difference along rows (H direction)
    % Combine results along 3rd dimension
    y = cat(3, u1);
    % u = cat(3, [diff(v, 1, 1); zeros(1, C)]);  % Difference along the rows

    s = sqrt(sum(y.^2, 3));

    val = sum(s(:));
end