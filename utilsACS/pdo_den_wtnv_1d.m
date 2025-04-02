function [x, cost, error, fid, reg] = pdo_den_wtnv_1d(y, lambda, tau, maxIter, tol, stableIter, weightEstimators)
%PDO_DEN_WTNV_1D Applies Weighted TNV Denoising to 2D Array of Signals
%
%   Inputs:
%       - y:              2D array (num_points, num_channels)
%       - lambda:         Regularization parameter
%       - tau:            Step size parameter
%       - maxIter:        Maximum number of iterations
%       - tol:            Convergence tolerance
%       - stableIter:     Number of iterations before checking tolerance
%       - weightEstimators: Weights for each channel (vector of length num_channels)
%
%   Outputs:
%       - x:              Denoised signals (same size as y)
%       - cost:           Cost function value per iteration
%       - error:          Error per iteration
%       - fid:            Data fidelity term
%       - reg:            Regularization term (nuclear norm)

    rho = 1.99;                % Relaxation parameter, in [1,2]
    sigma = 1/tau/8;           % Proximal parameter
    [H, C] = size(y);           % H = num_points, C = num_channels
    
    global weights
    weights = weightEstimators; % Weights for each channel (same length as num_channels)
 
    x2 = y;                     % Initialization of the solution
    u2 = zeros([size(y), 2]);    % Initialization of the dual solution

    cy = sum(y(:).^2) / 2;
    primalcostlowerbound = 0;
    
    ee = 1;  
    iter = 1;

    error(1) = 1;
    cost(1) = 100;
    fid(1) = 1;
    reg(1) = 1;
    
    while (iter < maxIter) && (ee > tol)  
        % Primal update
        x = prox_tau_f(x2 - tau * opDadj(u2), y, tau);
        
        % Dual update
        u = prox_sigma_g_conj(u2 + sigma * opD(2 * x - x2), lambda);

        % Update variables
        x2 = x2 + rho * (x - x2);
        u2 = u2 + rho * (u - u2);

        % Cost Calculation
        fid(iter+1)  = 0.5 * sum((x(:) - y(:)).^2);
        reg(iter+1)  = lambda * TVnuclear_EMZ(x);
        cost(iter+1) = fid(iter+1) + reg(iter+1);
        
        % Check convergence
        ee = abs(cost(iter+1) - cost(iter)) / cost(iter+1);
        error(iter+1) = ee;
        
        if (iter < stableIter) % Avoid stopping early in the first iterations
            ee = 1;
        end

        if mod(iter,5)==0
            fprintf('%4d  | %f | %e\n',iter, cost(iter+1), error(iter+1));
        end
        
        % figure(19), 
        % plot(x(:,1:10:50), 'LineWidth', 1);
        % title(['TNV v2 Iter: ', num2str(iter)]);
        % xlabel('x-axis');
        % ylabel('y-axis');
        % ylim([-5 25]);
        % grid on;
        % pause(0.01)

        iter = iter + 1;
    end
end

function [u] = prox_tau_f(v, B, tau)
    u = (v + tau * B) / (1 + tau);
end

function [u] = prox_sigma_g_conj(v, lambda)
    u = v - prox_nuc_weight(v, lambda);
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

function [u] = opDadj(v)
    % Computes the adjoint (divergence) of the operator opD
    % Inputs:
    % - v: 3D array of size [H, C, 2] (two gradient components)
    % Outputs:
    % - u: 2D array of size [H, C] (reconstructed signal)
    
    [H, C, ~] = size(v);
    
    % Adjoint operation along the first direction (rows - H)
    u1 = -[v(1,:,1); diff(v(:,:,1),1,1)];
    
    % Adjoint operation along the second direction (columns - C)
    u2 = -[v(:,1,2), diff(v(:,:,2),1,2)];
    
    % Combine the results
    u = u1 + u2;
end

function u = TVnuclear_EMZ(v)
    Jacobian = opD(v);  % Compute the gradient of the signal
    u = NuclearNorm(Jacobian);
end

function val = NuclearNorm(y)
    global weights;
    eps = 0.00001;

    s = sqrt(sum(y.^2, 3));
    % val = sum(s(:) ./ (weights(:) + eps));  % Weighted nuclear norm
    val = sum(s(:));
end

function xx = prox_nuc_weight(y, lambda)
    global weights;
    eps = 0.00001;
    s = sqrt(sum(y.^2, 3));
    scaling_factor = max(s - lambda, 0) ./ (s + eps);  % Soft-thresholding
    xx = y .* scaling_factor;  % Apply scaling to each element
end
