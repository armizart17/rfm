function [u_opt, ite] = ADMM_3dof_WeightedTV(A1, A2, A3, b, mu1, mu2, mu3, m, n, tol, mask,W)

    p = length(mask)/(m*n);
    minimask = reshape(mask, [m n p]);
    minimask = minimask(:, :, 1);
    minimask = minimask(:);
    b = mask .* b;
    b(isnan(b)) = 0;
    A = [A1 A2 A3];
    AtA = A' * A;
    Atb = A' * b;
    [u, ~] = cgs(AtA, Atb);
    B = reshape(u(1:end/3), m, n);
    C = reshape(u(end/3+1:2*end/3), m, n);
    N = reshape(u(2*end/3+1:end), m, n);

    % B = reshape(u(1:m*n), m, n);
    % C = reshape(u(m*n+1:2*m*n), m, n);
    % N = reshape(u(2*m*n+1:end), m, n);

    B = B(:);
    C = C(:);
    N = N(:);
    D = 0;
    v = 0;

    Wdiag = spdiags(W(:),0,m*n,m*n);

    F(1) = 1/2 * (norm(b - A1*B - A2*C - A3*N))^2 + ...
           mu1 * TVcalc_isotropic(B, m, n, minimask) + ...
           mu2 * TVcalc_isotropic(C, m, n, minimask) + ...
           mu3 * TVcalc_isotropic(N, m, n, minimask);

    ite = 0;
    error = 1;
    Bprev = B; 
    Cprev = C; 
    Nprev = N;
    
    while abs(error) > tol && ite < 50
        ite = ite + 1;

        rho = 1;
        % Update B
        B = IRLS_TV_weighted(b - A2*C - A3*N - D - v, A1, mu1/rho, m, n, tol, mask, minimask);

        % Update C
        C = IRLS_TV_weighted(b - A1*B - A3*N - D - v, A2, mu2/rho, m, n, tol, mask, minimask);

        % Update N
        N = IRLS_TV_weighted(b - A1*B - A2*C - D - v, A3, mu3/rho, m, n, tol, mask, minimask);

        % Update D
        w = b - A1*B - A2*C - A3*N - v;
        D = (rho / (rho + 1)) * w;

        % Update v
        v = v + A1*B + A2*C + A3*N + D - b;

        % Calculate objective function
        F(ite+1, 1) = 1/2 * (norm(b - A1*B - A2*C - A3*N))^2 + ...
                      mu1 * TVcalc_isotropic(B, m, n, minimask) + ...
                      mu2 * TVcalc_isotropic(C, m, n, minimask) + ...
                      mu3 * TVcalc_isotropic(N, m, n, minimask);

        % Check convergence
        error = sqrt(norm(B - Bprev)^2 + norm(C - Cprev)^2 + norm(N - Nprev)^2);
        Bprev = B;
        Cprev = C;
        Nprev = N;
    end

    u_opt = [B; C; N];
    disp('Number of iterations: ')
    disp(ite)

end
