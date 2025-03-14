function [B,C,ite] = AlterOptiAdmmAnisWeighted(A1,A2,b,mu1,mu2,m,n,tol,mask,W)
% function [B,C,ite] = AlterOptiAdmmAnisWeighted(A1,A2,b,mu1,mu2,m,n,tol,mask,W)
% Solver for SLD with Anisotropic TV regularization for ACS and weighted 
% TV regularization for BSC.
%
% Approach by British columbia,
% SWTV-ACE: Spatially Weighted Regularization Based Attenuation Coefficient
% Estimation Method for Hepatic Steatosis Detection

    p = length(mask)/(m*n);
    minimask = reshape(mask,[m n p]);
    minimask = minimask(:,:,1);
    minimask = minimask(:);
    b = mask.*b;
    b(isnan(b)) = 0;
    A = [A1 A2];
    AtA = A'*A;
    Atb = A'*b;
    [u,~] = cgs(AtA,Atb);
    B = reshape(u(1:end/2),m,n);
    C = reshape(u(end/2+1:end),m,n);
    B = B(:);
    C = C(:);
    D = 0;
    v = 0;
    
    F(1) = 1/2*(norm( b - A1*B - A2*C ))^2 + mu1*TVcalc_anisotropic(B,m,n,minimask) + ...
        mu2*TVcalc_anisotropic(C,m,n,W);
    change(1) = NaN;
    
    ite  = 0;
    error = 1;
    
    Bprev = B;
    while abs(error) > tol && ite < 100
        ite = ite + 1;
        
        rho = 1;
        % First part of ADMM algorithm: ATTENUATION
        B = IRLS_ANIS_TV(b-A2*C-D-v,A1,mu1/rho,m,n,tol,mask,minimask);
        
        % Second part of ADMM algorithm: BACKSCKATTER
        C = IRLS_ANIS_TV_weighted(b-A1*B-D-v,A2,mu2/rho,m,n,tol,mask,minimask,W);
    
        % Third part of ADMM algorithm: D
        % Least squares: 1/2*||D||_2^2 + rho/2*||D-w||_2^2
        w = b - A1*B - A2*C - v;
        D = (rho/(rho+1))*w;
        
        % Fourth part of ADMM algorithm: v
        v = v + A1*B + A2*C + D - b;
        F(ite+1,1) = 1/2*(norm( b - A1*B - A2*C ))^2 + mu1*TVcalc_anisotropic(B,m,n,minimask) + ...
            mu2*TVcalc_anisotropic(C,m,n,W);
    
        error = norm(B - Bprev)/norm(Bprev);
        Bprev = B;
        % change(ite+1) = error;
    end
    % figure,tiledlayout(2,1)
    % nexttile,
    % plot(F)
    % grid on
    % title('Loss function')
    % nexttile,
    % plot(change)
    % grid on
    % title('change in solution')
end