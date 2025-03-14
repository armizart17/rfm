function [B,C,ite] = AlterOpti_ADMM(A1,A2,b,mu1,mu2,m,n,tol,mask)

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
    %figure(109); imagesc(8.686*reshape(B,m,n)); colormap pink; caxis([0 1.2])
    B = B(:);
    C = C(:);
    D = 0;
    v = 0;
    
    F(1) = 1/2*(norm( b - A1*B - A2*C ))^2 + mu1*TVcalc_isotropic(B,m,n,minimask) + mu2*TVcalc_isotropic(C,m,n,minimask);
    
    ite  = 0;
    error = 1;
    Bprev = B; Cprev = C;
    while abs(error) > tol && ite < 50
        ite = ite + 1;
        
        rho = 1;
        % First part of ADMM algorithm: B
        B = IRLS_TV(b-A2*C-D-v,A1,mu1/rho,m,n,tol,mask,minimask);
        
        % Second part of ADMM algorithm: C
        C = IRLS_TV(b-A1*B-D-v,A2,mu2/rho,m,n,tol,mask,minimask);
        
        % Third part of ADMM algorithm: D
        % Least squares: 1/2*||D||_2^2 + rho/2*||D-w||_2^2
        w = b - A1*B - A2*C - v;
        D = (rho/(rho+1))*w;
        
        % Fourth part of ADMM algorithm: v
        v = v + A1*B + A2*C + D - b;
        F(ite+1,1) = 1/2*(norm( b - A1*B - A2*C ))^2 + mu1*TVcalc_isotropic(B,m,n,minimask) + mu2*TVcalc_isotropic(C,m,n,minimask);
        error = sqrt(norm(B - Bprev).^2 + norm(C - Cprev).^2);
        Cprev = C; Bprev = B;
    end
    disp('Number of iterations: ')
    disp(ite)

end

% Total Variation: 0.5*||A*u(:)-b||_2^2 + lambda*TV(u)
function u = IRLS_TV(b,A,mu,M,N,tol,mask,minimask)

    AtA = A'*A;
    Atb = A'*b;
    
    %figure(109); imagesc(8.686*reshape(u,M,N)); colormap pink; caxis([0 1.2])
    
    D = spdiags([-ones(M,1) ones(M,1)], [0 1], M,M+1);
    D(:,end) = [];
    D(M,M) = 0;
    Dx = kron(speye(N),D);
    
    D = spdiags([-ones(N,1) ones(N,1)], [0 1], N,N+1);
    D(:,end) = [];
    D(N,N) = 0;
    Dy = kron(D,speye(M));
    
    D = [Dx' Dy']';
    
    ite_irls = 0;
    error = 1;
    
    %[u,~] = cgs(AtA+mu*(D')*D,Atb);
    [u,~] = cgs(AtA,Atb);
    G(1) = 1/2*(norm( (b - A*u) ))^2 + mu*TVcalc_isotropic(u,M,N,minimask);
    
    while error > tol && ite_irls < 200
        
        X = reshape(u,M,N);
        ite_irls = ite_irls + 1;
        Dh = diff(X,[],1);
        Dh = [Dh;zeros(1,N)];
        Dv = diff(X,[],2);
        Dv = [Dv zeros(M,1)];
        
        vksquare = Dh.^2 + Dv.^2;
        vksquare = vksquare(:);
        
        eps = 0.3;
        P = sqrt(vksquare + eps^2);
        P = 1./P;
        %Huber seems not to work;
        
        %P = sqrt(P.^2 + eps^2);
        P = P(:).*minimask;
        omega = speye(M*N);   % sparse diagonal of just ones
        omega = spdiags(P,0,omega);   % sparse diagonal of P values instead of ones.
        W = kron(speye(2),omega);
        
        [u,~] = cgs(AtA + mu*D'*W*D, Atb, 1e-6, 200);
        G(ite_irls+1,1) = 1/2*(norm( (b - A*u) ))^2 + mu*TVcalc_isotropic(u,M,N,minimask);
        error = abs(G(ite_irls+1) - G(ite_irls));
        %error = (G(ite_irls+1) - G(ite_irls)).^2/G(ite_irls).^2;
        
    end

%figure(909); plot(1:length(G),G);

end

% TV Andres Leonel Coila
function [TV] = TVcalc_isotropic(B,M,N,mask)

    mask(isnan(mask)) = 0;
    mask = mask(:);
    
    X = reshape(B,M,N);
    Dh = diff(X,[],1);
    Dh = [Dh;zeros(1,N)];
    Dv = diff(X,[],2);
    Dv = [Dv zeros(M,1)];
    
    P = Dh.^2 + Dv.^2;
    %eps = 0.1;
    %eps = 0.05;
    %P = sqrt(P+eps^2);
    P = sqrt(P);
    TV = norm(P(:).*mask,1);

end


% TV Andres Leonel Coila - ANISOTROPIC
function [TV] = TVcalc_anisotropic(B,M,N,mask)

    mask(isnan(mask)) = 0;
    mask = mask(:);
    
    X = reshape(B,M,N);
    Dh = diff(X,[],1);
    Dh = [Dh;zeros(1,N)];
    Dv = diff(X,[],2);
    Dv = [Dv zeros(M,1)];
    
    P = abs(Dh) + abs(Dv);
    %P = sqrt(P);
    TV = norm(P(:).*mask,1);

end

