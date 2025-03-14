function u_0 = initialize_rpl_b_prior(Y,Z,W,mu,params)

    if nargin < 5
        ini_method = 0;
    else
        ini_method = params.ini_method;
        ini_tol = params.ini_tol;
    end

    y = Y(:);

    [p,q] = size(Y,[2,3]); % old
    % [p,q,r] = size(Y);
    
    mu_n = mu(2); mu_a = mu(3);

    dx = diag(ones(q-1,1),1) - diag([ones(q-1,1);0]);
    dy = diag(ones(p-1,1),1) - diag([ones(p-1,1);0]);
    Dx = sparse(kron(dx,speye(p)));
    Dy = sparse(kron(speye(q),dy));
    Rx = Dx*Dy;
    Ry = Dy*Dy;

    D = [Dx;Dy];
    R = [Rx;Ry];

    M = [Z,W]; MtM = M'*M; Mty = M'*y;

    switch ini_method
        case 0
            u_0 = sparse([],[],[],2*p*q,1);
        case 1
            [u_0,~] = cgs(MtM,Mty,ini_tol,10000);
        case 2
            L = [sqrt(mu_n)*D,sqrt(mu_a)*R]; LtL = L'*L;
            [u_0,~] = cgs(MtM+LtL,Mty,ini_tol,10000);
        case 3
            [u_0,~] = cgs(MtM,Mty,ini_tol,10000);
            [u_0,~] = rpl_tv_b_prior(Y,X,W,mu,u_0,params);
        otherwise
            u_0 = sparse([],[],[],2*p*q,1);
    end

end