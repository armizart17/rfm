function u_opt = rpl_tv_ini(Y,X,Z,W,mu,par)

    if nargin < 6
        tol   = 1e-8;
        kmax  = 100;
        eps_f = 1e-8;
    else
        tol   = par.tol;
        kmax  = par.kmax;
        eps_f = par.eps_f;
    end

    y = Y(:);
    [p,q] = size(Y,[2,3]); % old
    % [p,q,r] = size(Y);

    mu_b = mu(1);
    mu_n = mu(2);
    mu_a = mu(3);

    dx = diag(ones(q-1,1),1) - diag([ones(q-1,1);0]);
    dy = diag(ones(p-1,1),1) - diag([ones(p-1,1);0]);
    Dx = sparse(kron(dx,speye(p)));
    Dy = sparse(kron(speye(q),dy));
    Rx = Dx*Dy;
    Ry = Dy*Dy;

    D = [Dx;Dy];
    R = [Rx;Ry];

    Dt = D'; Rt = R';

    XtX = X'*X;
    ZtZ = Z'*Z;
    WtW = W'*W;

    XtZ = X'*Z;
    XtW = X'*W;
    ZtW = Z'*W;
    ZtX = XtZ';
    WtX = XtW';
    WtZ = ZtW';

    Xty = X'*y;
    Zty = Z'*y;
    Wty = W'*y;

    fid  = zeros(kmax,1);
    tv_b = zeros(kmax,1);
    tv_n = zeros(kmax,1);
    tv_a = zeros(kmax,1);
    f    = zeros(kmax,1);
    err  = ones(kmax,1);

    e = err(1);
    k = 1;
    eps_s = 1e-4;

%     Q_b = XtX + mu_b*Dt*D;
%     Q_n = ZtZ + mu_n*Dt*D;
%     Q_a = WtW + mu_a*Rt*R;
% 
%     y_b   = Xty;
%     [b,~] = cgs(Q_b,y_b,tol,10000);
% 
%     y_n   = Zty - XtZ*b;
%     [n,~] = cgs(Q_n,y_n,tol,10000);
% 
%     y_a   = Wty - (XtW*b+ZtW*n);
%     [a,~] = cgs(Q_a,y_a,tol,10000);

    M = sparse([X,Z,W]); MtM = M'*M;
    Mty = M'*y(:);
    u = cgs(MtM,Mty,1e-5,100);
    b = u(1:p*q);
    n = u(p*q+1:2*p*q);
    a = u(2*p*q+1:end);

    e_fid = X*b+Z*n+W*a - y;
    e_b   = norm2_fun(D*b,p*q);
    e_n   = norm2_fun(D*n,p*q);
    e_a   = norm2_fun(R*a,p*q);

    fid(1)  = 0.5*(norm(e_fid))^2;
    tv_b(1) = mu_b*sum(e_b);
    tv_n(1) = mu_n*sum(e_n);
    tv_a(1) = mu_a*sum(e_a);
    f(1)    = fid(1)+tv_b(1)+tv_n(1)+tv_a(1);
    err(1)  = abs(f(1));

    while e > tol

        s_b = f_R(norm2_fun(D*b,p*q).^2,eps_s,eps_f);
        s_n = f_R(norm2_fun(D*n,p*q).^2,eps_s,eps_f);
        s_a = f_R(norm2_fun(R*a,p*q).^2,eps_s,eps_f);

        omega_b = spdiags(kron(ones(2,1),s_b),0,speye(2*p*q));
        omega_n = spdiags(kron(ones(2,1),s_n),0,speye(2*p*q));
        omega_a = spdiags(kron(ones(2,1),s_a),0,speye(2*p*q));
        
        Q_b = XtX + mu_b*Dt*omega_b*D;
        Q_n = ZtZ + mu_n*Dt*omega_n*D;
        Q_a = WtW + mu_a*Rt*omega_a*R;

        y_b   = Xty - (XtZ*n+XtW*a);
        [b,~] = cgs(Q_b,y_b,tol,10000);

        y_n   = Zty - (ZtX*b+WtZ*a);
        [n,~] = cgs(Q_n,y_n,tol,10000);

        y_a   = Wty - (WtX*b+WtZ*n);
        [a,~] = cgs(Q_a,y_a,tol,10000);

        e_fid = X*b+Z*n+W*a - y;
        e_b   = (omega_b.^0.5)*D*b;
        e_n   = (omega_n.^0.5)*D*n;
        e_a   = (omega_a.^0.5)*R*a;

        fid(k+1)  = 0.5*(norm(e_fid))^2;
        tv_b(k+1) = 0.5*mu_b*(norm(e_b))^2;
        tv_n(k+1) = 0.5*mu_n*(norm(e_n))^2;
        tv_a(k+1) = 0.5*mu_a*(norm(e_a))^2;
        f(k+1)    = fid(k+1)+tv_b(k+1)+tv_n(k+1)+tv_a(k+1);
        err(k+1)  = abs((f(k+1)-f(k))/f(k+1));

        e = err(k+1);
        k = k+1;

        if kmax == k
            break;
        end

    end

    u_opt = [b;n;a];
%     figure;
%     plot(10*log10(fid)); xlabel("iter"); ylabel("fid"); ylim([30 40]);
end

function s = f_R(d,eps_s,eps_f)
    s = zeros(length(d),1);
    d_s = sign(d).*(abs(d)+eps_s);
    s(d > eps_f) = 2*(d_s(d > eps_f)).^(-0.5);
end

function d = norm2_fun(m,n)
    l = length(m)/n;
    t = reshape(m,n,l);
    d = sum(t.^2,2).^0.5;
end