function [u_opt,fcosts] = rpl_wRegtv_n_prior(Y,X,W,mu,u_0,params)
% function [u_opt,fcosts] = rpl_wRegtv_n_prior(Y,X,W,mu,u_0,params)
    
    if nargin < 6
        tol   = 1e-8;
        kmax  = 100;
        eps_f = 1e-8;
        m_est = 0;
    else
        tol   = params.tol;
        kmax  = params.kmax;
        eps_f = params.eps_f;
        m_est = params.m_est;
        df_op = params.df_op;
        c     = params.c;
        sigma = params.sigma;
        weight= params.weight; 
    end

    % spectrum ratio
    y = Y(:);

    [r,p,q] = size(Y); % old
    % [p,q,r] = size(Y);
    
    % hyperparameters
    c_fid = c(1);
    c_b   = c(2);
    c_a   = c(4);

    sigma_fid = sigma(1);
    sigma_b   = sigma(2);
    sigma_a   = sigma(4);

    mu_b = mu(1); 
    mu_a = mu(3);

    % differential operators
    if df_op == 1
        dx = 0.5*(diag(ones(q-1,1),1) - diag(ones(q-1,1),-1));
        dx(1,1) = -1; dx(1,2) = 1; dx(end,end) = 1; dx(end,end-1) = -1;
        dy = 0.5*(diag(ones(p-1,1),1) - diag(ones(p-1,1),-1));
        dy(1,1) = -1; dy(1,2) = 1; dy(end,end) = 1; dy(end,end-1) = -1;
    else
        dx = diag(ones(q-1,1),1) - diag([ones(q-1,1);0]);
        dy = diag(ones(p-1,1),1) - diag([ones(p-1,1);0]);
    end

    % differential matrices for b and n
    Dx = sparse(kron(dx,speye(p)));
    Dy = sparse(kron(speye(q),dy));

    % differential matrices for a
    Rx = Dx*Dy;
    Ry = Dy*Dy;

    % differential matrix-based operators
    D = [Dx;Dy]; Dt = D';
    R = [Rx;Ry]; Rt = R';

    % precomputation of matrices
    Xt = X';
    Wt = W';

    omega_q_fid = speye(length(y));

    omega_r_b   = speye(2*p*q);
    omega_r_a   = speye(2*p*q);

    % preallocation of vector for reporting functional costs
    fid  = zeros(kmax,1);
    tv_b = zeros(kmax,1);
    tv_a = zeros(kmax,1);
    f    = zeros(kmax,1);
    err  = ones(kmax,1);

    % initializing iterations
    k = 1;

    % tolerance for f_R function
    eps_s = 1e-4;

    % initial solution
    b = u_0(1:p*q);
    a = u_0(p*q+1:end);

    % errors for each term in the optimization problem
    e_fid = X*b + W*a - y;

    % e_b   = norm2_fun(D*b,p*q); % €
    % e_a   = norm2_fun(R*a,p*q); % €

    e_b   = D*b;
    e_a   = R*a;

    % costs at k = 1
    fid(1)  = 0.5*(norm(e_fid))^2;
    tv_b(1) = 0.5*mu_b*(norm(e_b))^2;
    tv_a(1) = 0.5*mu_a*(norm(e_a))^2;
    f(1)    = fid(1)+tv_b(1)+tv_a(1);
    err(1)  = abs(f(1));

    % error
    e = err(1);

    s_fid_mad = sigma_fid*median(abs(e_fid - median(e_fid)));
    s_b_mad   = sigma_b  *median(abs(e_b   - median(e_b)));
    s_a_mad   = sigma_a  *median(abs(e_a   - median(e_a)));

    % if m_est == 1
    %     w_fid = (abs(e_fid) <= c_fid*s_fid_mad).*(1 - (e_fid/(c_fid*s_fid_mad)).^2).^2;
    %     w_b   = (abs(e_b)   <= c_b  *s_b_mad)  .*(1 - (e_b  /(c_b  *  s_b_mad)).^2).^2;
    %     w_a   = (abs(e_a)   <= c_a  *s_a_mad)  .*(1 - (e_a  /(c_a  *  s_a_mad)).^2).^2;
    % else
    %     w_fid = 1./(1 + (e_fid/(c_fid*s_fid_mad + eps)).^2);
    %     w_b   = 1./(1 + (e_b/(c_b*s_b_mad       + eps)).^2);
    %     w_a   = 1./(1 + (e_a/(c_a*s_a_mad       + eps)).^2);
    % end

    % weights in regularization
    w_flatten = weight(:);

    w_b = repmat(w_flatten, 2, 1);
    w_a = repmat(w_flatten, 2, 1);

    omega_r_b   = spdiags(w_b,0,omega_r_b);
    omega_r_a   = spdiags(w_a,0,omega_r_a);

    % No Weights in Fidelity

    o_flatten = ones(size(w_flatten));
    w_fid = repmat(o_flatten, r, 1); % expand for each channel
    
    omega_q_fid = spdiags(w_fid,0,omega_q_fid);
    
    

    while e > tol

        s_b = f_R(norm2_fun(D*b,p*q).^2,eps_s,eps_f);
        s_a = f_R(norm2_fun(R*a,p*q).^2,eps_s,eps_f);

        omega_s_b = spdiags(kron(ones(2,1),s_b),0,speye(2*p*q));
        omega_s_a = spdiags(kron(ones(2,1),s_a),0,speye(2*p*q));

        omega_q_b = omega_s_b*omega_r_b;   
        omega_q_a = omega_s_a*omega_r_a;

        Q_b = Xt*omega_q_fid*X + mu_b*Dt*omega_q_b*D;
        Q_a = Wt*omega_q_fid*W + mu_a*Rt*omega_q_a*R;

        y_b   = (Xt*omega_q_fid)*y - (Xt*omega_q_fid*W)*a;
        [b,~] = cgs(Q_b,y_b,tol,10000);

        y_a   = (Wt*omega_q_fid)*y - (Wt*omega_q_fid*X)*b;
        [a,~] = cgs(Q_a,y_a,tol,10000);

        e_fid = X*b + W*a - y;

        % e_b   = norm2_fun((omega_s_b.^0.5)*D*b,p*q); % €
        % e_a   = norm2_fun((omega_s_a.^0.5)*R*a,p*q); % €
        
        e_b   = (omega_s_b.^0.5)*D*b;
        e_a   = (omega_s_a.^0.5)*R*a;

        fid(k+1)  = 0.5*(norm(e_fid))^2;

        % tv_b(k+1) = mu_b*sum(e_b); % € 
        % tv_a(k+1) = mu_a*sum(e_a); % €
        
        tv_b(k+1) = 0.5*mu_b*(norm((omega_r_b.^0.5)*e_b))^2; 
        tv_a(k+1) = 0.5*mu_a*(norm((omega_r_a.^0.5)*e_a))^2;
       
        f(k+1)    = fid(k+1) + tv_b(k+1) + tv_a(k+1);
        err(k+1)  = abs((f(k+1)-f(k))/f(k+1));

        e = err(k+1);
        k = k+1;

        if kmax == k
            break;
        end

        % s_fid_mad = sigma_fid*median(abs(e_fid - median(e_fid)));
        % s_b_mad   = sigma_b  *median(abs(e_b   - median(e_b)));
        % s_a_mad   = sigma_a  *median(abs(e_a   - median(e_a)));

        % if m_est == 1
        %     w_fid = (abs(e_fid) <= c_fid*s_fid_mad).*(1 - (e_fid/(c_fid*s_fid_mad)).^2).^2;
        %     w_b   = (abs(e_b)   <= c_b  *s_b_mad)  .*(1 - (e_b  /(c_b  *  s_b_mad)).^2).^2;
        %     w_a   = (abs(e_a)   <= c_a  *s_a_mad)  .*(1 - (e_a  /(c_a  *  s_a_mad)).^2).^2;
        % else
        %     w_fid = 1./(1 + (e_fid/(c_fid*s_fid_mad + eps)).^2);
        %     w_b   = 1./(1 + (e_b/(c_b*s_b_mad       + eps)).^2);
        %     w_a   = 1./(1 + (e_a/(c_a*s_a_mad       + eps)).^2);
        % end

        % omega_q_fid = spdiags(w_fid,0,omega_q_fid);

        % omega_r_b   = spdiags(kron(ones(2,1),w_b),0,omega_r_b); % €
        % omega_r_a   = spdiags(kron(ones(2,1),w_a),0,omega_r_a); % €
        
        % omega_r_b   = spdiags(w_b,0,omega_r_b);
        % omega_r_a   = spdiags(w_a,0,omega_r_a);


    end

    u_opt = [b;a];

    % functional costs
    fcosts.fid   = fid(1:k);
    fcosts.reg_b = tv_b(1:k);
    fcosts.reg_a = tv_a(1:k);
    fcosts.f     = f(1:k);

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

% function u_opt = rpl_robust_tv(Y,X,Z,W,mu,c,sigma,params)
% 
%     if nargin < 6
%         tol   = 1e-8;
%         kmax  = 100;
%         eps_f = 1e-8;
%         m_est = 0;
%     else
%         tol   = params.tol;
%         kmax  = params.kmax;
%         eps_f = params.eps_f;
%         m_est = params.m_est;
%     end
% 
%     y = Y(:);
%     p = size(Y,2);
%     q = size(Y,3);
% 
%     c_fid = c(1);
%     c_b   = c(2);
%     c_n   = c(3);
%     c_a   = c(4);
% 
%     sigma_fid = sigma(1);
%     sigma_b   = sigma(2);
%     sigma_n   = sigma(3);
%     sigma_a   = sigma(4);
% 
%     mu_b = mu(1);
%     mu_n = mu(2);
%     mu_a = mu(3);
% 
%     dx = diag(ones(q-1,1),1) - diag([ones(q-1,1);0]);
%     dy = diag(ones(p-1,1),1) - diag([ones(p-1,1);0]);
%     Dx = sparamsse(kron(dx,speye(p)));
%     Dy = sparamsse(kron(speye(q),dy));
%     Rx = Dx*Dy;
%     Ry = Dy*Dy;
% 
%     D = [Dx;Dy];
%     R = [Rx;Ry];
% 
%     Dt = D'; Rt = R';
% 
%     Xt = X';
%     Zt = Z';
%     Wt = W';
% 
%     omega_q_fid = speye(length(y));
% 
%     omega_r_b   = speye(2*p*q);
%     omega_r_n   = speye(2*p*q);
%     omega_r_a   = speye(2*p*q);
% 
%     fid  = zeros(kmax,1);
%     tv_b = zeros(kmax,1);
%     tv_n = zeros(kmax,1);
%     tv_a = zeros(kmax,1);
%     f    = zeros(kmax,1);
%     err  = ones(kmax,1);
% 
%     Q_b = Xt*X + mu_b*Dt*D;
%     Q_n = Zt*Z + mu_n*Dt*D;
%     Q_a = Wt*W + mu_a*Rt*R;
% 
%     y_b   = Xt*y;
%     [b,~] = cgs(Q_b,y_b,tol,10000);
% 
%     y_n   = Zt*y - (Xt*Z)*b;
%     [n,~] = cgs(Q_n,y_n,tol,10000);
% 
%     y_a   = Wt*y - ((Xt*W)*b+(Zt*W)*n);
%     [a,~] = cgs(Q_a,y_a,tol,10000);
% 
%     e_fid = X*b+Z*n+W*a - y;
%     e_b   = D*b;
%     e_n   = D*n;
%     e_a   = R*a;
% 
%     fid(1)  = 0.5*(norm(e_fid))^2;
%     tv_b(1) = 0.5*mu_b*(norm(e_b))^2;
%     tv_n(1) = 0.5*mu_n*(norm(e_n))^2;
%     tv_a(1) = 0.5*mu_a*(norm(e_a))^2;
%     f(1)    = fid(1)+tv_b(1)+tv_n(1)+tv_a(1);
%     err(1)  = abs(f(1));
% 
%     e = err(1);
%     k = 1;
%     eps_s = 1e-4;
% 
%     s_fid_mad = sigma_fid*median(abs(e_fid - median(e_fid)));
%     s_b_mad   = sigma_b  *median(abs(e_b   - median(e_b)));
%     s_n_mad   = sigma_n  *median(abs(e_n   - median(e_n)));
%     s_a_mad   = sigma_a  *median(abs(e_a   - median(e_a)));
% 
%     if m_est == 1
%         w_fid = (abs(e_fid) <= c_fid*s_fid_mad).*(1 - (e_fid/(c_fid*s_fid_mad)).^2).^2;
%         w_b   = (abs(e_b)   <= c_b  *s_b_mad)  .*(1 - (e_b  /(c_b  *  s_b_mad)).^2).^2;
%         w_n   = (abs(e_n)   <= c_n  *s_n_mad)  .*(1 - (e_n  /(c_n  *  s_n_mad)).^2).^2;
%         w_a   = (abs(e_a)   <= c_a  *s_a_mad)  .*(1 - (e_a  /(c_a  *  s_a_mad)).^2).^2;
%     else
%         w_fid = 1./(1 + (e_fid/(c_fid*s_fid_mad + eps)).^2);
%         w_b   = 1./(1 + (e_b/(c_b*s_b_mad       + eps)).^2);
%         w_n   = 1./(1 + (e_n/(c_n*s_n_mad       + eps)).^2);
%         w_a   = 1./(1 + (e_a/(c_a*s_a_mad       + eps)).^2);
%     end
% 
%     omega_q_fid = spdiags(w_fid,0,omega_q_fid);
% 
%     omega_r_b   = spdiags(w_b,0,omega_r_b);
%     omega_r_n   = spdiags(w_n,0,omega_r_n);
%     omega_r_a   = spdiags(w_a,0,omega_r_a);
% 
%     while e > tol
% 
%         s_b = f_R(norm2_fun(D*b,p*q).^2,eps_s,eps_f);
%         s_n = f_R(norm2_fun(D*n,p*q).^2,eps_s,eps_f);
%         s_a = f_R(norm2_fun(R*a,p*q).^2,eps_s,eps_f);
% 
%         omega_s_b = spdiags(kron(ones(2,1),s_b),0,speye(2*p*q));
%         omega_s_n = spdiags(kron(ones(2,1),s_n),0,speye(2*p*q));
%         omega_s_a = spdiags(kron(ones(2,1),s_a),0,speye(2*p*q));
% 
%         omega_q_b = omega_s_b*omega_r_b;
%         omega_q_n = omega_s_n*omega_r_n;
%         omega_q_a = omega_s_a*omega_r_a;
% 
%         Q_b = Xt*omega_q_fid*X + mu_b*Dt*omega_q_b*D;
%         Q_n = Zt*omega_q_fid*Z + mu_n*Dt*omega_q_n*D;
%         Q_a = Wt*omega_q_fid*W + mu_a*Rt*omega_q_a*R;
% 
%         y_b   = (Xt*omega_q_fid)*y - ((Xt*omega_q_fid*Z)*n+(Xt*omega_q_fid*W)*a);
%         [b,~] = cgs(Q_b,y_b,tol,10000);
% 
%         y_n   = (Zt*omega_q_fid)*y - ((Zt*omega_q_fid*X)*b+(Zt*omega_q_fid*W)*a);
%         [n,~] = cgs(Q_n,y_n,tol,10000);
% 
%         y_a   = (Wt*omega_q_fid)*y - ((Wt*omega_q_fid*X)*b+(Wt*omega_q_fid*Z)*n);
%         [a,~] = cgs(Q_a,y_a,tol,10000);
% 
%         e_fid = X*b+Z*n+W*a - y;
%         e_b   = (omega_s_b.^0.5)*D*b;
%         e_n   = (omega_s_n.^0.5)*D*n;
%         e_a   = (omega_s_a.^0.5)*R*a;
% 
%         fid(k+1)  = 0.5*(norm(e_fid))^2;
%         tv_b(k+1) = 0.5*mu_b*(norm((omega_r_b.^0.5)*e_b))^2;
%         tv_n(k+1) = 0.5*mu_n*(norm((omega_r_n.^0.5)*e_n))^2;
%         tv_a(k+1) = 0.5*mu_a*(norm((omega_r_a.^0.5)*e_a))^2;
%         f(k+1)    = fid(k+1)+tv_b(k+1)+tv_n(k+1)+tv_a(k+1);
%         err(k+1)  = abs((f(k+1)-f(k))/f(k+1));
% 
%         e = err(k+1);
%         k = k+1;
% 
%         if kmax == k
%             break;
%         end
% 
%         s_fid_mad = sigma_fid*median(abs(e_fid - median(e_fid)));
%         s_b_mad   = sigma_b  *median(abs(e_b   - median(e_b)));
%         s_n_mad   = sigma_n  *median(abs(e_n   - median(e_n)));
%         s_a_mad   = sigma_a  *median(abs(e_a   - median(e_a)));
% 
%         if m_est == 1
%             w_fid = (abs(e_fid) <= c_fid*s_fid_mad).*(1 - (e_fid/(c_fid*s_fid_mad)).^2).^2;
%             w_b   = (abs(e_b)   <= c_b  *s_b_mad)  .*(1 - (e_b  /(c_b  *  s_b_mad)).^2).^2;
%             w_n   = (abs(e_n)   <= c_n  *s_n_mad)  .*(1 - (e_n  /(c_n  *  s_n_mad)).^2).^2;
%             w_a   = (abs(e_a)   <= c_a  *s_a_mad)  .*(1 - (e_a  /(c_a  *  s_a_mad)).^2).^2;
%         else
%             w_fid = 1./(1 + (e_fid/(c_fid*s_fid_mad + eps)).^2);
%             w_b   = 1./(1 + (e_b/(c_b*s_b_mad       + eps)).^2);
%             w_n   = 1./(1 + (e_n/(c_n*s_n_mad       + eps)).^2);
%             w_a   = 1./(1 + (e_a/(c_a*s_a_mad       + eps)).^2);
%         end
% 
%         omega_q_fid = spdiags(w_fid,0,omega_q_fid);
% 
%         omega_r_b   = spdiags(w_b,0,omega_r_b);
%         omega_r_n   = spdiags(w_n,0,omega_r_n);
%         omega_r_a   = spdiags(w_a,0,omega_r_a);
% 
%     end
% 
%     u_opt = [b;n;a];
% 
% end
% 
% function s = f_R(d,eps_s,eps_f)
%     s = zeros(length(d),1);
%     d_s = sign(d).*(abs(d)+eps_s);
%     s(d > eps_f) = 2*(d_s(d > eps_f)).^(-0.5);
% end
% 
% function d = norm2_fun(m,n)
%     l = length(m)/n;
%     t = reshape(m,n,l);
%     d = sum(t.^2,2).^0.5;
% end