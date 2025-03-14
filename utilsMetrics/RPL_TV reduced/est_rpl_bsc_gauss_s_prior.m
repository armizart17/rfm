close all; clc;
warning('off');

%%

addpath('./est_methods/');
addpath('./utils/');

[test_params,rpl_params,hyp_params] = set_params(config_filename);

% reading US data
data = read_USdata(test_params);

% compute spectral data
spectral_data = compute_spectral_data(data,test_params);

band  = spectral_data.band;
depth = spectral_data.depth;
SR    = spectral_data.SR;

[r,p,q] = size(SR);

j_sam = test_params.j_sam;
j_ref = test_params.j_ref;
delta_s = test_params.n_sam - test_params.n_ref;
a_ref = test_params.a_ref;

comp_ref    = comp_ref_s_bsc(delta_s,band,p,q);
comp_freq_a = comp_mod_freq_a(a_ref,j_sam,j_ref,band,depth,q);

SR_comp = SR.*comp_ref.*comp_freq_a;

%%

est_method = test_params.est_method;

% log-spectrum Ratio
Y = log(SR_comp);

% matrices for RPL-based algorithms
X = kron(speye(p*q),ones(r,1));
W = -4*kron(speye(p*q),band.^j_sam);

% implementation parameters
rpl_params.df_op = test_params.df_op;

% parameters for RPL estimation algorithms
mu     = hyp_params.mu;
alpha  = hyp_params.alpha;
c_rtv  = hyp_params.c_rtv;
c_rtgv = hyp_params.c_rtgv;
sigma  = hyp_params.sigma_rtv;

% initialization for RPL-based methods
u_0 = initialize_rpl_n_prior(Y,X,W,mu,rpl_params);

% RPL estimation
switch(est_method)
    case 0 % do nothing
        u_opt = u_0;
    case 1 % RPL-TV
        tic;
        [u_opt,fcosts] = rpl_tv_n_prior(Y,X,W,mu,u_0,rpl_params);
        t = toc;
    otherwise
end

%%

dy = 0.5*(diag(ones(p-1,1),1) - diag(ones(p-1,1),-1));
dy(1,1) = -1; dy(1,2) = 1; dy(end,end) = 1; dy(end,end-1) = -1;
Dy = sparse(kron(speye(q),dy));

k = u_opt(1:p*q);
a = u_opt(p*q+1:2*p*q);

z = 1e+2*repmat(depth,1,q);
dz = reshape(Dy*z(:),p,q);
dz(end,:) = dz(end-1,:);

frec = (3:0.05:8)';

k_est = mean(exp(k));
s_est = -0.0245;

bsc_est_gauss = k_est.*exp(s_est.*frec.^2);

M = [ones(length(frec),1),log(frec)];
r = log(bsc_est_gauss);

x = cgs(M'*M,M'*r,1e-16,100);

b = exp(x(1));
n = x(2);

bsc_fit_powlaw = b*frec.^n;

k_teo = 1;
s_teo = -0.0245;

bsc_teo = k_teo.*exp(s_teo.*frec.^2);

e_r = norm(bsc_est_gauss - bsc_teo,2)/norm(bsc_teo,2);

figure;
plot(frec,10*log10(bsc_est_gauss),'r'); hold on;
plot(frec,10*log10(bsc_fit_powlaw),'g'); hold on;
plot(frec,10*log10(bsc_teo),'b'); grid on;
legend('gaussian-est','powlaw-fit','gaussian-teo');
set(gca,'FontSize',20);

figure; imagesc(reshape(8.686*Dy*a./dz(:),p,q)); colorbar; colormap('jet');
set(gca,'FontSize',20); %16
xlabel('Lateral distance [cm]','fontsize',20);
ylabel('Axial distance [cm]','fontsize',20);
set(gca,'FontSize',20);