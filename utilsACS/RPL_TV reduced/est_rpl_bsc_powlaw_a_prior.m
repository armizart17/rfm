close all; clc;
warning('off');

%%

addpath('./est_methods/');
addpath('./utils/');
config_filename = './configs/config-sm-kw-powlaw-sam1.xlsx';
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
delta_a = test_params.a_sam - test_params.a_ref;
a_ref = test_params.a_ref;

comp_ref    = comp_ref_a(delta_a,j_ref,band,depth,q);
comp_freq_a = comp_mod_freq_a(a_ref,j_sam,j_ref,band,depth,q);

SR_comp = SR.*comp_ref.*comp_freq_a;

%%

est_method = test_params.est_method;

% log-spectrum Ratio
Y = log(SR_comp);

% matrices for RPL-based algorithms
X = kron(speye(p*q),ones(r,1));
Z = kron(speye(p*q),log(band));

% implementation parameters
rpl_params.df_op = test_params.df_op;

% parameters for RPL estimation algorithms
mu     = hyp_params.mu;
alpha  = hyp_params.alpha;
c_rtv  = hyp_params.c_rtv;
c_rtgv = hyp_params.c_rtgv;
sigma  = hyp_params.sigma_rtv;

% initialization for RPL-based methods
u_0 = initialize_rpl_a_prior(Y,X,Z,mu,rpl_params);

% RPL estimation
switch(est_method)
    case 0 % do nothing
        u_opt = u_0;
    case 1 % RPL-TV
        tic;
        [u_opt,fcosts] = rpl_tv_a_prior(Y,X,Z,mu,u_0,rpl_params);
        t = toc;
    otherwise
end

%%

b = exp(u_opt(1:p*q));
n = u_opt(p*q+1:2*p*q);

figure; imagesc(reshape(b,p,q)); colorbar; colormap('jet');
set(gca,'FontSize',20);
xlabel('Lateral distance [cm]','fontsize',20);
ylabel('Axial distance [cm]','fontsize',20);
    
figure; imagesc(reshape(n,p,q)); colorbar; colormap('jet');
set(gca,'FontSize',20)
xlabel('Lateral distance [cm]','fontsize',20);
ylabel('Axial distance [cm]','fontsize',20);