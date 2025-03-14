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

delta_n = 0;

comp_ref = comp_ref_n_bsc(delta_n,band,p,q);

SR_comp = SR.*comp_ref;

%%

est_method = test_params.est_method;

% indices initialization
f = band;
[r,p,q] = size(SR_comp);

% log-spectrum Ratio
Y = log(SR_comp);

% matrices for RPL-based algorithms
X = kron(speye(p*q),ones(r,1));
W = -4*kron(speye(p*q),f');

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

b = u_opt(1:p*q);
a = u_opt(p*q+1:2*p*q);

z = 1e+2*repmat(depth,1,q);
dz = reshape(Dy*z(:),p,q);
dz(end,:) = dz(end-1,:);

% please, add the reference values

figure; imagesc(reshape(exp(b),p,q)); colorbar; colormap('jet');
set(gca,'FontSize',20);
xlabel('Lateral distance [cm]','fontsize',20);
ylabel('Axial distance [cm]','fontsize',20);
    
figure; imagesc(reshape(8.686*Dy*a./dz(:),p,q)); colorbar; colormap('jet');
set(gca,'FontSize',20)
xlabel('Lateral distance [cm]','fontsize',20);
ylabel('Axial distance [cm]','fontsize',20);