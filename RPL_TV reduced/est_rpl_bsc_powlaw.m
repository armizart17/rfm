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

j_sam = 1;
j_ref = 1;
a_ref = 0;

comp_freq_a = comp_mod_freq_a(a_ref,j_sam,j_ref,band,depth,q);

SR_comp = SR.*comp_freq_a;

%%

est_method = test_params.est_method;

% log-spectrum Ratio
Y = log(SR_comp);

% matrices for RPL-based algorithms
X = kron(speye(p*q),ones(r,1));
Z = kron(speye(p*q),log(band));
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
u_0 = initialize_rpl(Y,X,Z,W,mu,rpl_params);

% RPL estimation
switch(est_method)
    case 0 % do nothing
        u_opt = u_0;
    case 1 % RPL-TV
        tic;
        [u_opt,fcosts] = rpl_tv(Y,X,Z,W,mu,u_0,rpl_params);
        t = toc;
    case 2 % RPL-TGV
        tic;
        [u_opt,fcosts] = rpl_tgv(Y,X,Z,W,alpha,u_0,rpl_params);
        t = toc;
    case 3 % Robust RPL-TV
        tic;
        u_opt = rpl_robust_tv(Y,X,Z,W,mu,c_rtv,sigma,u_0,rpl_params);
        t = toc;
    case 4 % Robust RPL-TGV
        tic;
        sigma_rtgv = 1.4484*ones(7,1);
        u_opt = rpl_robust_tgv(Y,X,Z,W,alpha,c_rtgv,sigma_rtgv,u_0,rpl_params);
        t = toc;
    case 5 % RPL-VTV
        tic;
        mu = 2.5e+5;
        [u_opt,fcosts] = rpl_robust_wvtv(Y,X,Z,W,mu,c_rtv,sigma,u_0,rpl_params);
        t = toc;
    otherwise
end

%%

ref_vals  = gen_refvals(band,data,test_params);
est_maps  = gen_estmaps(u_opt,p,q,depth,band,ref_vals,test_params);
plot_maps = set_plotmaps(est_maps,test_params);

%plot_maps{3}.param_roi = delta_SNR;

figure(1); QUS_parameter_display(data,plot_maps{1},test_params,0);

for k = 1:length(plot_maps)
   figure(k+1); QUS_parameter_display(data,plot_maps{k},test_params,1);
end

%%

src        = test_params.src;
stats_flag = src < 3;

if stats_flag == 1
    sam_maps = data.sam_maps;
else
    sam_maps = cell(3,1);
end

stats = cell(3+2*(src == 3),1);

for k = 1:length(sam_maps)
    stats{k} = get_metrics(est_maps{k},sam_maps{k},stats_flag,spectral_data,test_params);
end

if src > 2
    stats{k+1} = get_metrics(est_maps{k+1},0,0,spectral_data,test_params);
    stats{k+2} = get_metrics(est_maps{k+2},0,0,spectral_data,test_params);
end

%%

frec = (3:0.05:8)';

b = est_maps{1};
n = est_maps{2};

bsc_fit_powlaw = mean(b(:))*frec.^mean(n(:));

figure;
plot(frec,10*log10(bsc_fit_powlaw),'r');

% figure; %(k+2);
% subplot(221); plot(fcosts.fid,'r');
% subplot(222); plot(fcosts.reg_b,'g');
% subplot(223); plot(fcosts.reg_n,'b');
% subplot(224); plot(fcosts.reg_a,'m');