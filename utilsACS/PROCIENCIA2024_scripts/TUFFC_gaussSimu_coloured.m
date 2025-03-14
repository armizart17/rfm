% Test Gaussian BSC coloured by line
% Dec 2024 for TUFFC
% By EAMZ

% clear all; 
% clc;
warning('off');

Np2dB = 20*log10(exp(1));
dB2Np = 1/Np2dB;
range_bmode = [-60 0];
% Bmode = @(RF) 20*log10(abs(hilbert(RF))) - max (20*log10(abs(hilbert(RF))));
plotBSCdB = true;

% estim_method = '3DoF';
% estim_method = '2DoF-a';
% estim_method = '2DoF-g'; 
estim_method = '2DoF-s'; 


%% PATH DATA
pathData = '.\NewPL_data\data_ok\';

%%

numCase = 3;
switch numCase
    case 2  % VARIATION a, g, s
        a_sam_the = 0.6; a_ref_the = 0.4;
        g_sam_the = 1; g_ref_the = 0.01;
        s_sam_the = -0.20; s_ref_the = 0;
    case 3  % VARIATION a, g, s
        a_sam_the = 0.6; a_ref_the = 0.4;
        g_sam_the = 0.01; g_ref_the = 1;
        s_sam_the = 0.0; s_ref_the = -0.1;
end

% AUX = load([pathData,'ref0.mat']); % JUST AUXILIAR TO GIVE AXIS (GAUSSIAN SIMU HAVE THIS, NO NECESSARY)
name_sam= sprintf('ref_a%.1fg%.2fs%.4f', a_sam_the, g_sam_the, s_sam_the);
name_ref= sprintf('ref_a%.1fg%.2fs%.4f', a_ref_the, g_ref_the, s_ref_the);

SAM     = load(fullfile(pathData, name_sam + ".mat"));
REF     = load(fullfile(pathData, name_ref + ".mat"));

j_sam = 1;
j_ref = 1;
alpha_ref = a_ref_the; % [dB/cm/MHz]
alpha_sam = a_sam_the; % [dB/cm/MHz]

% Priors
delta_alpha_prior = alpha_sam - alpha_ref; % [dB/cm/MHz]
delta_g_prior     = log(g_sam_the/g_ref_the);  % ln (g_s / g_ref)
delta_s_prior     = s_sam_the - s_ref_the;   % s_s - s_ref



%% SET DATA


% B-MODE CHECK
bmode_sam = db(abs(hilbert(SAM.rf)));
bmode_sam = bmode_sam - max(bmode_sam(:));

bmode_ref = db(abs(hilbert(REF.rf)));
bmode_ref = bmode_ref - max(bmode_ref(:));

% 
% figure,
% 
% subplot(121), 
% imagesc(SAM.x*1E2, SAM.z*1E2, bmode_sam, range_bmode), axis("image");
% xlabel('Lateral [cm]'), ylabel('Axial [cm]');
% title('SAM')
% colormap('gray')
% 
% subplot(122), 
% imagesc(REF.x*1E2, REF.z*1E2, bmode_ref,range_bmode), axis("image");
% xlabel('Lateral [cm]'), ylabel('Axial [cm]');
% title('REF')
% colormap('gray')

%% POWER LAW METHOD PARAMETERS

pars.bw          = [3 9]; % [MHz]
pars.overlap     = 0.8;
pars.blocksize   = 10; % wavelengths

pars.z_roi       = [1.0 4.5]*1E-2; % [m] 
pars.x_roi       = [-1.5 1.5]*1E-2; % [m] 

pars.saran_layer = false;

pars.window_type = 3; %  (1) Hanning, (2) Tuckey, (3) Hamming, (4) Tchebychev

%% REGULARIZATION PARAMETERS

% Implementation parameters
par_rpl.ini_method = 1; % METHOD LEAST SQUARES INITIALIZATION 
par_rpl.ini_tol = 1e-16;
par_rpl.df_op = 0; %(1) 1/2 f(x+1)-f(x-1), (0) f(x+1)-f(x)

par_rpl.tol   = 1e-16;
par_rpl.kmax  = 100;
par_rpl.eps_f = 1e-16;

% Robust (optional)
par_rpl.m_est   = 1;
par_rpl.c       = [1, 100, NaN, 100]; % c_fid, c_b, c_n, c_a
par_rpl.sigma   = [0.1, 0.1, NaN, 0.1]; % sigma_fid, sigma_b, sigma_n, sigma_a

% Parameters for RPL-TV
%     mu_b  = 1E0-1E1; % ORIGINALS
%     mu_n  = 1E3-1E5; % ORIGINALS
%     mu_a  = 1E3-1E5; % ORIGINALS

% mu_rpl_tv    = [10^1.75, 1E4, 1E4]; % [mu_b, mu_n, mu_a] % Test regular    
% mu_rpl_tv    = [1E3, 1E4, 1E4]; % [mu_b, mu_n, mu_a] % Test regular    
mu_rpl_tv    = [1e3, 1e5, 1e5]; % test PowerLaw simu results %EMZ

%% POWER SPECTRA ESTIMATION

% spectralData_sam = calc_powerSpectra(SAM, pars);
spectralData_sam = calc_powerSpectra_vSimple(SAM, pars); % @
S_sam = spectralData_sam.powerSpectra;

% spectralData_ref = calc_powerSpectra(REF, pars);
spectralData_ref = calc_powerSpectra_vSimple(REF, pars); % @
S_ref = spectralData_ref.powerSpectra;


SR_emz = S_sam ./ S_ref;

SR = permute(SR_emz,[3,1,2]);

%% COMPENSATE GAUSS ATTEMPT 2 DoF-a
if strcmp( estim_method, '2DoF-a')

band    = spectralData_sam.band;
depth   = spectralData_sam.depth;
[r,p,q] = size(SR);
  
comp_ref    = comp_ref_a(delta_alpha_prior,j_ref,band,depth,q);
comp_freq_a = comp_mod_freq_a(alpha_ref,j_sam,j_ref,band,depth,q);

SR_comp = SR .* comp_ref .* comp_freq_a;

% indices initialization
f = band(:); % always column vector
% [r,p,q] = size(SR_comp);

% log-spectrum Ratio Y_a = X.g + Z.s 
Y = log(SR_comp);

% matrices for RPL-based algorithms
X = kron( speye(p*q), ones(size(f)) );
Z = kron( speye(p*q), f.^2 );

% initialization for RPL-based methods
u_0 = initialize_rpl_a_prior(Y, X, Z, mu_rpl_tv, par_rpl);
    
% RPL estimation
[u_opt,~] = rpl_tv_a_prior(Y, X, Z, mu_rpl_tv, u_0, par_rpl);

dy = 0.5*(diag(ones(p-1,1),1) - diag(ones(p-1,1),-1));
dy(1,1) = -1; dy(1,2) = 1; dy(end,end) = 1; dy(end,end-1) = -1;
Dy = sparse(kron(speye(q),dy));

% \Deltas
g = u_opt(1:p*q);
s = u_opt(p*q+1:2*p*q);

% Prior "a"
a = delta_alpha_prior*ones(p*q, 1);

% utils 
z = 1E2*repmat(depth,1,q); % 1E2*spectralData_sam.depth * ones(1, q); % 2d array
dz = reshape(Dy*z(:),p,q);
dz(end,:) = dz(end-1,:);

%% COMPENSATE GAUSS ATTEMPT 2 DoF-s
elseif strcmp( estim_method, '2DoF-s')

band    = spectralData_sam.band;
depth   = spectralData_sam.depth;
[r,p,q] = size(SR);
  
comp_ref    = comp_ref_s_bsc(delta_s_prior, band, p, q);
comp_freq_a = comp_mod_freq_a(alpha_ref,j_sam,j_ref,band,depth,q);

SR_comp = SR .* comp_ref .* comp_freq_a;

% indices initialization
f = band(:); % always column vector
% [r,p,q] = size(SR_comp);

% log-spectrum Ratio Y_s = X.g + W.a 
Y = log(SR_comp);

% matrices for RPL-based algorithms
X = kron( speye(p*q), ones(size(f)) );
W = kron( speye(p*q), -4*f );

% initialization for RPL-based methods
u_0 = initialize_rpl_n_prior(Y, X, W, mu_rpl_tv, par_rpl);
    
% RPL estimation
[u_opt,~] = rpl_tv_n_prior(Y, X, W, mu_rpl_tv, u_0, par_rpl);

dy = 0.5*(diag(ones(p-1,1),1) - diag(ones(p-1,1),-1));
dy(1,1) = -1; dy(1,2) = 1; dy(end,end) = 1; dy(end,end-1) = -1;
Dy = sparse(kron(speye(q),dy));

% \Deltas
g = u_opt(1:p*q);
a = u_opt(p*q+1:2*p*q);

% Prior "s"
s = delta_s_prior*ones(p*q, 1);

% utils 
z = 1E2*repmat(depth,1,q); % 1E2*spectralData_sam.depth * ones(1, q); % 2d array
dz = reshape(Dy*z(:),p,q);
dz(end,:) = dz(end-1,:);    

%% COMPENSATE GAUSS ATTEMPT 2 DoF-g
elseif strcmp( estim_method, '2DoF-g')

band    = spectralData_sam.band;
depth   = spectralData_sam.depth;
[r,p,q] = size(SR);
  
comp_ref    = comp_ref_k_bsc(delta_g_prior);
comp_freq_a = comp_mod_freq_a(alpha_ref,j_sam,j_ref,band,depth,q);

SR_comp = SR .* comp_ref .*comp_freq_a;

% indices initialization
f = band(:); % always column vector
% [r,p,q] = size(SR_comp);

% log-spectrum Ratio Y_g = Z.s + W.a 
Y = log(SR_comp);

% matrices for RPL-based algorithms
Z = kron( speye(p*q), f.^2 );
W = kron( speye(p*q), -4*f );

% initialization for RPL-based methods
u_0 = initialize_rpl_b_prior(Y, Z, W, mu_rpl_tv, par_rpl);
    
% RPL estimation
[u_opt,~] = rpl_tv_b_prior(Y, Z, W, mu_rpl_tv, u_0, par_rpl);

dy = 0.5*(diag(ones(p-1,1),1) - diag(ones(p-1,1),-1));
dy(1,1) = -1; dy(1,2) = 1; dy(end,end) = 1; dy(end,end-1) = -1;
Dy = sparse(kron(speye(q),dy));

% \Deltas
s = u_opt(1:p*q);
a = u_opt(p*q+1:2*p*q);

% Prior "g"
g = delta_g_prior*ones(p*q, 1);

% utils 
z = 1E2*repmat(depth,1,q); % 1E2*spectralData_sam.depth * ones(1, q); % 2d array
dz = reshape(Dy*z(:),p,q);
dz(end,:) = dz(end-1,:);        

%% COMPENSATE GAUSS ATTEMPT 3 DoF
elseif strcmp( estim_method, '3DoF')

band    = spectralData_sam.band;
depth   = spectralData_sam.depth;
[r,p,q] = size(SR);

comp_freq_a = comp_mod_freq_a(alpha_ref,j_sam,j_ref,band,depth,q);

SR_comp = SR.*comp_freq_a;

% indices initialization
f = band(:); % always column vector
% [r,p,q] = size(SR_comp);

% log-spectrum Ratio Y = X.g + Z.s + W.a
Y = log(SR_comp);

% matrices for RPL-based algorithms
X = kron( speye(p*q), ones(size(f)) );
Z = kron( speye(p*q), f.^2 );
W = kron( speye(p*q), -4*f );

% initialization for RPL-based methods
u_0 = initialize_rpl(Y, X, Z, W, mu_rpl_tv, par_rpl);
    
% RPL estimation
[u_opt,~] = rpl_tv(Y, X, Z, W, mu_rpl_tv, u_0, par_rpl);

dy = 0.5*(diag(ones(p-1,1),1) - diag(ones(p-1,1),-1));
dy(1,1) = -1; dy(1,2) = 1; dy(end,end) = 1; dy(end,end-1) = -1;
Dy = sparse(kron(speye(q),dy));

% \Deltas
g = u_opt(1:p*q);
s = u_opt(p*q+1:2*p*q);
a = u_opt(2*p*q+1:3*p*q);

% utils 
z = 1E2*repmat(depth,1,q); % 1E2*spectralData_sam.depth * ones(1, q); % 2d array
dz = reshape(Dy*z(:),p,q);
dz(end,:) = dz(end-1,:);


end
%%

%% QUS PARAMETERS 

plotBSCdB   = true;

g_ratio     = reshape(exp(g), p, q);
g_ratio_dB  = 10*log10(g_ratio);
alpha_ratio = reshape(a, p, q);
s_ratio     = reshape(s, p, q); 

% METRICS 3 DOF
acs_sam   = alpha_ratio + alpha_ref;

mean2d = @(x) mean(x(:));
std2d = @(x) std(x(:));
cv2d = @(x) 100*std(x(:))/mean(x(:));

calc2dStats = {@(x) mean(x(:)), @(x) std(x(:)), @(x) 100 * std(x(:)) / mean(x(:))};

[m_g, s_g, cv_g] = deal(calc2dStats{1}(g_ratio), calc2dStats{2}(g_ratio), calc2dStats{3}(g_ratio));
if plotBSCdB 
    [m_g, s_g, cv_g] = deal(calc2dStats{1}(g_ratio_dB), calc2dStats{2}(g_ratio_dB), calc2dStats{3}(g_ratio_dB));
end

[m_s, s_s, cv_s] = deal(calc2dStats{1}(s_ratio), calc2dStats{2}(s_ratio), calc2dStats{3}(s_ratio));
[m_a, s_a, cv_a] = deal(calc2dStats{1}(acs_sam), calc2dStats{2}(acs_sam), calc2dStats{3}(acs_sam));

fprintf('----%s----\n', estim_method);
fprintf('Δα         : %.3f ± %.4f, %%CV = %.4f\n', round(mean2d(alpha_ratio), 3), round(std2d(alpha_ratio), 4), round(cv2d(alpha_ratio), 4));
fprintf('α_s        : %.3f ± %.4f, %%CV = %.4f\n', round(m_a, 3), round(s_a, 4), round(cv_a, 4));
    if plotBSCdB 
fprintf('g_s/g_r[dB]: %.3f ± %.4f, %%CV = %.4f\n', round(m_g, 3), round(s_g, 4), round(cv_g, 4));
    else
fprintf('g_s/g_r    : %.3f ± %.4f, %%CV = %.4f\n', round(m_g, 3), round(s_g, 4), round(cv_g, 4));
    end

fprintf('Δs         : %.4f ± %.4f, %%CV = %.4f\n', round(m_s, 4), round(s_s, 4), round(cv_s, 4));
fprintf('------------\n');

%% IMAGESC PLOTS
Xaxis = spectralData_ref.lateral;
Zaxis = spectralData_ref.depth;
cm = 1e2;

axis_n = [0 1.2];
axis_a = [0 3];
axis_b = [-60 0]; % dB
fontSize = 16;

figure, 
set(gcf,'units','normalized','outerposition',[0 0.15 1 0.75]); box on;
sgtitle(estim_method, 'FontSize', fontSize+2, 'FontWeight', 'bold');

subplot(1,3,1)
imagesc(Xaxis*cm, Zaxis*cm, acs_sam), colorbar
axis("image");
xlabel('Lateral [cm]'), ylabel('Depth [cm]'), colormap turbo;
% title('\Delta \alpha ');
% title(['ACS: ', num2str(round(m_a, 3)), ' \pm ', num2str(round(s_a, 2)), ', CV = ', num2str(round(cv_a, 3))])
title({'$\alpha_s$:', ...
       [num2str(round(m_a, 3)), ' $\pm$ ', num2str(round(s_a, 3)), ', CV = ', num2str(round(cv_a, 3))]}, ...
      'Interpreter', 'latex');

h2 = colorbar; 
ylabel(h2,'dB\cdotcm^{-1}\cdotMHz^{-1}','FontSize', fontSize);
set(gca,'fontsize',fontSize)

subplot(1,3,2)
imagesc(Xaxis*cm, Zaxis*cm, g_ratio)
h2 = colorbar; 
if plotBSCdB 
   imagesc(Xaxis*cm, Zaxis*cm, g_ratio_dB)
   h2 = colorbar;
   ylabel(h2,'dB','FontSize', fontSize);
end
axis("image");
xlabel('Lateral [cm]'), colormap turbo;
% title(['$\frac{g_s}{g_r}: ', num2str(round(m_g, 3)), ' \pm ', num2str(round(s_g, 2)), ', CV = ', num2str(round(cv_g, 3)), '$'], ...
%       'Interpreter', 'latex')
title({'$\frac{g_s}{g_r}$:', ...
       [num2str(round(m_g, 3)), ' $\pm$ ', num2str(round(s_g, 3)), ', CV = ', num2str(round(cv_g, 3))]}, ...
      'Interpreter', 'latex');

set(gca,'fontsize',fontSize)

subplot(1,3,3)
imagesc(Xaxis*cm, Zaxis*cm, s_ratio), colorbar
axis("image");
xlabel('Lateral [cm]'), colormap turbo;
% title(['$\Delta s$: ', num2str(round(m_s, 3)), ' $\pm$ ', num2str(round(s_s, 3)), ', CV = ', num2str(round(cv_s, 3))], ...
%       'Interpreter', 'latex');
title({'$\Delta s$:', ...
       [num2str(round(m_s, 3)), ' $\pm$ ', num2str(round(s_s, 3)), ', CV = ', num2str(round(cv_s, 3))]}, ...
      'Interpreter', 'latex');
set(gca,'fontsize',fontSize)

%% FIGURE INTERP OVERLAY BMODE, DELTA SNR, ACS, DELTA BSC, DELTA N

fontSize = 16;

figure,
set(gcf,'units','normalized','outerposition',[0 0.1 1 0.8]); box on;

tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
sgtitle(estim_method, 'FontSize', fontSize+2, 'FontWeight', 'bold');

%%%%%%%%%%%%%%%%%%%%%%%%%% alpha_s (ACS) %%%%%%%%%%%%%%%%%%%%%%%%%%
acs_sam = alpha_ratio + alpha_ref;

units           = 1E3;
bmodeFull       = bmode_sam;
colorImg        = acs_sam;
range_bmode     = [-60 0];
range_img       = [0.1 1.2];
transparency    = 0.65;
x_img           = spectralData_sam.lateral*units;
z_img           = spectralData_sam.depth*units;
xFull           = SAM.x*units;
zFull           = SAM.z*units;
[X, Z] = meshgrid(xFull, zFull);
roi = and(X >= x_img(1), X <= x_img(end)) & ...
      and(Z >= z_img(1), Z <= z_img(end));

t = nexttile;
    [~,hB,hColor] = imOverlayInterp(bmodeFull, colorImg, range_bmode, range_img, ...
                        transparency, x_img, z_img, roi, xFull, zFull);   
    hold on;
    contour(xFull, zFull, roi, 1,'w--', 'LineWidth', 2)
    hold off;
    xlabel('Lateral'), ylabel('Depth');
    hColor.Label.String = 'dB\cdotcm^{-1}\cdotMHz^{-1}';
    title('$\alpha_s$', 'Interpreter', 'latex')
    set(gca,'fontsize',fontSize)
%%%%%%%%%%%%%%%%%%%%%%%%%% Delta alpha (ACS) %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% Delta g (BSC) %%%%%%%%%%%%%%%%%%%%%%%%%%

units           = 1E3;
bmodeFull       = bmode_sam;
colorImg        = g_ratio;
range_bmode     = [-60 0];
range_img       = [];
transparency    = 0.65;
x_img           = spectralData_sam.lateral*units;
z_img           = spectralData_sam.depth*units;
xFull           = SAM.x*units;
zFull           = SAM.z*units;
[X, Z] = meshgrid(xFull, zFull);
roi = and(X >= x_img(1), X <= x_img(end)) & ...
      and(Z >= z_img(1), Z <= z_img(end));

    if plotBSCdB 
       colorImg = g_ratio_dB;
    end

t = nexttile;
    [~,hB,hColor] = imOverlayInterp(bmodeFull, colorImg, range_bmode, range_img, ...
                        transparency, x_img, z_img, roi, xFull, zFull);
    hold on;
    contour(xFull, zFull, roi, 1,'w--', 'LineWidth', 2)
    hold off;
    xlabel('Lateral'), ylabel('Depth');
    hColor.Label.String = '';
        if plotBSCdB 
            hColor.Label.String ='dB';
        end
    title('$\frac{g_s}{g_r}$', 'Interpreter','latex')
    set(gca,'fontsize',fontSize)
%%%%%%%%%%%%%%%%%%%%%%%%%% Delta g ratio in dB %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% Delta s %%%%%%%%%%%%%%%%%%%%%%%%%%

units           = 1E3;
bmodeFull       = bmode_sam;
colorImg        = s_ratio;
range_bmode     = [-60 0];
range_img       = [];
transparency    = 0.65;
x_img           = spectralData_sam.lateral*units;
z_img           = spectralData_sam.depth*units;
xFull           = SAM.x*units;
zFull           = SAM.z*units;
[X, Z] = meshgrid(xFull, zFull);
roi = and(X >= x_img(1), X <= x_img(end)) & ...
      and(Z >= z_img(1), Z <= z_img(end));


t = nexttile;
    [~,hB,hColor] = imOverlayInterp(bmodeFull, colorImg, range_bmode, range_img, ...
                        transparency, x_img, z_img, roi, xFull, zFull);
    hold on;
    contour(xFull, zFull, roi, 1,'w--', 'LineWidth', 2)
    hold off;
    xlabel('Lateral'), ylabel('Depth');
    hColor.Label.String = '';
    title('$\Delta s$', 'Interpreter','latex')
    set(gca,'fontsize',fontSize)
%%%%%%%%%%%%%%%%%%%%%%%%%% Delta s %%%%%%%%%%%%%%%%%%%%%%%%%%