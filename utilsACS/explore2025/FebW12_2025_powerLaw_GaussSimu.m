% Test of calc_powerSpectra function (simplified Chahuara code)
% USING DATA SIMU KWAVE GAUSS ATTEMPT TO CHANGE BSC CHECK KERNEL VARIABLES
% IN SIMULATION
% ATTENUATION alpha_power 1.001, 1.005, 1.05, 1.1
% OPTION AVAILABLE 2DoF with priors and 3-DoF by changing estim_method
% Dec 2025 EMZ

%%
% clear all, 
% clc, 
warning('off');
% close all;

addpath(genpath(pwd))

Np2dB       = 20*log10(exp(1));
dB2Np       = 1/Np2dB;
range_bmode = [-60 0];

plotBmode   = false;
plotBSCdB   = true;


methods = {'3-DoF', '2-DoF-a', '2-DoF-g', '2-DoF-s'};
% methods = {'2-DoF-s'};

% First row for headers, second for data
bsc_results = cell(2, length(methods)); 
acs_results = cell(2, length(methods));

% Store headers
bsc_results(1, :) = methods;
acs_results(1, :) = methods;
%% LOAD DATA
pathData = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\dataACS_kwave';

% DATA NEW AMZ
alpha_sam = 0.7; 
j_sam = 1.1;

alpha_ref = 0.5;
j_ref = j_sam;

folderDataSam = sprintf('sim_TUFFC25\\acs%.1f_pow%.1f', alpha_sam, j_sam);
folderDataSam = strrep(folderDataSam, '.', 'p');
% folderDataSam = 'sim_TUFFC25\acs1_pow1p1';

rf_sam_name = strcat('rfsam_', sprintf('%.3f', j_sam));
rf_sam_name = strrep(rf_sam_name, '.', 'p');
rf_sam_name = strcat(rf_sam_name, '.mat');
% SAM = load(fullfile(pathData, folderData, rf_sam_name)); %*old*
SAM = load(fullfile(pathData, folderDataSam, rf_sam_name));

SAM.alpha_power = j_sam;
SAM.acs = alpha_sam; % [dB/cm/MHz] 

folderDataRef = sprintf('sim_TUFFC25\\acs%.1f_pow%.1f', alpha_ref, j_sam);
folderDataRef = strrep(folderDataRef, '.', 'p');
% folderDataRef = 'sim_TUFFC25\acs0p5_pow1p1';

rf_ref_name = strcat('rfref_', sprintf('%.3f', j_ref));
rf_ref_name = strrep(rf_ref_name, '.', 'p');
rf_ref_name = strcat(rf_ref_name, '.mat');
% REF = load(fullfile(pathData, folderData, rf_ref_name)); %*old*
REF = load(fullfile(pathData, folderDataRef, rf_ref_name)); 

REF.alpha_power = j_ref; 
REF.acs = alpha_ref; % [dB/cm/MHz]

% FROM RPM (GROUNTRUTH) ESTIMATE DELTA_S_PRIOR
switch j_sam 

    case 1.1
        delta_alpha_prior = alpha_sam - alpha_ref; % [dB/cm/MHz]
        delta_g_prior     = log(1);  % ln (g_sam / g_ref) 0dB
        delta_s_prior     = 0.0245; % s_sam - s_ref og BSC
        % delta_s_prior     = 0.022; % s_sam - s_ref og ACS
    case 1.3 % TBD
        delta_alpha_prior = alpha_sam - alpha_ref; % [dB/cm/MHz]
        delta_g_prior     = log(1);  % ln (g_sam / g_ref) 0.9
        delta_s_prior     = 0.021; % s_sam - s_ref og
end

% B-MODE CHECK

bmode_sam = db(abs(hilbert(SAM.rf)));
bmode_sam = bmode_sam - max(bmode_sam(:));

bmode_ref = db(abs(hilbert(REF.rf)));
bmode_ref = bmode_ref - max(bmode_ref(:));

%% SPECTRAL METHOD PARAMETERS
pars.P = 4096; % NFFT only for calculate BSC_RPM_ok 10wl
% pars.P = 8192; % 15wl
pars.bw          = [3 8.5]; % [MHz]
pars.overlap     = 0.8;
pars.blocksize   = 12; % wavelengths
% pars.z_roi       = [5 25]*1E-3; % [m]  % half
pars.z_roi       = [5 40]*1E-3; % all
pars.z_roi       = [5 45]*1E-3; % all **
pars.x_roi       = [-18 18]*1E-3; % [m] % all

pars.window_type = 3; %  (1) Hanning, (2) Tuckey, (3) Hamming, (4) Tchebychev
pars.saran_layer = true;

if (plotBmode)
figure,

subplot(121), 
imagesc(SAM.x*1E3, SAM.z*1E3, bmode_sam, range_bmode), axis("image"), hold on;
rectangle('Position', 1E3*[pars.x_roi(1) pars.z_roi(1) pars.x_roi(2)-pars.x_roi(1) pars.z_roi(2)-pars.z_roi(1)], ...
        'EdgeColor','w', 'LineWidth', 2, 'LineStyle','--'), hold off;
clim(range_bmode)
cb = colorbar;
cb.Label.String = 'dB'; % Add the label "dB"
xlabel('Lateral [mm]'), ylabel('Depth [mm]');
title('SAM')
colormap('gray')

subplot(122), 
imagesc(REF.x*1E3, REF.z*1E3, bmode_ref, range_bmode), axis("image");
rectangle('Position', 1E3*[pars.x_roi(1) pars.z_roi(1) pars.x_roi(2)-pars.x_roi(1) pars.z_roi(2)-pars.z_roi(1)], ...
        'EdgeColor','w', 'LineWidth', 2, 'LineStyle','--'), hold off;
clim(range_bmode)
cb = colorbar;
cb.Label.String = 'dB'; % Add the label "dB"
xlabel('Lateral [mm]'), ylabel('Depth [mm]');
title('REF')
colormap('gray')
end

%%

%% POWER SPECTRA ESTIMATION
% spectralData_sam = calc_powerSpectra(SAM, pars);
spectralData_sam = calc_powerSpectra_vSimple(SAM, pars); % @
S_sam = spectralData_sam.powerSpectra;

num_ref = 1;

% spectralData_ref = calc_powerSpectra(REF, pars);
spectralData_ref = calc_powerSpectra_vSimple(REF, pars); % @
S_ref = spectralData_ref.powerSpectra;

% if num_ref > 1
%     S_ref_aux = zeros([size(S_sam), num_ref]);
% 
%     for iRef = 1:num_ref
%         spectralData_ref = calc_powerSpectra(REF, pars);
%         S_ref_aux(:,:,:,iRef) = spectralData_ref.powerSpectra;    
%     end
%     S_ref = mean(S_ref_aux, 4);
% end

%%
SR_emz = S_sam ./ S_ref;

SR = permute(SR_emz,[3,1,2]);

%% EXTRA $ EMZ TO PLOT Spectrum at given position
% z_target = 1.5E-3;
% x_target = 0E-3;
% 
% [~, i_idx] = min(abs(spectralData_sam.depth   - z_target)); % $
% [~, j_idx] = min(abs(spectralData_sam.lateral - x_target)); % $
% 
%  SS_idx = squeeze(S_sam(i_idx,j_idx,:));
%  SR_idx = squeeze(S_ref(i_idx,j_idx,:));
% 
% figure, 
% plot(spectralData_sam.band, SS_idx, 'DisplayName', 'Ssam'), hold on, grid on
% plot(spectralData_sam.band, SR_idx, 'DisplayName', 'Sref'), hold off
% legend('Location','best');

%% GENERAL REGULARIZTATION SETTINGS
% Implementation parameters
par_rpl.tol        = 1e-16;
par_rpl.kmax       = 100;
par_rpl.eps_f      = 1e-16;
par_rpl.m_est      = 0; %Robust
par_rpl.ini_tol    = 1e-16;
par_rpl.df_op      = 1;
par_rpl.ini_method = 1; % METHOD LEAST SQUARES INITIALIZATION 

% Parameters for RPL-TV
mu_rpl_tv    = [1E3; 1E3; 1E4]; % [mu_g, mu_s, mu_a]
% mu_rpl_tv    = [0.01; 0.01; 0.01]; % [mu_g, mu_s, mu_a]
% mu_rpl_tv    = [0.001; 0.001; 0.001]; % [mu_g, mu_s, mu_a]

%% FOR BUCLE
for iMet = 1:length(methods)

estim_method = methods{iMet};

%% COMPENSATE GAUSS ATTEMPT 2-DoF-a
if strcmp( estim_method, '2-DoF-a')

band    = spectralData_sam.band;
depth   = spectralData_sam.depth;
[r,p,q] = size(SR);
  
comp_ref    = comp_ref_a(-delta_alpha_prior,j_ref,band,depth,q);
comp_freq_a = comp_mod_freq_a(alpha_ref,j_sam,j_ref,band,depth,q);

SR_comp = SR .* comp_ref .* comp_freq_a;

% indices initialization
f = band(:); % always column vector
% [r,p,q] = size(SR_comp);

% log-spectrum Ratio Y_a = X.g + Z.s 
Y = log(SR_comp);

% matrices for RPL-based algorithms
X = kron( speye(p*q), ones(size(f)) );
Z = kron( speye(p*q), -f.^2 );

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

% Prior "a" known
a_Np2dB = delta_alpha_prior*ones(p*q, 1);

% utils 
z = 1E2*repmat(depth,1,q); % 1E2*spectralData_sam.depth * ones(1, q); % 2d array
dz = reshape(Dy*z(:),p,q);
dz(end,:) = dz(end-1,:);

%% COMPENSATE GAUSS ATTEMPT 2-DoF-s
elseif strcmp( estim_method, '2-DoF-s')

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
s = +delta_s_prior*ones(p*q, 1);

% utils 
z = 1E2*repmat(depth,1,q); % 1E2*spectralData_sam.depth * ones(1, q); % 2d array
dz = reshape(Dy*z(:),p,q);
dz(end,:) = dz(end-1,:);  

a_Np2dB = Np2dB*Dy*a./dz(:);

%% COMPENSATE GAUSS ATTEMPT 2-DoF-g
elseif strcmp( estim_method, '2-DoF-g')

band    = spectralData_sam.band;
depth   = spectralData_sam.depth;
[r,p,q] = size(SR);
  
comp_ref    = comp_ref_g_bsc(delta_g_prior);
comp_freq_a = comp_mod_freq_a(alpha_ref,j_sam,j_ref,band,depth,q);

SR_comp = SR .* comp_ref .*comp_freq_a;

% indices initialization
f = band(:); % always column vector
% [r,p,q] = size(SR_comp);

% log-spectrum Ratio Y_g = Z.s + W.a 
Y = log(SR_comp);

% matrices for RPL-based algorithms
Z = kron( speye(p*q), -f.^2 ); % EMZ
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

a_Np2dB = Np2dB*Dy*a./dz(:);

%% COMPENSATE GAUSS ATTEMPT 3-DoF
elseif strcmp( estim_method, '3-DoF')

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
Z = kron( speye(p*q), -f.^2 ); % EMZ
% W = kron( speye(p*q), -4*f.^j_sam );
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

a_Np2dB = Np2dB*Dy*a./dz(:);

end
%%
%% QUS PARAMETERS 

g_ratio     = reshape(exp(g), p, q);
g_ratio_dB  = 10*log10(g_ratio);
alpha_ratio = reshape(a_Np2dB, p, q);
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

fprintf('-----%s---\n', estim_method);
fprintf('α_s        : %.3f ± %.4f, %%CV = %.4f\n', round(m_a, 3), round(s_a, 4), round(cv_a, 4));
    if plotBSCdB 
fprintf('g_s/g_r[dB]: %.3f ± %.4f, %%CV = %.4f\n', round(m_g, 3), round(s_g, 4), round(cv_g, 4));
    else
fprintf('g_s/g_r    : %.3f ± %.4f, %%CV = %.4f\n', round(m_g, 3), round(s_g, 4), round(cv_g, 4));
    end

fprintf('Δs         : %.4f ± %.4f, %%CV = %.4f\n', round(m_s, 4), round(s_s, 4), round(cv_s, 4));
fprintf('--------\n');
 
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
if plotBmode
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
range_bmode     = [-100 0];
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
end
%%
% keyboard

%% BSC RECONSTRUCTION DATA AC (change later)
freq = spectralData_sam.band;

% g_est = mean(g_ratio(:));
% s_est = mean(s_ratio(:));

g_est = median(g_ratio(:));
s_est = median(s_ratio(:));

bsc_est_gauss = g_est.*exp(-s_est.*freq.^2);

% Perform linear regression  BSC GAUSS TO POWER LAW 
% ln(bsc) = d_n . ln(f) + ln(d_b) 


% % y = [log(freq), I]. [d_n, ln(d_b)]
% M_pl = [ones(length(freq),1),log(freq')];
% r_pl = log(bsc_est_gauss);
% x_pl = cgs(M_pl'*M_pl,M_pl'*r_pl,1e-16,10000); % coeffs_pl
% d_b_pl = exp(x_pl(1));
% d_n_pl = x_pl(2);


coeffs_pl   = polyfit(log(freq), log(bsc_est_gauss), 1); % Fit y = m.x + c
d_n_pl      = coeffs_pl(1); % Slope = d_n  (mean),  (median)
ln_d_b_pl   = coeffs_pl(2); % Intercept = ln(d_b) 
d_b_pl      = exp(ln_d_b_pl); % 1.0917 (mean), 0.9079(median)

bsc_fit_powlaw = d_b_pl*freq.^d_n_pl;

% Display results
fprintf('-----PowerLaw (b.f^n)-----\n')
fprintf('Δn           = %f\n', d_n_pl);
fprintf('d_b          = %f\n', d_b_pl);
fprintf('b_s/b_r [dB] = %f\n', 10*log10(d_b_pl));
fprintf('---------\n')
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RPM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% NAME BSC
nameBSC = strcat('BSC_', sprintf('sam%.2f_ref%.2f', alpha_sam, alpha_ref ));
nameBSC = strrep(nameBSC, '.', 'p');
% BSC = calculateBSC_RPM_ok(SAM, REF, pars); % slow only once**
% BSC = calculateBSC_RPM_fast(SAM, REF, pars); % fast TBD**
% save(fullfile(pathData, folderDataSam, nameBSC), "BSC")
load(fullfile(pathData, folderDataSam, nameBSC));
%%freq = BSC.band;
% bsc_rpm = BSC.BSCcurve_Uni(:,1); % mean
bsc_rpm = BSC.BSCcurve_Uni(:,2); % median

% MM = [freq.^2, ones(size(freq))];
% rr = log(bsc_rpm);
% pp = cgs(MM'*MM,MM'*rr,1e-16,10000); % coeffs

% Perform linear regression  bsc = d_s . f^2 + ln(d_g) 
coeffs   = polyfit(-freq.^2, log(bsc_rpm), 1); % Fit y = mx + c
d_s      = coeffs(1); % Slope = d_s -0.0319 (mean), -0.0317 (median)
ln_dg    = coeffs(2); % Intercept = ln(d_g) 
d_g      = exp(ln_dg); % 1.0917 (mean), 0.9079(median)

% Display results
fprintf('-----RPM Gauus (g.exp(-s.f^2))-----\n')
fprintf('Δs           = %f\n', d_s);
fprintf('d_g          = %f\n', d_g);
fprintf('g_s/g_r [dB] = %f\n', 10*log10(d_g));
fprintf('---------\n')

bsc_rpm_gauss = d_g*exp(-d_s* freq.^2);

%% FIGURES

% figure
% % figure (101), 
% plot(freq, bsc_rpm, 'b-', 'DisplayName', 'RPM'), 
% hold on;
% plot(freq, bsc_rpm_gauss, 'b--', 'DisplayName', 'RPM Fit Gauss');
% % plot(freq, bsc_fit_powlaw, 'g--', 'DisplayName', 'Power Law Fit');
% plot(freq, bsc_est_gauss, 'r--', 'DisplayName', estim_method);
% grid on;
% % hold off;
% xlabel('Frequency [MHz]'), 
% ylabel('BSC [cm^{-1}\cdot sr^{-1}]')
% title('BSC')
% legend ('Location', 'Best')
% set(gca,'fontsize',fontSize)

figure,
% figure (102), 
loglog(freq, bsc_rpm, 'b-', 'DisplayName', 'RPM'), 
hold on;
loglog(freq, bsc_rpm_gauss, 'b--', 'DisplayName', 'RPM Fit Gauss');
% loglog(freq, bsc_fit_powlaw, 'g--', 'DisplayName', 'Power Law Fit');
loglog(freq, bsc_est_gauss, 'r--', 'DisplayName', estim_method);
grid on;
% hold off;
xlabel('Frequency [MHz]'), 
ylabel('BSC [cm^{-1}\cdot sr^{-1}]')
ylim([0.1 1.1])
title('BSC')
legend ('Location', 'Best')
set(gca,'fontsize',fontSize)

% figure,
% % figure(103),
% plot(freq, 10*log10(bsc_rpm), 'b-', 'DisplayName', 'RPM'), 
% hold on;
% plot(freq, 10*log10(bsc_rpm_gauss), 'b--', 'DisplayName', 'RPM Fit Gauss');
% % plot(freq, 10*log10(bsc_fit_powlaw),'g--', 'DisplayName', 'Power Law Fit');
% plot(freq, 10*log10(bsc_est_gauss), 'r--', 'DisplayName', estim_method);
% grid on;
% % hold off;
% xlabel('Frequency [MHz]'), 
% ylabel('BSC [dB]')
% title('BSC')
% legend ('Location', 'Best')
% set(gca,'fontsize',fontSize)
%% SAVE ALL BSC EST GAUSS

bsc_results{2, iMet} = bsc_est_gauss;

%% SAVE ALL ACS MAPS

acs_results{2, iMet} = acs_sam;


end
keyboard
%% PLOTS BSC TOGETHER
xlim_range = [3.0 8.51];
bsc_results{1, 1} = '3-DoF';
bsc_results{1, 2} = '2-DoF "a"';
bsc_results{1, 3} = '2-DoF "b"';
bsc_results{1, 4} = '2-DoF "n" ';

% Define properties for customization
line_width = 2.85; % Set line width
font_size = 30; % Adjust font size
ylim_range = [0.1 1.05]; % Y-axis limits

% Convert hexadecimal colors to RGB (MATLAB requires values between 0 and 1)
color_rpm = '#000000'; % Black
color_1 = '#FF0000';  % 3dof
color_2 = '#D95319';  % 2dof a
color_3 = '#0072BD';  % 2dof b
color_4 = '#77AC30';  % 2dof n

% Convert hex to RGB format
color_rpm = hex2rgb(color_rpm);
color_1 = hex2rgb(color_1); % a
color_2 = hex2rgb(color_2); % b
color_3 = hex2rgb(color_3); % n
color_4 = hex2rgb(color_4);

% Create figure and plot data
figure, 
loglog(freq, bsc_rpm, '-', 'Color', color_rpm, 'LineWidth', line_width+0.5, 'DisplayName', 'RPM');
hold on;
loglog(freq, bsc_results{2, 1}, '--', 'Color', color_1, 'LineWidth', line_width, 'DisplayName', bsc_results{1, 1});
loglog(freq, bsc_results{2, 2}, '--', 'Color', color_2, 'LineWidth', line_width, 'DisplayName', bsc_results{1, 2});
loglog(freq, bsc_results{2, 3}, '--', 'Color', color_3, 'LineWidth', line_width, 'DisplayName', bsc_results{1, 3});
loglog(freq, bsc_results{2, 4}, '--', 'Color', color_4, 'LineWidth', line_width, 'DisplayName', bsc_results{1, 4});

% Customize plot
grid on;
xlabel('Frequency [MHz]', 'FontSize', font_size);
ylabel('BSC [cm^{-1}\cdot sr^{-1}]', 'FontSize', font_size);
ylim(ylim_range);
xlim(xlim_range);
title('Backscatter Coefficient (BSC)', 'FontSize', font_size + 2);
legend('Location', 'southwest');
set(gca, 'FontSize', font_size);
hold off;

% Function to convert hexadecimal to RGB
function rgb = hex2rgb(hex)
    hex = char(hex); % Ensure it's a string
    if hex(1) == '#'
        hex(1) = []; % Remove the '#' symbol
    end
    rgb = reshape(sscanf(hex, '%2x') / 255, 1, 3);
end

%% PLOTS ACS MAPS TOGETHER

acs_results{1, 1} = '3-DoF';
acs_results{1, 2} = '2-DoF "a"';
acs_results{1, 3} = '2-DoF "b"';
acs_results{1, 4} = '2-DoF "n" ';

fontSize = 16;

figure,
% set(gcf,'units','normalized','outerposition',[0 0.1 1 0.8]); box on;

tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
sgtitle(estim_method, 'FontSize', fontSize+2, 'FontWeight', 'bold');

for iMet = 1:length(methods)
%%%%%%%%%%%%%%%%%%%%%%%%%% alpha_s (ACS) %%%%%%%%%%%%%%%%%%%%%%%%%%
acs_sam = alpha_ratio + alpha_ref;

units           = 1E3;
bmodeFull       = bmode_sam;
colorImg        = acs_results{2, iMet};
range_bmode     = [-100 0];
range_img       = [0.1 1.2];
transparency    = 0.65;
x_img           = spectralData_sam.lateral*units;
z_img           = spectralData_sam.depth*units;
xFull           = SAM.x*units;
zFull           = SAM.z*units;
[X, Z] = meshgrid(xFull, zFull);
z_img(end) = 45;
roi = and(X >= x_img(1), X <= x_img(end)) & ...
      and(Z >= z_img(1), Z <= z_img(end));

t = nexttile;
    [~,hB,hColor] = imOverlayInterp(bmodeFull, colorImg, range_bmode, range_img, ...
                        transparency, x_img, z_img, roi, xFull, zFull);   
    hold on;
    contour(xFull, zFull, roi, 1,'w--', 'LineWidth', 2)
    hold off;
    xlabel('Lateral [cm]'), ylabel('Depth [cm]');
    hColor.Label.String = 'dB\cdotcm^{-1}\cdotMHz^{-1}';
    title(bsc_results{1, iMet})
    set(gca,'fontsize',fontSize)

end
%%
