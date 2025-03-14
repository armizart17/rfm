% Test of calc_powerSpectra function (simplified Chahuara code)
% USING DATA SIMU KWAVE GAUSS ATTEMPT TO CHANGE BSC CHECK KERNEL VARIABLES
% IN SIMU
% THE REFERENCE WAS GAUSS MODIFIED TOO
% Nov 2024 EMZ

%%
clear all, clc, 
% close all;

Np2dB       = 20*log10(exp(1));
dB2Np       = 1/Np2dB;
range_bmode = [-60 0];

addpath(genpath(pwd))
%% LOAD DATA
pathData = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\dataACS_kwave\simGauss_v1';

% V2
% SAM = load(fullfile(pathData, 'rf_gauss_v2.mat'));
% REF = load(fullfile(pathData, 'rf_ref_v2.mat'));

% V3
SAM = load(fullfile(pathData, 'rf_gauss_v3.mat'));
REF = load(fullfile(pathData, 'rf_ref_v3.mat'));

% B-MODE CHECK

% bmode_sam = abs(hilbert(SAM.rf));
bmode_sam = db(abs(hilbert(SAM.rf)));
bmode_sam = bmode_sam - max(bmode_sam(:));

bmode_ref = db(abs(hilbert(REF.rf)));
bmode_ref = bmode_ref - max(bmode_ref(:));

%% POWER LAW METHOD PARAMETERS

pars.bw          = [3 9]; % [MHz]
pars.overlap     = 0.8;
pars.blocksize   = 10; % wavelengths
% pars.z_roi       = [5 25]*1E-3; % [m]  % half
pars.z_roi       = [5 50]*1E-3;
pars.x_roi       = [-18 18]*1E-3; % [m] 

pars.window_type = 3; %  (1) Hanning, (2) Tuckey, (3) Hamming, (4) Tchebychev
pars.saran_layer = true;

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
%% COMPENSATE GAUSS ATTEMPT 3 DoF
band    = spectralData_sam.band;
depth   = spectralData_sam.depth;
[r,p,q] = size(SR);
    
% j_sam = 1;
% j_ref = 1;
% a_ref = 0.4; % [dB/cm/MHz]
% 
% comp_freq_a = comp_mod_freq_a(a_ref, j_sam, j_ref, band, depth, q); % 1 when j_sam = j_ref;
% SR_comp = SR.*comp_freq_a;

SR_comp = SR;

% indices initialization
f = band(:); % always column vector
% [r,p,q] = size(SR_comp);

% log-spectrum Ratio Y = X.g + Z.s + W.a
Y = log(SR_comp);

% matrices for RPL-based algorithms
X = kron( speye(p*q), ones(size(f)) );
Z = kron( speye(p*q), f.^2 );
W = kron( speye(p*q), -4*f );

% Implementation parameters
par_rpl.tol        = 1e-16;
par_rpl.kmax       = 100;
par_rpl.eps_f      = 1e-16;
par_rpl.m_est      = 0; %Robust
par_rpl.ini_tol    = 1e-16;
par_rpl.df_op      = 1;
par_rpl.ini_method = 1; % METHOD LEAST SQUARES INITIALIZATION 

% Parameters for RPL-TV
%     mu_b  = 1e+0;%1e+0; %1e+0; % ORIGINALS
%     mu_n  = 1e+3;%1e+3; %1e+2; % ORIGINALS
%     mu_a  = 1e+3;%1e+3; %1e+3; % ORIGINALS

mu_rpl_tv    = [1E1; 1E3; 1E4]; % [mu_k, mu_s, mu_a]

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

%% QUS PARAMETERS

plotBSCdB   = true;

g_ratio     = reshape(exp(g), p, q);
g_ratio_dB  = 10*log10(g_ratio);
alpha_ratio = reshape(Np2dB*Dy*a./dz(:), p, q);
s_ratio     = reshape(s, p, q); 

% METRICS 3 DOF
alpha_ref = 0.4;
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

fprintf('-----3DoF----\n');
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
sgtitle('\bf3DOF', 'FontSize', fontSize+2)

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
set(gcf,'units','normalized','outerposition',[0 0.15 1 0.75]); box on;
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');


%%%%%%%%%%%%%%%%%%%%%%%%%% alpha_s (ACS) %%%%%%%%%%%%%%%%%%%%%%%%%%
acs_sam = alpha_ratio + alpha_ref;

units           = 1E3;
bmodeFull       = bmode_sam;
colorImg        = acs_sam;
range_bmode     = [-60 0];
range_img       = [0.8 1.2];
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



%% BSC RECONSTRUCTION

freq = spectralData_sam.band;

g_est = mean(g_ratio(:));
s_est = mean(s_ratio(:));

bsc_est_gauss = g_est.*exp(s_est.*freq.^2);

% Fit BSC to powerLaw model

% HC WAY
% M = [ones(length(freq),1), log(freq)];
% r = log(bsc_est_gauss);
% 
% x = cgs(M'*M,M'*r,1e-16,10000);
% 
% b = exp(x(1));
% n = x(2);

% EMZ Way
% Transform to logarithmic scale lo linearize
% log_f   = log(freq);
% log_eta = log(bsc_est_gauss);

% Perform linear regression
coeffs  = polyfit(log(freq(:)), log(bsc_est_gauss), 1); % Fit y = mx + c
n       = coeffs(1); % Slope = n
ln_b    = coeffs(2); % Intercept = ln(b)
b       = exp(ln_b);

bsc_fit_powlaw = b*freq.^n;

% Display results
fprintf('-----PowerLaw-----\n')
fprintf('n = %f\n', n);
fprintf('b = %f\n', b);
fprintf('b = %-E\n', b);
fprintf('b_dB = %.2f\n', 10*log10(b))
fprintf('-----PowerLaw-----\n')


%
% figure, 
% sgtitle('\bfTarget 7')
% plot(freqs,  bsc, 'b', 'DisplayName', 'CIRS'), hold on
% grid on;
% plot(freqs,  b*f.^n, 'k--', 'DisplayName', 'Fit'), hold off
% xlabel('Freq [MHz]'), ylabel('BSC')
% 
% legend('Location','best')


% Theretical
g_the = 1;
s_the = -0.0245;
bsc_the = g_the.*exp(s_the.*freq.^2);

figure;
plot(freq,10*log10(bsc_est_gauss),'r.-', 'DisplayName','gaussian-est'); hold on;
plot(freq,10*log10(bsc_fit_powlaw),'g.-', 'DisplayName', 'powlaw-fit'); hold on;
% plot(freq,10*log10(bsc_the),'b.-', 'DisplayName', 'gaussian-teo'); grid on;
grid;
legend('Location', 'best');
title('BSC')
ylabel('BSC [dB]')
set(gca,'FontSize',fontSize);

figure;
loglog(freq,(bsc_est_gauss),'r.-', 'DisplayName','gaussian-est'); hold on;
loglog(freq,(bsc_fit_powlaw),'g.-', 'DisplayName', 'powlaw-fit'); hold on;
% plot(freq,(bsc_the),'b.-', 'DisplayName', 'gaussian-teo'); grid on;
xlim([3 9])
grid;
legend('Location', 'best');
title('BSC')
ylabel('BSC [cm^{-1}\cdot sr^{-1}]')
set(gca,'FontSize', fontSize);


%% HELP GAUSSIAN KERNELS

% sizeKernel  = [7 7];
% sigmaKernel = 5;
% 
% dens_old    = medium.density;
% 
% kernelREF   = fspecial('gaussian', sizeKernel, sigmaKernel); % Kernel1: Size 15x15, Sigma 2
% 
% figure, 
% subplot(1,2,1)
% imagesc(kernelREF), axis('image'), colorbar
% title(['Kernel Gauss w = ', num2str(sizeKernel), ', \sigma=', num2str(sigmaKernel)])
% 
% 
% subplot(1,2,2)
% plot(kernelREF(:, ceil(sizeKernel(1) / 2) ))
% title(['Kernel Gauss w = ', num2str(sizeKernel), ', \sigma=', num2str(sigmaKernel)])
% 
% dens_gauss  = imfilter(dens_old, kernelREF, 'same', 'replicate');
% 
% %
% clim = [900 1100];
% nBins = 100;
% 
% figure, 
% 
% subplot(2,2,1)
% imagesc(dens_old, clim), colorbar, axis("image");
% title('Density old')
% 
% subplot(2,2,2)
% imagesc(dens_gauss, clim), colorbar, axis("image");
% title(['Density Gauss w = ', num2str(sizeKernel), ', \sigma=', num2str(sigmaKernel)])
% 
% subplot(2,2,[3 4])
% h_old = histogram(dens_old(:), nBins, 'FaceColor', 'b'); % Blue color for dens_old
% hold on
% h_gauss = histogram(dens_gauss(:), nBins, 'FaceColor', 'r'); % Red color for dens_gauss
% hold off
% legend('Density old', ['Density Gauss \sigma=', num2str(sigmaKernel)])


