% Test of calc_powerSpectra function (simplified Chahuara code)
% USING DATA SIMU KWAVE GAUSS ATTEMPT TO CHANGE BSC CHECK KERNEL VARIABLES
% IN SIMULATION
% ATTENUATION alpha_power  1.1
% OPTION AVAILABLE 2DoF with priors and 3-DoF by changing estim_method
% FebW3 2025 EMZ
% Main difference: New simulations, other Kgrid, better visualization in
% Bmode (pas très obsucre comment avant AC)

%%
% clear all, 
% clc, 
warning('off');
% close all;

methodsRegu = true;
% methodsRegu = false;

addpath(genpath(pwd))

Np2dB       = 20*log10(exp(1));
dB2Np       = 1/Np2dB;
log2log10   = log10(exp(1)); % pow2db(exp(1))
range_bmode = [-60 0];

plotBmode   = false;
plotBSCdB   = true;
plotMaps    = false;

methods = {'3-DoF', '2-DoF-a', '2-DoF-g', '2-DoF-s'};
% methods = {'2-DoF-s'};

% First row for headers, second for data
bsc_results = cell(2, length(methods)); 
maps_results = cell(4, length(methods));

% Store headers
bsc_results(1, :) = {sprintf('3-DoF'), sprintf('2-DoF "a"'), sprintf('2-DoF "b"'),  sprintf('2-DoF "n"')};
maps_results(1, :) = {sprintf('3-DoF'), sprintf('2-DoF "a"'), sprintf('2-DoF "b"'),  sprintf('2-DoF "n"')};

caseDoF = 1; %1, 21, 22, 23: ok, error a, error g, error s
%% LOAD DATA
pathData = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\dataACS_kwave\simFeb2025';

% DATA NEW AMZ
alpha_sam = 0.5; 
j_sam = 1.1;

alpha_ref = 0.5;
j_ref = j_sam;

folderDataSam = sprintf('Gauss1_a%.1f_pow%.1f', alpha_sam, j_sam);
folderDataSam = strrep(folderDataSam, '.', 'p');

rf_sam_name = strcat('rfsam_', sprintf('%.3f', j_sam));
rf_sam_name = strrep(rf_sam_name, '.', 'p');
rf_sam_name = strcat(rf_sam_name, '.mat');
SAM = load(fullfile(pathData, folderDataSam, rf_sam_name));

SAM.alpha_power = j_sam;
SAM.acs = alpha_sam; % [dB/cm/MHz] 

folderDataRef = sprintf('Gauss1_a%.1f_pow%.1f', alpha_ref, j_sam);
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
switch caseDoF 

    case 1 % original
        delta_alpha_prior = alpha_sam - alpha_ref; % [dB/cm/MHz]
        delta_g_prior     = log(1);  % ln (g_sam / g_ref) 0dB
        delta_s_prior     = 0.09; % s_sam - s_ref og BSC
        % delta_s_prior     = 0.022; % s_sam - s_ref og ACS

    case 21 % 2-DoF, incorrect alpha prior
        delta_alpha_prior = -0.1; % [dB/cm/MHz]

        delta_g_prior     = log(1);  % ln (g_sam / g_ref) 0dB
        delta_s_prior     = 0.0245; % s_sam - s_ref og BSC

    case 22 % 2-DoF, incorrect g prior
        delta_alpha_prior = alpha_sam - alpha_ref; % [dB/cm/MHz]

        % delta_g_prior     = log(0.5);  % ln (g_sam / g_ref) eq log(db2pow(3dB))
        delta_g_prior     = log(db2pow(-5));
        delta_s_prior     = 0.0245; % s_sam - s_ref og BSC

    case 23 % 2-DoF, incorrect s prior
        delta_alpha_prior = alpha_sam - alpha_ref; % [dB/cm/MHz]
        delta_g_prior     = log(1);  % ln (g_sam / g_ref) 0dB

        delta_s_prior     = -0.05; % s_sam - s_ref og BSC

end

% switch j_sam 
% 
%     case 1.1
%         delta_alpha_prior = alpha_sam - alpha_ref; % [dB/cm/MHz]
%         delta_g_prior     = log(1);  % ln (g_sam / g_ref) 0dB
%         delta_s_prior     = 0.0245; % s_sam - s_ref og BSC
%         % delta_s_prior     = 0.022; % s_sam - s_ref og ACS
% 
% end

% B-MODE CHECK

bmode_sam = db(abs(hilbert(SAM.rf)));
bmode_sam = bmode_sam - max(bmode_sam(:));

bmode_ref = db(abs(hilbert(REF.rf)));
bmode_ref = bmode_ref - max(bmode_ref(:));

%% SPECTRAL METHOD PARAMETERS
pars.P = 4096; % NFFT only for calculate BSC_RPM_ok 10wl
% pars.P = 8192; % 15wl
pars.bw          = [3 8.5]; % [MHz] % old
pars.bw          = [3 9]; % [MHz] % I think better for BSC performance
pars.bw          = [3 8.75]; % [MHz] % new
pars.overlap     = 0.8;
pars.blocksize   = 12; % wavelengths
% pars.z_roi       = [5 25]*1E-3; % [m]  % half
% pars.z_roi       = [5 40]*1E-3; % all

% new
pars.z_roi       = [5 45]*1E-3; % all
pars.x_roi       = [-18 18]*1E-3;
pars.bw          = [4 8];

pars.window_type = 3; %  (1) Hanning, (2) Tuckey, (3) Hamming, (4) Tchebychev
pars.saran_layer = false;

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

%% POWER SPECTRA ESTIMATION
% spectralData_sam = calc_powerSpectra(SAM, pars);
spectralData_sam = calc_powerSpectra_vSimple(SAM, pars); % @
S_sam = spectralData_sam.powerSpectra;

num_ref = 1;

% spectralData_ref = calc_powerSpectra(REF, pars);
spectralData_ref = calc_powerSpectra_vSimple(REF, pars); % @
S_ref = spectralData_ref.powerSpectra;

%%
%% PLOT FULL Spectral
% Get size of Sp
[m, n, pFull] = size(spectralData_sam.Sfull);

nLines = 5;
lin_cen = round(n / 2); 
lat_range = max(1, lin_cen-fix(nLines/2)):min(n, lin_cen+fix(nLines/2)); 

Sd_matrix = squeeze(mean(spectralData_sam.Sfull(:, lat_range, :), 2)); % Mean over 2nd dim (lateral)

freq_Pos   = (0:(pFull/2-1))*SAM.fs/pFull*1e-6; % freq_rangePos = size ( (0:(1/2-1/pFull))'*SAM.fs) ;
Sd_Pos     = Sd_matrix(:,1:pFull/2);

Sd_Pos_dB  = pow2db(Sd_Pos ./ max(Sd_Pos, [], 2));
% Plot the FULL extracted spectra
figure;
set(gcf,'units','normalized','outerposition',[0 0.1 0.75 0.5]); box on;
subplot(1,3,1)
imagesc(freq_Pos, spectralData_sam.depth*1e3, Sd_Pos);
xlabel('Frequency [MHz]');
ylabel('Depth [mm]');
colorbar;
title('Power Spectrum');

subplot(1,3,2)
imagesc(freq_Pos, spectralData_sam.depth*1e3, Sd_Pos_dB);
xlabel('Frequency [MHz]');
ylabel('Depth [mm]');
h2 = colorbar; 
ylabel(h2,'dB');
title('Power Spectrum');

subplot(1,3,3)
plot(freq_Pos, Sd_Pos_dB(1, :), 'DisplayName', 'Top')
hold on, grid on;
plot(freq_Pos, Sd_Pos_dB(round(m/2), :), 'DisplayName', 'Half')
plot(freq_Pos, Sd_Pos_dB(end, :), 'DisplayName', 'Bottom')
yline(-20, 'k--')
xlim([0 15]);
hold off;
xlabel('Frequency [MHz]');
ylabel('NormMax [dB]');
title('Power Spectrum');
legend('Location', 'Best');

% FIND PEAKS at all depths

mDepths = size(Sd_Pos, 1);
f = freq_Pos;
ratio = db2mag(-20);
arrayBands = zeros(mDepths,  2);
for dd = 1:mDepths
    y = Sd_Pos(dd, :);
    [fLeft,fRight] = findFreqBand(f, y, ratio);
    arrayBands(dd,1) = fLeft;
    arrayBands(dd,2) = fRight;
end

figure,
imagesc(freq_Pos, spectralData_sam.depth*1e3, Sd_Pos_dB);
hold on
plot(arrayBands(:,1), spectralData_sam.depth * 1e3, 'r', 'LineWidth', 2); % Left boundary
plot(arrayBands(:,2), spectralData_sam.depth * 1e3, 'r', 'LineWidth', 2); % Right boundar
hold off
xlabel('Frequency [MHz]');
ylabel('Depth [mm]');
h2 = colorbar; 
ylabel(h2,'dB');
title('Norm Power Spectrum');

%%
SR_emz = S_sam ./ S_ref;

SR = permute(SR_emz,[3,1,2]);

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
% mu_rpl_tv    = [1E3; 1E3; 1E4]; % [mu_g, mu_s, mu_a]
% mu_rpl_tv    = [0.01; 0.01; 0.01]; % [mu_g, mu_s, mu_a]
% mu_rpl_tv    = [0.001; 0.001; 0.001]; % [mu_g, mu_s, mu_a]

%% FOR BUCLE
for iMet = 1:length(methods)

estim_method = methods{iMet};

%% COMPENSATE GAUSS ATTEMPT 2-DoF-a
if strcmp( estim_method, '2-DoF-a')

    if (methodsRegu); mu_rpl_tv    = [1E3; 1E3; 1E4]; % [mu_b, mu_n, mu_a]
    else              mu_rpl_tv    = [0.001 0.001 0.001];
    end

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

    if (methodsRegu); mu_rpl_tv    = [1E3; 1E3; 1E4]; % [mu_b, mu_n, mu_a]
    else              mu_rpl_tv    = [0.001 0.001 0.001];
    end

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

    if (methodsRegu); mu_rpl_tv    = [1E3; 1E3; 1E4]; % [mu_b, mu_n, mu_a]
    else              mu_rpl_tv    = [0.001 0.001 0.001];
    end

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

    if (methodsRegu); mu_rpl_tv    = [1E3; 1E3; 1E4]; % [mu_b, mu_n, mu_a]
    else              mu_rpl_tv    = [0.001 0.001 0.001];
    end

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
fprintf('b_s/b_r[dB]: %.3f ± %.4f, %%CV = %.4f\n', round(m_g, 3), round(s_g, 4), round(cv_g, 4));
    else
fprintf('b_s/b_r    : %.3f ± %.4f, %%CV = %.4f\n', round(m_g, 3), round(s_g, 4), round(cv_g, 4));
    end

fprintf('Δn         : %.4f ± %.4f, %%CV = %.4f\n', round(m_s, 4), round(s_s, 4), round(cv_s, 4));
fprintf('--------\n');
 
%% IMAGESC PLOTS
Xaxis = spectralData_ref.lateral;
Zaxis = spectralData_ref.depth;
cm = 1e2;

axis_n = [0 1.2];
axis_a = [0 3];
axis_b = [-60 0]; % dB
fontSize = 16;

if plotMaps    
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

end
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

g_est = median(g_ratio(:));
s_est = median(s_ratio(:));

bsc_est_gauss = g_est.*exp(-s_est.*freq.^2);

% Power Law Comparison: (bsc_est_gauss = b. f^n)
% Linearize: log (bsc_est_Gauss = n.log(f) + log(b) 

% Option 1: Polyfit
coeffs_pl   = polyfit(log(freq), log(bsc_est_gauss), 1); % Fit y = m.x + c
d_n_pl      = coeffs_pl(1); % Slope = d_n  (mean),  (median)
ln_d_b_pl   = coeffs_pl(2); % Intercept = ln(d_b) 
d_b_pl      = exp(ln_d_b_pl); % 1.0917 (mean), 0.9079(median)

% Option 2: Matrix [log(BSC)] = [ log(f) | 1] * [ n ; log(b) ] rr = MM * qq
% MM           = [ log(freq), ones( size(freq) ) ]; 
% rr           = log(bsc_est_gauss);
% qq           = cgs(MM' * MM, MM' * rr, 1e-16, 100); % qq = MM \ rr;
% d_n_pl       = qq(1);
% d_b_pl       = exp(qq(2));


bsc_fit_powlaw = d_b_pl*freq.^d_n_pl;

% % Display results
% fprintf('-----PowerLaw (b.f^n)-----\n')
% fprintf('Δn           = %f\n', d_n_pl);
% fprintf('d_b          = %f\n', d_b_pl);
% fprintf('b_s/b_r [dB] = %f\n', 10*log10(d_b_pl));
% fprintf('---------\n')

%% SAVE ALL BSC EST GAUSS

bsc_results{2, iMet} = bsc_est_gauss;
bsc_results{3, iMet} = bsc_fit_powlaw;

%% SAVE ALL MAPS

maps_results{2, iMet} = acs_sam; 
maps_results{3, iMet} = g_ratio_dB; 
maps_results{4, iMet} = s_ratio; 


end
%%
keyboard

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BSC RPM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NAME BSC
nameBSC = strcat('BSC_', sprintf('sam%.2f_ref%.2f', alpha_sam, alpha_ref ));
nameBSC = strrep(nameBSC, '.', 'p');
% BSC = calculateBSC_RPM_ok(SAM, REF, pars); % slow only once**
BSC = calculateBSC_RPM_fast(SAM, REF, pars); % fast TBD**
% save(fullfile(pathData, folderDataSam, nameBSC), "BSC")
% load(fullfile(pathData, folderDataSam, nameBSC));
%%freq = BSC.band;
% bsc_rpm = BSC.BSCcurve_Uni(:,1); % mean
bsc_rpm = BSC.BSCcurve_Uni(:,2); % median

% MM = [freq.^2, ones(size(freq))];
% rr = log(bsc_rpm);
% pp = cgs(MM'*MM,MM'*rr,1e-16,10000); % coeffs

freq = BSC.band;

% Perform linear regression  bsc = d_s . f^2 + ln(d_g) 
coeffs   = polyfit(-freq.^2, log(bsc_rpm), 1); % Fit y = mx + c
d_s      = coeffs(1); % Slope = d_s -0.0319 (mean), -0.0317 (median)
ln_dg    = coeffs(2); % Intercept = ln(d_g) 
d_g      = exp(ln_dg); % 1.0917 (mean), 0.9079(median)

% Display results
fprintf('-----RPM Gauss (g.exp(-s.f^2))-----\n')
fprintf('d_g          = %f\n', d_g);
fprintf('g_s/g_r [dB] = %f\n', 10*log10(d_g));
fprintf('Δs           = %f\n', d_s);
fprintf('---------\n')

bsc_rpm_gauss = d_g*exp(-d_s* freq.^2);
%% Box Plot distribution DELTAs **
% Define method labels

delta_b_theo = d_g;
delta_n_theo = d_s;
numMethods = size(maps_results, 2); % Number of methods (iMet values)

% Extract Data
acs_data = cell(1, numMethods);
g_ratio_data = cell(1, numMethods);
s_ratio_data = cell(1, numMethods);

for iMet = 1:numMethods
    acs_data{iMet}     = maps_results{2, iMet}(:);  % Flatten to column
    g_ratio_data{iMet} = maps_results{3, iMet}(:);
    s_ratio_data{iMet} = maps_results{4, iMet}(:);
end

% Convert to matrix for plotting (ensuring correct format)
acs_mat     = padconcatenation(acs_data, NaN, 1); % Pad with NaN for different lengths
g_ratio_mat = padconcatenation(g_ratio_data, NaN, 1);
s_ratio_mat = padconcatenation(s_ratio_data, NaN, 1);

font_size = 18;
method_labels = string({maps_results{1, :}}); % Convert first row to string array

% Exclude the second column in plot a
acs_mat_filtered = acs_mat(:, [1, 3, 4]);
method_labels_a = method_labels([1, 3, 4]);

% Box Plot a
figure;
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 700, 700]); % [x, y, width, height] in pixels
box on;
boxplot(acs_mat_filtered, 'Labels', method_labels_a);
% axis("image")
yline(alpha_sam, 'k--')

if (methodsRegu); ylim([0.05 1.01])
else              ylim([-80 80])
end
title('\alpha');
ylabel('\alpha [dB\cdotcm^{-1}\cdotMHz^{-1}]');
set(gca, 'FontSize', font_size);

% Exclude the third column in plot b
g_ratio_mat_filtered = g_ratio_mat(:, [1, 2, 4]);
method_labels_b = method_labels([1, 2, 4]);

% Box Plot b
figure;
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 700, 700]); % [x, y, width, height] in pixels
box on;
boxplot(g_ratio_mat_filtered, 'Labels', method_labels_b);
% axis("image")
yline(10*log10(delta_b_theo), 'k--')
yline(0, 'k--')
if (methodsRegu); %ylim([-0.5 0.5])
else              ylim([-60 20])
end

title('\Deltab');
ylabel('\Deltab [dB]');
set(gca, 'FontSize', font_size);

% Exclude the fourth column in plot n
s_ratio_mat_filtered = s_ratio_mat(:, [1, 2, 3]);
method_labels_n = method_labels([1, 2, 3]);

% Box Plot n
figure;
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 700, 700]); % [x, y, width, height] in pixels
box on;
boxplot(s_ratio_mat_filtered, 'Labels', method_labels_n);
% axis("image")
yline(delta_n_theo, 'k--')
% if (methodsRegu); ylim([0 0.05]); yticks([0:0.01:0.05])
% else              ylim([-15 15])
% end
title('\Deltan');
ylabel('\Deltan [a.u.]');
set(gca, 'FontSize', font_size);

function M = padconcatenation(C, padval, dim)
    % C: Cell array to concatenate
    % padval: Value used to pad (e.g., NaN)
    % dim: Dimension along which to concatenate (1 = rows, 2 = columns)
    max_length = max(cellfun(@numel, C));
    M = cellfun(@(x) padarray(x, [max_length - numel(x), 0], padval, 'post'), C, 'UniformOutput', false);
    M = cell2mat(M);
end

%%
keyboard;



%% PLOTS BSC TOGETHER
xlim_range = [3 8.51];
ylim_range = [0.1 1.1]; % Y-axis limits

% Define properties for customization
line_width = 2.85; % Set line width
font_size = 16; % Adjust font size

% BSC THEORETICAL
bsc_delta_theo = bsc_rpm;

m_3dof.mpe   = mean( (bsc_results{2, 1} - bsc_delta_theo)./bsc_delta_theo );
m_3dof.mae   = mean( abs(bsc_results{2, 1} - bsc_delta_theo)./bsc_delta_theo );
m_3dof.rmse  = rmse(bsc_results{2, 1}, bsc_delta_theo);
m_3dof.nrmse = m_3dof.rmse / mean(bsc_results{2, 1});

m_2dofa.mpe   = mean( (bsc_results{2, 2} - bsc_delta_theo)./bsc_delta_theo );
m_2dofa.mae   = mean( abs(bsc_results{2, 2} - bsc_delta_theo)./bsc_delta_theo );
m_2dofa.rmse  = rmse(bsc_results{2, 2}, bsc_delta_theo);
m_2dofa.nrmse = m_2dofa.rmse / mean(bsc_results{2, 2});

m_2dofb.mpe   = mean( (bsc_results{2, 3} - bsc_delta_theo)./bsc_delta_theo );
m_2dofb.mae   = mean( abs(bsc_results{2, 3} - bsc_delta_theo)./bsc_delta_theo );
m_2dofb.rmse  = rmse(bsc_results{2, 3}, bsc_delta_theo);
m_2dofb.nrmse = m_2dofb.rmse / mean(bsc_results{2, 3});

m_2dofn.mpe   = mean( (bsc_results{2, 4} - bsc_delta_theo)./bsc_delta_theo );
m_2dofn.mae   = mean( abs(bsc_results{2, 4} - bsc_delta_theo)./bsc_delta_theo );
m_2dofn.rmse  = rmse(bsc_results{2, 4}, bsc_delta_theo);
m_2dofn.nrmse = m_2dofn.rmse / mean(bsc_results{2, 4});

% Extract field names
fields = fieldnames(m_3dof);

% Create a table
Tbsc = table(struct2cell(m_3dof), struct2cell(m_2dofa), ...
    struct2cell(m_2dofb), struct2cell(m_2dofn), 'RowNames', fields, 'VariableNames', methods);


bsc_results{1, 1} = sprintf('3-DoF      (RMSE = %.2f%%) \n', 100*m_3dof.rmse);
bsc_results{1, 2} = sprintf('2-DoF "a" (RMSE = %.2f%%) \n', 100*m_2dofa.rmse);
bsc_results{1, 3} = sprintf('2-DoF "b" (RMSE = %.2f%%) \n', 100*m_2dofb.rmse);
bsc_results{1, 4} = sprintf('2-DoF "n" (RMSE = %.2f%%) \n', 100*m_2dofn.rmse);

% Convert hexadecimal colors to RGB (MATLAB requires values between 0 and 1)
colot_gt = '#000000'; % Black
color_1  = '#FF0000';  % 3dof
color_2  = '#D95319';  % 2dof a
color_3  = '#0072BD';  % 2dof b
color_4  = '#77AC30';  % 2dof n

% Create figure and plot data
figure, 
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 700, 700]); % [x, y, width, height] in pixels

semilogy(freq, bsc_delta_theo, '-', 'Color', hex2rgb(colot_gt), 'LineWidth', line_width+0.5, 'DisplayName', 'GT');
hold on;
semilogy(freq, bsc_results{2, 1}, '--', 'Color', hex2rgb(color_1), 'LineWidth', line_width, 'DisplayName', bsc_results{1, 1});
semilogy(freq, bsc_results{2, 2}, '--', 'Color', hex2rgb(color_2), 'LineWidth', line_width, 'DisplayName', bsc_results{1, 2});
semilogy(freq, bsc_results{2, 3}, '--', 'Color', hex2rgb(color_3), 'LineWidth', line_width, 'DisplayName', bsc_results{1, 3});
semilogy(freq, bsc_results{2, 4}, '--', 'Color', hex2rgb(color_4), 'LineWidth', line_width, 'DisplayName', bsc_results{1, 4});
hold off;

% Customize plot
grid on;
xlabel('Frequency [MHz]', 'FontSize', font_size);
ylabel('BSC [cm^{-1}\cdot sr^{-1}]', 'FontSize', font_size);
% ylim(ylim_range);
xlim(xlim_range);
title('Backscatter Coefficient (BSC)', 'FontSize', font_size + 2);
legend('Location', 'southwest', 'FontSize', font_size + 2);
set(gca, 'FontSize', font_size);
hold off;


%%
% keyboard
% %% PLOTS ACS MAPS TOGETHER
% 
% acs_results{1, 1} = '3-DoF';
% acs_results{1, 2} = '2-DoF "a"';
% acs_results{1, 3} = '2-DoF "b"';
% acs_results{1, 4} = '2-DoF "n" ';
% 
% fontSize = 16;
% 
% figure,
% % set(gcf,'units','normalized','outerposition',[0 0.1 1 0.8]); box on;
% 
% tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
% sgtitle(estim_method, 'FontSize', fontSize+2, 'FontWeight', 'bold');
% 
% for iMet = 1:length(methods)
% %%%%%%%%%%%%%%%%%%%%%%%%%% alpha_s (ACS) %%%%%%%%%%%%%%%%%%%%%%%%%%
% acs_sam = alpha_ratio + alpha_ref;
% 
% units           = 1E3;
% bmodeFull       = bmode_sam;
% colorImg        = acs_results{2, iMet};
% range_bmode     = [-100 0];
% range_img       = [0.1 1.2];
% transparency    = 0.65;
% x_img           = spectralData_sam.lateral*units;
% z_img           = spectralData_sam.depth*units;
% xFull           = SAM.x*units;
% zFull           = SAM.z*units;
% [X, Z] = meshgrid(xFull, zFull);
% z_img(end) = 45;
% roi = and(X >= x_img(1), X <= x_img(end)) & ...
%       and(Z >= z_img(1), Z <= z_img(end));
% 
% t = nexttile;
%     [~,hB,hColor] = imOverlayInterp(bmodeFull, colorImg, range_bmode, range_img, ...
%                         transparency, x_img, z_img, roi, xFull, zFull);   
%     hold on;
%     contour(xFull, zFull, roi, 1,'w--', 'LineWidth', 2)
%     hold off;
%     xlabel('Lateral [cm]'), ylabel('Depth [cm]');
%     hColor.Label.String = 'dB\cdotcm^{-1}\cdotMHz^{-1}';
%     title(bsc_results{1, iMet})
%     set(gca,'fontsize',fontSize)
% 
% end
%%
