% Test of calc_powerSpectra function (simplified Chahuara code)
% USING DATA SIMU KWAVE POWER LAW COLOURATION
% THIS TIME 5 MEDIUMS GENERATED TO SEE AVERAGE ESTIMATIONS IN ALL CASES
% a = 0.6, b = 0.01, n = 1.5
% a = 0.4, b = 1, n = 0
% files are tipically name "rfa0.4b0.01n1.5.mat"
% FebW3 TUFFC 25

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
range_bmode = [-60 0];

plotBmode   = false;
plotBSCdB   = true;
plotMaps    = false;


methods = {'3-DoF', '2-DoF-a', '2-DoF-b', '2-DoF-n'};

% methods = {sprintf('2-DoF "a"'), sprintf('2-DoF "b"'),  sprintf('2-DoF "n"')};
% methods = {'3-DoF'};
% methods = {'2-DoF-b'};

% First row for headers, second for data
bsc_results = cell(3, length(methods)); 
maps_results = cell(2, length(methods));

% Store headers
bsc_results(1, :) = {sprintf('3-DoF'), sprintf('2-DoF "a"'), sprintf('2-DoF "b"'),  sprintf('2-DoF "n"')};
maps_results(1, :) = {sprintf('3-DoF'), sprintf('2-DoF "a"'), sprintf('2-DoF "b"'),  sprintf('2-DoF "n"')};
%% LOAD DATA
pathData = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\dataACS_kwave\simFeb2025\PowerLaw1p';
folderDataSam = 'a0p6_pow1p';
folderDataRef = 'a0p4_pow1p';

% DATA NEW AMZ
alpha_sam = 0.6;
b_sam = 0.01;
n_sam = 1.5;
j_sam = 1.0;

alpha_ref = 0.4;
b_ref = 1;
n_ref = 0;
j_ref = j_sam;

%% ONE CASE
% rf_sam_name = ['rf5_a',num2str(alpha_sam),'b',num2str(b_sam),'n',num2str(n_sam)];
% SAM = load(fullfile(pathData, folderDataSam, rf_sam_name));
% SAM.alpha_power = j_sam;
% SAM.acs = alpha_sam; % [dB/cm/MHz] 
% 
% rf_ref_name = ['rf5_a',num2str(alpha_ref),'b',num2str(b_ref),'n',num2str(n_ref)];
% REF = load(fullfile(pathData, folderDataRef, rf_ref_name)); 
% REF.alpha_power = j_ref; 
% REF.acs = alpha_ref; % [dB/cm/MHz]
%% ALL CASES AVERAGE
%%%%%%%%%%%%%%%%%%%%%%%% SAM %%%%%%%%%%%%%%%%%%%%%%%%
rfIndices = [1,2,3,4,5];  % Example: Load rf1, rf3, and rf5
numRF = length(rfIndices);  % Number of RF files selected

fileName = ['rf',num2str(rfIndices(1)),'_a',num2str(alpha_sam),'b',num2str(b_sam),'n',num2str(n_sam)];

% Initialize SAM structure
SAM = load(fullfile(pathData, folderDataSam, fileName));
rfSize = size(SAM.rf);  % Get dimensions from the first file

% Preallocate the RF 3D array
newrf = nan([rfSize, numRF], 'like', SAM.rf);  % Use 'like' for consistency

% Load and store RF data into the 3D array
for i = 1:numRF

    fileName = ['rf',num2str(rfIndices(i)),'_a',num2str(alpha_sam),'b',num2str(b_sam),'n',num2str(n_sam)];
    filePath = fullfile(pathData, folderDataSam, fileName);
    
    fprintf("Loading SAM: %s\n", filePath);
    
    newrf(:,:,i) = load(filePath, 'rf').rf;  % Load only 'rf' variable
end

% Store in the structure
SAM.acs = alpha_sam;
SAM.rf = newrf;
%%%%%%%%%%%%%%%%%%%%%%%% SAM %%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% REF %%%%%%%%%%%%%%%%%%%%%%%%
% rfIndices = [1,2, 3, 4,5];  % Example: Load rf1, rf3, and rf5
numRF = length(rfIndices);  % Number of RF files selected

fileName = ['rf',num2str(rfIndices(1)),'_a',num2str(alpha_sam),'b',num2str(b_sam),'n',num2str(n_sam)];

% Initialize SAM structure
REF = load(fullfile(pathData, folderDataSam, fileName));
rfSize = size(REF.rf);  % Get dimensions from the first file

% Preallocate the RF 3D array
newrf = nan([rfSize, numRF], 'like', REF.rf);  % Use 'like' for consistency

% Load and store RF data into the 3D array
for i = 1:numRF

    fileName = ['rf',num2str(rfIndices(i)),'_a',num2str(alpha_ref),'b',num2str(b_ref),'n',num2str(n_ref)];
    filePath = fullfile(pathData, folderDataRef, fileName);
    
    fprintf("Loading REF: %s\n", filePath);
    
    newrf(:,:,i) = load(filePath, 'rf').rf;  % Load only 'rf' variable
end

% Store in the structure
REF.acs = alpha_ref;
REF.rf = newrf;
%%%%%%%%%%%%%%%%%%%%%%%% REF %%%%%%%%%%%%%%%%%%%%%%%%

%%

% FROM RPM (GROUNTRUTH) ESTIMATE delta_n_prior
switch j_sam 
    case 1.0
        delta_alpha_prior = alpha_sam - alpha_ref; % [dB/cm/MHz]
        delta_b_prior     = log(b_sam / b_ref);  % ln (g_sam / g_ref) ln(0.01/1)
        delta_n_prior     = n_sam - n_ref; % s_sam - s_ref og ACS
end

% B-MODE CHECK

bmode_sam = db(abs(hilbert(SAM.rf(:,:,1))));
bmode_sam = bmode_sam - max(bmode_sam(:));

bmode_ref = db(abs(hilbert(REF.rf(:,:,1))));
bmode_ref = bmode_ref - max(bmode_ref(:));

%% SPECTRAL METHOD PARAMETERS
pars.P = 2048; % NFFT only for calculate BSC_delta_b_priorRPM_ok 10wl
pars.bw          = [3 9]; % [MHz]
pars.overlap     = 0.8;
pars.blocksize   = 10; % wavelengths

% % try to get before PLN
% pars.blocksize   = 10; % wavelengths
% pars.z_roi       = [10 45]*1E-3; 
% pars.x_roi       = [-15 15]*1E-3; 
% pars.bw          = [3 9];

% % TUFFC
% pars.P = 1024;
pars.blocksize   = 12; % wavelengths
pars.z_roi       = [10 45]*1E-3; 
pars.x_roi       = [-15 15]*1E-3; 
pars.bw          = [3 8.75];

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
%% POWER SPECTRA ESTIMATION
% spectralData_sam = calc_powerSpectra(SAM, pars);
spectralData_sam = calc_powerSpectra_vSimple(SAM, pars); % @
S_sam = spectralData_sam.powerSpectra;

num_ref = 1;

% spectralData_ref = calc_powerSpectra(REF, pars);
spectralData_ref = calc_powerSpectra_vSimple(REF, pars); % @
S_ref = spectralData_ref.powerSpectra;

% Ratio Computation
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

%% GENERAL REGULARIZATION SETTINGS
% Implementation parameters
par_rpl.tol        = 1e-16;
par_rpl.kmax       = 100;
par_rpl.eps_f      = 1e-16;
par_rpl.m_est      = 0; %Robust
% par_rpl.ini_tol    = 1e-16;
par_rpl.ini_tol    = 1e-5;
% par_rpl.df_op      = 1;
par_rpl.df_op      = 0;
par_rpl.ini_method = 1; % METHOD LEAST SQUARES INITIALIZATION 

%% FOR BUCLE
for iMet = 1:length(methods)

estim_method = methods{iMet};

%% COMPENSATE GAUSS ATTEMPT 2-DoF-a
if strcmp( estim_method, '2-DoF-a')

    if (methodsRegu); mu_rpl_tv    = [1E3; 10^3.5; 1E4]; % [mu_b, mu_n, mu_a]
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
Z = kron( speye(p*q), log(f) );

% initialization for RPL-based methods
u_0 = initialize_rpl_a_prior(Y, X, Z, mu_rpl_tv, par_rpl);
    
% RPL estimation
[u_opt,~] = rpl_tv_a_prior(Y, X, Z, mu_rpl_tv, u_0, par_rpl);

if par_rpl.df_op == 1
dy = 0.5*(diag(ones(p-1,1),1) - diag(ones(p-1,1),-1));
dy(1,1) = -1; dy(1,2) = 1; dy(end,end) = 1; dy(end,end-1) = -1;
Dy = sparse(kron(speye(q),dy));
else
dy = diag(ones(p-1,1),1) - diag([ones(p-1,1);0]);
Dy = sparse(kron(speye(q),dy)); %diag(ones(p*q -1,1),1) - diag(ones(p*q,1));    
end

% \Deltas
g = u_opt(1:p*q);
s = u_opt(p*q+1:2*p*q);

% Prior "a" known
a_Np2dB = delta_alpha_prior*ones(p*q, 1);

% utils 
z = 1E2*repmat(depth,1,q); % 1E2*spectralData_sam.depth * ones(1, q); % 2d array
dz = reshape(Dy*z(:),p,q);
dz(end,:) = dz(end-1,:);

%% COMPENSATE GAUSS ATTEMPT 2-DoF-n
elseif strcmp( estim_method, '2-DoF-n')

    if (methodsRegu); mu_rpl_tv    = [1E3; 1E3; 1E4]; % [mu_b, mu_n, mu_a]
    else              mu_rpl_tv    = [0.001 0.001 0.001];
    end


band    = spectralData_sam.band;
depth   = spectralData_sam.depth;
[r,p,q] = size(SR);
  
comp_ref    = comp_ref_n_bsc(delta_n_prior, band, p, q);
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
s = delta_n_prior*ones(p*q, 1);

% utils 
z = 1E2*repmat(depth,1,q); % 1E2*spectralData_sam.depth * ones(1, q); % 2d array
dz = reshape(Dy*z(:),p,q);
dz(end,:) = dz(end-1,:);  

a_Np2dB = Np2dB*Dy*a./dz(:);

%% COMPENSATE GAUSS ATTEMPT 2-DoF-b
elseif strcmp( estim_method, '2-DoF-b')

    if (methodsRegu); mu_rpl_tv    = [1E3; 1E3; 1E4]; % [mu_b, mu_n, mu_a]
    else              mu_rpl_tv    = [0.001 0.001 0.001];
    end

band    = spectralData_sam.band;
depth   = spectralData_sam.depth;
[r,p,q] = size(SR);
  
comp_ref    = comp_ref_b_bsc(delta_b_prior);
comp_freq_a = comp_mod_freq_a(alpha_ref,j_sam,j_ref,band,depth,q);

comp_ref = 102;
SR_comp = SR .* comp_ref .*comp_freq_a;

% indices initialization
f = band(:); % always column vector
% [r,p,q] = size(SR_comp);

% log-spectrum Ratio Y_g = Z.s + W.a 
Y = log(SR_comp);

% matrices for RPL-based algorithms
Z = kron( speye(p*q), log(f) ); % EMZ PowLaw  Size: [p*q*r, p*q] 
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
g = delta_b_prior*ones(p*q, 1);

% utils 
z = 1E2*repmat(depth,1,q); % 1E2*spectralData_sam.depth * ones(1, q); % 2d array
dz = reshape(Dy*z(:),p,q);
dz(end,:) = dz(end-1,:);       

a_Np2dB = Np2dB*Dy*a./dz(:);

%% COMPENSATE GAUSS ATTEMPT 3-DoF
elseif strcmp( estim_method, '3-DoF')

    if (methodsRegu); mu_rpl_tv    = [10^2; 10^2.5; 10^3.885]; % [mu_b, mu_n, mu_a]
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
Z = kron( speye(p*q), log(f) ); % EMZ PowLaw  Size: [p*q*r, p*q] 
% Z = kron( speye(p*q), -f.^2 ); % EMZ Gauss  Size: [p*q*r, p*q] 
% W = kron( speye(p*q), -4*f.^j_sam );
W = kron( speye(p*q), -4*f );

% initialization for RPL-based methods
u_0 = initialize_rpl(Y, X, Z, W, mu_rpl_tv, par_rpl);
    
% RPL estimation
[u_opt,~] = rpl_tv(Y, X, Z, W, mu_rpl_tv, u_0, par_rpl);

if par_rpl.df_op == 1
dy = 0.5*(diag(ones(p-1,1),1) - diag(ones(p-1,1),-1));
dy(1,1) = -1; dy(1,2) = 1; dy(end,end) = 1; dy(end,end-1) = -1;
Dy = sparse(kron(speye(q),dy));
else
dy = diag(ones(p-1,1),1) - diag([ones(p-1,1);0]);
Dy = sparse(kron(speye(q),dy)); %diag(ones(p*q -1,1),1) - diag(ones(p*q,1));    
end

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

% Median
g_est = median(g_ratio(:));
s_est = median(s_ratio(:));
% Reconstruct BSC POW LAW
bsc_est_powlaw = g_est.*(freq.^s_est);
% SAVE ALL BSC EST GAUSS (median)
bsc_results{2, iMet} = bsc_est_powlaw;

% Average
g_est = mean(g_ratio(:));
s_est = mean(s_ratio(:));
% Reconstruct BSC POW LAW
bsc_est_powlaw = g_est.*(freq.^s_est);
% SAVE ALL BSC EST GAUSS (average)
bsc_results{3, iMet} = bsc_est_powlaw;



%% SAVE ALL MAPS

maps_results{2, iMet} = acs_sam; 
maps_results{3, iMet} = g_ratio_dB; 
maps_results{4, iMet} = s_ratio; 

end
%%
% keyboard

%% NEW DISTRIBUTION DELTAs**
% Define method labels

delta_b_theo = b_sam / b_ref;
delta_n_theo = n_sam - n_ref;
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

if (methodsRegu); ylim([0.25 0.95])
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
if (methodsRegu); ylim([-22 -19])
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
if (methodsRegu); ylim([1.36 1.79])
else              ylim([-15 15])
end
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
%% METRICS TABLE FORM (ACS)

m_3dof  = get_metrics_homo_gt(maps_results{2, 1}, logical(ones(size(acs_sam))), alpha_sam, '3-DoF');
m_2dofa = get_metrics_homo_gt(maps_results{2, 2}, logical(ones(size(acs_sam))), alpha_sam, '2-DoF-a');
m_2dofb = get_metrics_homo_gt(maps_results{2, 3}, logical(ones(size(acs_sam))), alpha_sam, '2-DoF-b');
m_2dofn = get_metrics_homo_gt(maps_results{2, 4}, logical(ones(size(acs_sam))), alpha_sam, '2-DoF-n');

% Extract field names
fields = fieldnames(m_3dof);

% Create a table
Ta = table(struct2cell(m_3dof), struct2cell(m_2dofa), ...
    struct2cell(m_2dofb), struct2cell(m_2dofn), 'RowNames', fields, 'VariableNames', methods);

m_3dof  = get_metrics_homo_gt(maps_results{3, 1}, logical(ones(size(acs_sam))), alpha_sam, '3-DoF');
m_2dofa = get_metrics_homo_gt(maps_results{3, 2}, logical(ones(size(acs_sam))), alpha_sam, '2-DoF-a');
m_2dofb = get_metrics_homo_gt(maps_results{3, 3}, logical(ones(size(acs_sam))), alpha_sam, '2-DoF-b');
m_2dofn = get_metrics_homo_gt(maps_results{3, 4}, logical(ones(size(acs_sam))), alpha_sam, '2-DoF-n');

Tb = table(struct2cell(m_3dof), struct2cell(m_2dofa), ...
    struct2cell(m_2dofb), struct2cell(m_2dofn), 'RowNames', fields, 'VariableNames', methods);

m_3dof  = get_metrics_homo_gt(maps_results{4, 1}, logical(ones(size(acs_sam))), alpha_sam, '3-DoF');
m_2dofa = get_metrics_homo_gt(maps_results{4, 2}, logical(ones(size(acs_sam))), alpha_sam, '2-DoF-a');
m_2dofb = get_metrics_homo_gt(maps_results{4, 3}, logical(ones(size(acs_sam))), alpha_sam, '2-DoF-b');
m_2dofn = get_metrics_homo_gt(maps_results{4, 4}, logical(ones(size(acs_sam))), alpha_sam, '2-DoF-n');

Tn = table(struct2cell(m_3dof), struct2cell(m_2dofa), ...
    struct2cell(m_2dofb), struct2cell(m_2dofn), 'RowNames', fields, 'VariableNames', methods);

clear m_3dof m_2dofa m_2dofb m_2dofn
% T = [Ta; Tb; Tn]; clear Ta Tb Tn m_3dof m_2dofa m_2dofb m_2dofn

%% PLOTS BSC TOGETHER
xlim_range = [3 8.71];
ylim_range = [0.03 1]; % Y-axis limits

% Define properties for customization
line_width = 2.85; % Set line width
font_size = 16; % Adjust font size

% BSC THEORETICAL
delta_b_theo = b_sam / b_ref;
delta_n_theo = n_sam - n_ref;
bsc_delta_theo = delta_b_theo* (freq.^delta_n_theo);

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
ylim(ylim_range);
xlim(xlim_range);
title('Backscatter Coefficient (BSC)', 'FontSize', font_size + 2);
legend('Location', 'northwest', 'FontSize', font_size + 2);
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

%%
keyboard

%% PLOTS ACS MAPS TOGETHER
% EXCEP PRIOR a
fontSize = 16;

figure,
% set(gcf,'units','normalized','outerposition',[0 0.1 1 0.8]); box on;

% tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact'); *old
 tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
sgtitle('\alpha', 'FontSize', fontSize+2, 'FontWeight', 'bold');

idx_m = [1, 3, 4];
methods_alpha = methods(idx_m);
% for iMet = 1:length(methods)*old
for jjj = 1:length(methods_alpha)
iMet = idx_m(jjj);
%%%%%%%%%%%%%%%%%%%%%%%%%% alpha_s (ACS) %%%%%%%%%%%%%%%%%%%%%%%%%%
acs_sam = alpha_ratio + alpha_ref;

units           = 1E3;
bmodeFull       = bmode_sam;
colorImg        = maps_results{2, iMet};
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
    xlabel('Lateral [mm]'), ylabel('Depth [mm]');
    hColor.Label.String = 'dB\cdotcm^{-1}\cdotMHz^{-1}';
    title(maps_results{1, iMet})
    set(gca,'fontsize',fontSize)

end
%%

%% FIGURE DISTRIBUTION 

% delta_b_theo = b_sam / b_ref;
% delta_n_theo = n_sam - n_ref;
% numMethods = size(maps_results, 2); % Number of methods (iMet values)
% 
% % Extract Data
% acs_data = cell(1, numMethods);
% g_ratio_data = cell(1, numMethods);
% s_ratio_data = cell(1, numMethods);
% 
% for iMet = 1:numMethods
%     acs_data{iMet}     = maps_results{2, iMet}(:);  % Flatten to column
%     g_ratio_data{iMet} = maps_results{3, iMet}(:);
%     s_ratio_data{iMet} = maps_results{4, iMet}(:);
% end
% 
% % Convert to matrix for plotting (ensuring correct format)
% acs_mat     = padconcatenation(acs_data, NaN, 1); % Pad with NaN for different lengths
% g_ratio_mat = padconcatenation(g_ratio_data, NaN, 1);
% s_ratio_mat = padconcatenation(s_ratio_data, NaN, 1);
% 
% % Define method labels
% method_labels = string({maps_results{1, :}}); % Convert first row to string array
% 
% % Box Plot a
% figure;
% boxplot(acs_mat, 'Labels', method_labels);
% yline(alpha_sam, 'k--')
% ylim([0.25 0.95])
% title('ACS');
% ylabel('ACS [dB\cdotcm^{-1}\cdotMHz^{-1}]');
% set(gca, 'FontSize', font_size);
% 
% % Box Plot b
% figure;
% boxplot(g_ratio_mat, 'Labels', method_labels);
% yline(10*log10(delta_b_theo), 'k--')
% ylim([-22 -19])
% title('\Deltab');
% ylabel('\Deltab [dB]');
% set(gca, 'FontSize', font_size);
% 
% % Box Plot n
% figure;
% boxplot(s_ratio_mat, 'Labels', method_labels);
% yline(delta_n_theo, 'k--')
% ylim([1.36 1.79])
% title('\Deltan');
% ylabel('\Deltan [a.u.]');
% set(gca, 'FontSize', font_size);
% 
% function M = padconcatenation(C, padval, dim)
%     % C: Cell array to concatenate
%     % padval: Value used to pad (e.g., NaN)
%     % dim: Dimension along which to concatenate (1 = rows, 2 = columns)
%     max_length = max(cellfun(@numel, C));
%     M = cellfun(@(x) padarray(x, [max_length - numel(x), 0], padval, 'post'), C, 'UniformOutput', false);
%     M = cell2mat(M);
% end

%% SAVE FIGURES
dirFigout = './TUFFC25/simuPowLaw/v1';
if (~exist(dirFigout)); mkdir (dirFigout); end
titleFigout = 'Fig';
save_all_figures_to_directory(dirFigout, titleFigout)