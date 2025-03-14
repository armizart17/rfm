% Test of calc_powerSpectra function (simplified Chahuara code)
% USING DATA SIMU KWAVE GAUSS ATTEMPT TO CHANGE BSC CHECK KERNEL VARIABLES
% IN SIMULATION
% ATTENUATION alpha_power 1.001, 1.005, 1.05, 1.1
% OPTION AVAILABLE 2DoF with priors and 3-DoF by changing estim_method
% Dec 2025 EMZ

%% FULL
% P20 ACS sample = 0.7036 dB/cm/MHz
% P20 ACS sample = 0.3121 dB/cm/MHz
% -----P20 PowerLaw (b.f^n)-----
% n      = 3.864606
% b      = 5.705978e-06
% b [dB] = -52.436699
% ---------
% -----P22 PowerLaw (b.f^n)-----
% n      = 3.670544
% b      = 5.125957e-07
% b [dB] = -62.902251
% ---------
%%
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
% methods = {'3-DoF', '2-DoF-n'};
% methods = {'2-DoF-s'};

% First row for headers, second for data
bsc_results = cell(2, length(methods)); 
maps_results = cell(2, length(methods));

% Store headers
% bsc_results(1, :) = {sprintf('3-DoF'), sprintf('2-DoF "n"')};
% maps_results(1, :) = {sprintf('3-DoF'), sprintf('2-DoF "n"')};

bsc_results(1, :) = {sprintf('3-DoF'), sprintf('2-DoF "a"'), sprintf('2-DoF "b"'),  sprintf('2-DoF "n"')};
maps_results(1, :) = {sprintf('3-DoF'), sprintf('2-DoF "a"'), sprintf('2-DoF "b"'),  sprintf('2-DoF "n"')};

pathDataUIUC = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\Data UIUC\';
pathEMZdata = [pathDataUIUC,'SiemensProcData_EMZ\'];

%% LOAD DATA
operator = 'ZT'; % CHANGEABLE: AH or ZT
freqValue = 4; % 2.5, 4, 5.5 [MHZ]

numPhantomSam = 20;

numPhantomRef = 22;

if old 
j_sam = 1.0;
j_ref = 1.0;

% SAMPLE SPECS
if numPhantomSam == 20 % P20 ACS = 0.6851 dB/cm/MHz
    b_sam       = 5.705978E-06;
    n_sam       = 3.864606;       
    alpha_sam   = 0.6851; % [dB/cm/MHz]

elseif numPhantomSam == 22 % P22 ACS = 0.3007 dB/cm/MHz
    b_sam       = 5.125957E-07;
    n_sam       = 3.670544;
    alpha_sam   = 0.3007; % [dB/cm/MHz]
end

% REFERENCE SPECS
if numPhantomRef == 20 % P20 ACS = 0.6851 dB/cm/MHz    
    b_ref       = 5.705978E-06;
    n_ref       = 3.864606;    
    alpha_ref   = 0.6851; % [dB/cm/MHz]

elseif numPhantomRef == 22 % P22 ACS = 0.3007 dB/cm/MHz
    b_ref       = 5.125957E-07;
    n_ref       = 3.670544;
    alpha_ref   = 0.3007; % [dB/cm/MHz]
end

else
% SAMPLE SPECS
if numPhantomSam == 20 % P20 ACS = 0.6851 dB/cm/MHz
    b_sam       = 5.705978E-06;
    n_sam       = 3.864606;       
    alpha_sam   = 0.6851; % [dB/cm/MHz]

elseif numPhantomSam == 22 % P22 ACS = 0.3007 dB/cm/MHz
    b_sam       = 5.125957E-07;
    n_sam       = 3.670544;
    alpha_sam   = 0.3007; % [dB/cm/MHz]
end

% REFERENCE SPECS
if numPhantomRef == 20 % P20 ACS = 0.6851 dB/cm/MHz    
    b_ref       = 5.705978E-06;
    n_ref       = 3.864606;    
    alpha_ref   = 0.6851; % [dB/cm/MHz]

elseif numPhantomRef == 22 % P22 ACS = 0.3007 dB/cm/MHz
    b_ref       = 5.125957E-07;
    n_ref       = 3.670544;
    alpha_ref   = 0.3007; % [dB/cm/MHz]
end

end

phantomSamInfo = selectPhantom(operator, numPhantomSam, freqValue);
phantomRefInfo =  selectPhantom(operator, numPhantomRef, freqValue);

%%%%%%%%%%%%%%%%%%%%% RESUMEN %%%%%%%%%%%%%%%%%%%%%
% fprintf('Data resumen\n');

sam = {phantomSamInfo.operator, ['P', num2str(phantomSamInfo.numPhantom)], ...
    phantomSamInfo.freqChar, phantomSamInfo.numCase};

ref = {phantomRefInfo.operator, ['P', num2str(phantomRefInfo.numPhantom)], ...
    phantomRefInfo.freqChar, phantomRefInfo.numCase};

T = table(sam',ref', ...
    'VariableNames', {'SAMPLE','REFERENCE'}, ...
    'RowName',{'Operator','N° Phantom', 'Freq [MHz]', 'N° Case'}); 

disp(T)
clear sam ref

%%%%%%%%%%%%%%%%%%%%% FOLDER NAME %%%%%%%%%%%%%%%%%%%%%

folderNameSam =  ['P', num2str(phantomSamInfo.numPhantom), '_', ...
    'F', phantomSamInfo.freqChar, 'MHz_', ...
     num2str(phantomSamInfo.numCase)];

folderNameRef =  ['P', num2str(phantomRefInfo.numPhantom), '_', ...
    'F', phantomRefInfo.freqChar, 'MHz_', ...
     num2str(phantomRefInfo.numCase)];

%%%%%%%%%%%%%%%%%%%%% READ FOLDER %%%%%%%%%%%%%%%%%%%%%

pathData = [pathEMZdata, 'Operator_', phantomSamInfo.operator, ...
    '/', folderNameSam];
pathRef = [pathEMZdata, 'Operator_', phantomRefInfo.operator, ...
    '/', folderNameRef];

SAM = load([pathData, '/dataSiemens.mat']);
REF = load([pathRef, '/dataSiemens.mat']);

%% SET DATA
SAM.rf = SAM.RF_data(:,:,1);      
SAM.x  = SAM.lateral*1E-3; % was in [mm] set to [m]
SAM.z  = SAM.axial*1E-3;   % was in [mm] set to [m]
SAM.fs = 40E6; % [Hz]

REF.rf = REF.RF_data(:,:,1);      
REF.x  = REF.lateral*1E-3;  % was in [mm] set to [m]
REF.z  = SAM.axial*1E-3;  % was in [mm] set to [m]
REF.fs = 40E6; % [Hz]

% B-MODE CHECK
bmode_sam = db(hilbert(SAM.rf));
bmode_sam = bmode_sam - max(bmode_sam(:));

bmode_ref = db(hilbert(REF.rf));
bmode_ref = bmode_ref - max(bmode_ref(:));


%% SPECTRAL METHOD PARAMETERS
switch freqValue % 2.5, 4, 5.5 [MHZ]
    case 4 
        pars.bw         = [1.65 4.75]; % [MHz]
    case 5.5 
        pars.bw         = [1.6 5.75]; % [MHz]
end

pars.overlap     = 0.8;
pars.blocksize   = 10; % wavelengths
% pars.z_roi       = [3 10]*1E-2; % [m] 
pars.x_roi       = [-5 5]*1E-2; % [m] 

pars.z_roi       = [1 10]*1E-2; % [m] 

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
SR_emz = S_sam ./ S_ref;

SR = permute(SR_emz,[3,1,2]); clear SR_emz

%% GENERAL REGULARIZTATION SETTINGS
% Implementation parameters
par_rpl.tol        = 1e-16;
par_rpl.kmax       = 100;
par_rpl.eps_f      = 1e-16;
par_rpl.m_est      = 0; %Robust
par_rpl.ini_tol    = 1e-16;
par_rpl.df_op      = 0;
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
Z = kron( speye(p*q), log(f) ); % EMZ PowLaw  Size: [p*q*r, p*q] 

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

%% COMPENSATE GAUSS ATTEMPT 2-DoF-n
elseif strcmp( estim_method, '2-DoF-n')

    if (methodsRegu); mu_rpl_tv    = [1E3; 1E3; 1E5]; % [mu_b, mu_n, mu_a]
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
s = +delta_n_prior*ones(p*q, 1);

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

    if (methodsRegu); mu_rpl_tv    = [1E3; 1E3; 1E3]; % [mu_b, mu_n, mu_a]
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

%% SAVE ALL BSC EST GAUSS
% Delta_b*b_ref*f^(delta_n_prior +n_ref)

freq = spectralData_sam.band; % Given choosen BW

% OPTION A
% b_sam     = b_ratio * b_ref;
% n_sam     = n_ratio + n_ref;
% bsc_sam   = b_sam *(freq_bsc.^n_sam);
% med_bsc   = median(bsc_sam(:));

% OPTION B
b_sam_est     = median(g_ratio * b_ref, 'all');
n_sam_est     = median(s_ratio + n_ref, 'all');
bsc_est_powlaw   = b_sam_est *(freq.^n_sam_est);


bsc_results{2, iMet} = bsc_est_powlaw;

%% SAVE ALL MAPS

maps_results{2, iMet} = acs_sam; 
maps_results{3, iMet} = g_ratio_dB; 
maps_results{4, iMet} = s_ratio; 


end

%% BSC  GT

%
fileNameBSC_P20 = 'BSC_vs_freq_P20.txt';
fileNameBSC_P22 = 'BSC_vs_freq_P22.txt';

pathBSC = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\Data UIUC\Siemens_S3000_6C1HD_THI_OFF';

% freqs | bsc

file_BSC_P20 = readmatrix(fullfile(pathBSC,fileNameBSC_P20)); % Freq, BSC
file_BSC_P22 = readmatrix(fullfile(pathBSC,fileNameBSC_P22)); % Freq, BSC

freqs_P20 = file_BSC_P20(:,1);
bsc_P20   = file_BSC_P20(:,2);

freqs_P22 = file_BSC_P22(:,1);
bsc_P22   = file_BSC_P22(:,2);


if numPhantomSam == 20
    % interp measurements for visualization (only see in choosen BW)
    BSC_gt = interp1(freqs_P20, bsc_P20, freq); 

elseif numPhantomSam == 22
    % interp measurements for visualization (only see in choosen BW)
    BSC_gt = interp1(freqs_P22, bsc_P22, freq); 
end

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

if (methodsRegu); ylim([0 1.2])
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
if (methodsRegu)  ylim([9 11])
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
if (methodsRegu)  ylim([0 0.6]); %yticks([0:0.01:0.05])
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

%% PLOTS BSC TOGETHER
xlim_range = [1.5 5.72];
ylim_range = [10^-5 10^-2]; % Y-axis limits

% Define properties for customization
line_width = 2.85; % Set line width
font_size = 16; % Adjust font size

BSC_gt_dB  = 10*log10(BSC_gt);
% diff_fit_dB = @(bsc_pred, bsc_gt) mean ( abs ( 10*log10(bsc_pred) - 10*log10(bsc_gt) ) );
% BSC THEORETICAL

m_3dof.mpe   = mean( (bsc_results{2, 1} - BSC_gt)./BSC_gt );
m_3dof.mae   = mean( abs(bsc_results{2, 1} - BSC_gt)./BSC_gt );
% m_3dof.rmse_dB  = rmse(pow2db(bsc_results{2, 1}), BSC_gt_dB);
m_3dof.rmse  = rmse(bsc_results{2, 1}, BSC_gt);
m_3dof.nrmse = m_3dof.rmse / mean(BSC_gt);
m_3dof.diff_dB = mean ( abs(10*log10(bsc_results{2, 1}) - BSC_gt_dB) ); 

m_2dofa.mpe   = mean( (bsc_results{2, 2} - BSC_gt)./BSC_gt );
m_2dofa.mae   = mean( abs(bsc_results{2, 2} - BSC_gt)./BSC_gt );
% m_2dofa.rmse_dB  = rmse(pow2db(bsc_results{2, 2}), BSC_gt_dB);
m_2dofa.rmse  = rmse(bsc_results{2, 2}, BSC_gt);
m_2dofa.nrmse = m_2dofa.rmse / mean(BSC_gt);
m_2dofa.diff_dB = mean ( abs(10*log10(bsc_results{2, 2}) - BSC_gt_dB) ); 

m_2dofb.mpe   = mean( (bsc_results{2, 3} - BSC_gt)./BSC_gt );
m_2dofb.mae   = mean( abs(bsc_results{2, 3} - BSC_gt)./BSC_gt );
% m_2dofb.rmse_dB  = rmse(pow2db(bsc_results{2, 3}), BSC_gt_dB);
m_2dofb.rmse  = rmse(bsc_results{2, 3}, BSC_gt);
m_2dofb.nrmse = m_2dofb.rmse / mean(BSC_gt);
m_2dofb.diff_dB = mean ( abs(10*log10(bsc_results{2, 3}) - BSC_gt_dB) ); 

m_2dofn.mpe   = mean( (bsc_results{2, 4} - BSC_gt)./BSC_gt );
m_2dofn.mae   = mean( abs(bsc_results{2, 4} - BSC_gt)./BSC_gt );
% m_2dofn.rmse_dB  = rmse(pow2db(bsc_results{2, 4}), BSC_gt_dB);
m_2dofn.rmse  = rmse(bsc_results{2, 4}, BSC_gt);
m_2dofn.nrmse = m_2dofn.rmse / mean(BSC_gt);
m_2dofn.diff_dB = mean ( abs(10*log10(bsc_results{2, 4}) - BSC_gt_dB) ); 

% Extract field names
fields = fieldnames(m_3dof);

% Create a table
Tbsc = table(struct2cell(m_3dof), struct2cell(m_2dofa), ...
    struct2cell(m_2dofb), struct2cell(m_2dofn), 'RowNames', fields, 'VariableNames', methods);


bsc_results{1, 1} = sprintf('3-DoF      (RMSE = %.2f%%) \n', 100*m_3dof.rmse);
bsc_results{1, 2} = sprintf('2-DoF "a" (RMSE = %.2f%%) \n', 100*m_2dofa.rmse);
bsc_results{1, 3} = sprintf('2-DoF "b" (RMSE = %.2f%%) \n', 100*m_2dofb.rmse);
bsc_results{1, 4} = sprintf('2-DoF "n" (RMSE = %.2f%%) \n', 100*m_2dofn.rmse);

% bsc_results{1, 1} = sprintf('3-DoF      (dB = %.2f) \n', m_3dof.diff_dB);
% bsc_results{1, 2} = sprintf('2-DoF "a" (dB = %.2f) \n', m_2dofa.diff_dB);
% bsc_results{1, 3} = sprintf('2-DoF "b" (dB = %.2f) \n', m_2dofb.diff_dB);
% bsc_results{1, 4} = sprintf('2-DoF "n" (dB = %.2f) \n', m_2dofn.diff_dB);
% 


% Convert hexadecimal colors to RGB (MATLAB requires values between 0 and 1)
colot_gt = '#000000'; % Black
color_1  = '#FF0000';  % 3dof
color_2  = '#D95319';  % 2dof a
color_3  = '#0072BD';  % 2dof b
color_4  = '#77AC30';  % 2dof n

% Create figure and plot data
figure, 
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 700, 700]); % [x, y, width, height] in pixels

semilogy(freq, BSC_gt, '-', 'Color', hex2rgb(colot_gt), 'LineWidth', line_width+0.5, 'DisplayName', 'GT');
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
legend('Location', 'southeast', 'FontSize', font_size + 2);
set(gca, 'FontSize', font_size);
hold off;

% Write to Excel
if methodsRegu;   nameExcel = 'metricsBSCregu_phantom.xlsx'; 
else     nameExcel = 'metricsBSCnoregu_phantom.xlsx'; 
end

excelFile = fullfile(dirFigout, nameExcel);

writetable(Tbsc, excelFile, 'Sheet', 'Metrics', 'WriteRowNames', true);

fprintf('Table BSC saved to %s\n', excelFile);

clear m_3dof m_2dofa m_2dofb m_2dofn

%%
% keyboard
%% PLOTS ACS MAPS TOGETHER

% maps_results{1, 1} = '3-DoF';
% acs_results{1, 2} = '2-DoF "a"';
% acs_results{1, 3} = '2-DoF "b"';
% acs_results{1, 4} = '2-DoF "n" ';

fontSize = 16;

figure,
% set(gcf,'units','normalized','outerposition',[0 0.1 1 0.8]); box on;

tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
sgtitle('\alpha', 'FontSize', fontSize+2, 'FontWeight', 'bold');

for iMet = 1:length(methods)
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
    xlabel('Lateral [cm]'), ylabel('Depth [cm]');
    hColor.Label.String = 'dB\cdotcm^{-1}\cdotMHz^{-1}';
    title(maps_results{1, iMet})
    set(gca,'fontsize',fontSize)

end
%%
