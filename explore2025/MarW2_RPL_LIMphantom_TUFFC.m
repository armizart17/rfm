% Test of calc_powerSpectra function (simplified Chahuara code)
% USING EXPERIMENTAL PHANTOM DATA ACQ BY EMZ AND C. Soto
% Phantom: 261, acs 0.54 bsc 4.81E-4@3MHz
% Phantom: 544, acs 0.53 bsc 6.73E-4@3MHz

%%
% clear all, 
% clc, 
warning('off');
% close all;

methodsRegu = true;
% methodsRegu = false;
manualroi   = false;

addpath(genpath(pwd))

Np2dB       = 20*log10(exp(1));
dB2Np       = 1/Np2dB;
range_bmode = [-60 0];

plotBmode   = false;
plotBSCdB   = true;
plotMaps    = true;

methods = {'3-DoF', '2-DoF-a', '2-DoF-b', '2-DoF-n'};
% methods = {'3-DoF', '2-DoF-n'};

% First row for headers, second for data
bsc_results  = cell(2, length(methods)); 
maps_results = cell(2, length(methods));

bsc_results(1, :)  = {sprintf('3-DoF'), sprintf('2-DoF "a"'), sprintf('2-DoF "b"'),  sprintf('2-DoF "n"')};
maps_results(1, :) = {sprintf('3-DoF'), sprintf('2-DoF "a"'), sprintf('2-DoF "b"'),  sprintf('2-DoF "n"')};

% bsc_results(1, :)  = {sprintf('3-DoF'), sprintf('2-DoF "n"')};
% maps_results(1, :) = {sprintf('3-DoF'), sprintf('2-DoF "n"')};

pathData = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\SavedDataQUSPhantom\bf';

%% LOAD DATA

% folderDataSam = 'targets';numPhantomSam = 'T6';
% folderDataSam = '544'; numPhantomSam = '544';
folderDataSam = '261'; numPhantomSam = '261';

folderDataRef = '544'; numPhantomRef = '544';
% folderDataRef = '261'; numPhantomRef = '261';

j_sam = 1.0;
j_ref = 1.0;

% SAMPLE SPECS
switch str2double(numPhantomSam)
    case 261
        alpha_sam   = 0.47; % [dB/cm/MHz]
    case 544
        alpha_sam   = 0.53; % [dB/cm/MHz]
end

switch str2double(numPhantomRef)
    case 261
        % alpha_ref   = 0.54; % [dB/cm/MHz] % manufacturer
        alpha_ref   = 0.47; % [dB/cm/MHz] % real tested
    case 544
        alpha_ref   = 0.53; % [dB/cm/MHz]
end

% delta_alpha_prior   = alpha_sam - alpha_ref; % [dB/cm/MHz]
delta_alpha_prior   = 0; % [dB/cm/MHz] 
delta_b_prior       = log(0.165347); % from RPM method
delta_n_prior       = 0.684336; % from RPM method

% STABLE
samEnd = ["_F","_F_2","_F_3"];
samName = numPhantomSam + samEnd(1);
SAM = load (fullfile(pathData, folderDataSam, samName));
SAM.acs = alpha_sam;

% TBD
% filesSam = dir(fullfile(pathData, folderDataSam,'*.mat'));
% samName = filesSam(1).name;
% SAM     = load( fullfile(pathData, folderDataSam, samName) );

% SAM.acs = alpha_sam;

filesRef = dir(fullfile(pathData, folderDataRef,'*.mat'));
% filesRef = filesRef(x,:); % for select specific "x" filesRef
numRefs  = length(filesRef); 
REF     = load( fullfile(pathData, folderDataRef, filesRef(1).name ) );
newrf  = nan([size(REF.rf), numRefs], 'like', REF.rf); % Use 'like' for type consistency
for i = 1:numRefs
    newrf(:,:,i) = load(fullfile(pathData,folderDataRef,filesRef(i).name ), 'rf').rf(:,:,1); % Directly extract rf, avoiding redundant variables
end

REF.rf  = newrf;
REF.acs = alpha_ref;

% B-MODE CHECK
bmode_sam = db(hilbert(SAM.rf(:,:,1)));
bmode_sam = bmode_sam - max(bmode_sam(:));

bmode_ref = db(hilbert(REF.rf(:,:,1)));
bmode_ref = bmode_ref - max(bmode_ref(:));

%% ROI SELECTION

if ~manualroi

    pars.x_roi       = [-15 15]*1E-3; % [m] 
    pars.z_roi       = [6 33]*1E-3; % [m] % Same than sonix
    % pars.z_roi       = [30 50]*1E-3; % [m] 
    
else 

    figure('Units','centimeters', 'Position',[5 5 15 15]),
    imagesc(SAM.x*1E3, SAM.z*1E3,bmode_sam,range_bmode);
    colormap gray; clim(range_bmode);
    hb2=colorbar; ylabel(hb2,'dB')
    xlabel('Lateral [mm]'), ylabel('Depth [mm]'); 
    title('Bmode')
    
    confirmation = '';
    while ~strcmp(confirmation,'Yes')
        rect = getrect;
        confirmation = questdlg('Sure?');
        if strcmp(confirmation,'Cancel')
            disp(rect)
            break
        end
    end
    close,

    pars.x_roi     = [rect(1), rect(1)+rect(3)]*1E-3; % [m]
    pars.z_roi     = [rect(2), rect(2)+rect(4)]*1E-3; % [m]
end

%% SPECTRAL METHOD PARAMETERS

pars.P           = 256; % NFFT for 10wl is 256, 20wl 512
pars.bw          = [3.5 8.25]; % [MHz]
pars.overlap     = 0.8;
pars.blocksize   = 10; % wavelengths
pars.window_type = 3; %  (1) Hanning, (2) Tuckey, (3) Hamming, (4) Tchebychev
pars.saran_layer = false;

if (plotBmode)
deadBand = 0.1e-2;
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
% ylim([deadBand*1000 50])
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
% ylim([deadBand*1000 50])
end

%% POWER SPECTRA ESTIMATION
spectralData_sam = calc_powerSpectra(SAM, pars);
% spectralData_sam = calc_powerSpectra_vSimple(SAM, pars); % @
S_sam = spectralData_sam.powerSpectra;

num_ref = 1;

spectralData_ref = calc_powerSpectra(REF, pars);
% spectralData_ref = calc_powerSpectra_vSimple(REF, pars); % @
S_ref = spectralData_ref.powerSpectra;

%% PLOT ALL BW


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
title({'$\Delta b$:', ...
       [num2str(round(m_g, 3)), ' $\pm$ ', num2str(round(s_g, 3)), ', CV = ', num2str(round(cv_g, 3))]}, ...
      'Interpreter', 'latex');

set(gca,'fontsize',fontSize)

subplot(1,3,3)
imagesc(Xaxis*cm, Zaxis*cm, s_ratio), colorbar
axis("image");
xlabel('Lateral [cm]'), colormap turbo;
% title(['$\Delta s$: ', num2str(round(m_s, 3)), ' $\pm$ ', num2str(round(s_s, 3)), ', CV = ', num2str(round(cv_s, 3))], ...
%       'Interpreter', 'latex');
title({'$\Delta n$:', ...
       [num2str(round(m_s, 3)), ' $\pm$ ', num2str(round(s_s, 3)), ', CV = ', num2str(round(cv_s, 3))]}, ...
      'Interpreter', 'latex');
set(gca,'fontsize',fontSize)

end

%% FIGURE SIMPLE OVERLAY BMODE, ACS, DELTA BSC and DELTA N

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
colorImg        = bigImg(acs_sam, spectralData_sam.rf_roi);
range_bmode     = [-60 0];
range_img       = [0.1 1.2];
transparency    = 0.65;
x_img           = spectralData_sam.x_roi*units;
z_img           = spectralData_sam.z_roi*units;
xFull           = SAM.x*units;
zFull           = SAM.z*units;

t = nexttile;
    [~,hB,hColor] = imOverlaySimple(bmodeFull, colorImg, range_bmode, range_img, ...
                        transparency, x_img, z_img, xFull, zFull);   
    hold on;
    rectangle('Position', units*[pars.x_roi(1) pars.z_roi(1) pars.x_roi(2)-pars.x_roi(1) pars.z_roi(2)-pars.z_roi(1)], ...
        'EdgeColor','w', 'LineWidth', 2, 'LineStyle','--'), 
    hold off;
    xlabel('Lateral'), ylabel('Depth');
    hColor.Label.String = 'dB\cdotcm^{-1}\cdotMHz^{-1}';
    title('$\alpha_s$', 'Interpreter', 'latex')
    set(gca,'fontsize',fontSize)
%%%%%%%%%%%%%%%%%%%%%%%%%% Delta alpha (ACS) %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% Delta b (BSC) %%%%%%%%%%%%%%%%%%%%%%%%%%

units           = 1E3;
bmodeFull       = bmode_sam;
colorImg        = bigImg(g_ratio, spectralData_sam.rf_roi);
range_bmode     = [-60 0];
range_img       = [];
transparency    = 0.65;
x_img           = spectralData_sam.x_roi*units;
z_img           = spectralData_sam.z_roi*units;
xFull           = SAM.x*units;
zFull           = SAM.z*units;

    if plotBSCdB 
       colorImg = bigImg(g_ratio_dB, spectralData_sam.rf_roi);
    end

t = nexttile;
    [~,hB,hColor] = imOverlaySimple(bmodeFull, colorImg, range_bmode, range_img, ...
                        transparency, x_img, z_img, xFull, zFull);
    hold on;
    rectangle('Position', units*[pars.x_roi(1) pars.z_roi(1) pars.x_roi(2)-pars.x_roi(1) pars.z_roi(2)-pars.z_roi(1)], ...
        'EdgeColor','w', 'LineWidth', 2, 'LineStyle','--'), 
    hold off;
    xlabel('Lateral'), ylabel('Depth');
    hColor.Label.String = '';
        if plotBSCdB 
            hColor.Label.String ='dB';
        end
    % title('$\frac{b_s}{b_r}$', 'Interpreter','latex')
    title('$\Delta b$', 'Interpreter','latex')
    set(gca,'fontsize',fontSize)
%%%%%%%%%%%%%%%%%%%%%%%%%% Delta g ratio in dB %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% Delta n %%%%%%%%%%%%%%%%%%%%%%%%%%

units           = 1E3;
bmodeFull       = bmode_sam;
colorImg        = bigImg(s_ratio, spectralData_sam.rf_roi);
range_bmode     = [-60 0];
range_img       = [];
transparency    = 0.65;
x_img           = spectralData_sam.x_roi*units;
z_img           = spectralData_sam.z_roi*units;
xFull           = SAM.x*units;
zFull           = SAM.z*units;

t = nexttile;
    [~,hB,hColor] = imOverlaySimple(bmodeFull, colorImg, range_bmode, range_img, ...
                        transparency, x_img, z_img, xFull, zFull);
    hold on;
    rectangle('Position', units*[pars.x_roi(1) pars.z_roi(1) pars.x_roi(2)-pars.x_roi(1) pars.z_roi(2)-pars.z_roi(1)], ...
        'EdgeColor','w', 'LineWidth', 2, 'LineStyle','--'), 
    hold off;
    xlabel('Lateral'), ylabel('Depth');
    hColor.Label.String = '';
    title('$\Delta n$', 'Interpreter','latex')
    set(gca,'fontsize',fontSize)
%%%%%%%%%%%%%%%%%%%%%%%%%% Delta s %%%%%%%%%%%%%%%%%%%%%%%%%%
end

%% FIGURE INTERP OVERLAY BMODE, ACS, DELTA BSC, DELTA N
% if plotBmode
% fontSize = 16;
% 
% figure,
% set(gcf,'units','normalized','outerposition',[0 0.1 1 0.8]); box on;
% 
% tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
% sgtitle(estim_method, 'FontSize', fontSize+2, 'FontWeight', 'bold');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%% alpha_s (ACS) %%%%%%%%%%%%%%%%%%%%%%%%%%
% acs_sam = alpha_ratio + alpha_ref;
% 
% units           = 1E3;
% bmodeFull       = bmode_sam;
% colorImg        = acs_sam;
% range_bmode     = [-100 0];
% range_img       = [0.1 1.2];
% transparency    = 0.65;
% x_img           = spectralData_sam.lateral*units;
% z_img           = spectralData_sam.depth*units;
% xFull           = SAM.x*units;
% zFull           = SAM.z*units;
% [X, Z] = meshgrid(xFull, zFull);
% roi = and(X >= x_img(1), X <= x_img(end)) & ...
%       and(Z >= z_img(1), Z <= z_img(end));
% 
% t = nexttile;
%     [~,hB,hColor] = imOverlayInterp(bmodeFull, colorImg, range_bmode, range_img, ...
%                         transparency, x_img, z_img, roi, xFull, zFull);   
%     hold on;
%     contour(xFull, zFull, roi, 1,'w--', 'LineWidth', 2)
%     hold off;
%     xlabel('Lateral'), ylabel('Depth');
%     hColor.Label.String = 'dB\cdotcm^{-1}\cdotMHz^{-1}';
%     title('$\alpha_s$', 'Interpreter', 'latex')
%     set(gca,'fontsize',fontSize)
% %%%%%%%%%%%%%%%%%%%%%%%%%% Delta alpha (ACS) %%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%% Delta g (BSC) %%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% units           = 1E3;
% bmodeFull       = bmode_sam;
% colorImg        = g_ratio;
% range_bmode     = [-60 0];
% range_img       = [];
% transparency    = 0.65;
% x_img           = spectralData_sam.lateral*units;
% z_img           = spectralData_sam.depth*units;
% xFull           = SAM.x*units;
% zFull           = SAM.z*units;
% [X, Z] = meshgrid(xFull, zFull);
% roi = and(X >= x_img(1), X <= x_img(end)) & ...
%       and(Z >= z_img(1), Z <= z_img(end));
% 
%     if plotBSCdB 
%        colorImg = g_ratio_dB;
%     end
% 
% t = nexttile;
%     [~,hB,hColor] = imOverlayInterp(bmodeFull, colorImg, range_bmode, range_img, ...
%                         transparency, x_img, z_img, roi, xFull, zFull);
%     hold on;
%     contour(xFull, zFull, roi, 1,'w--', 'LineWidth', 2)
%     hold off;
%     xlabel('Lateral'), ylabel('Depth');
%     hColor.Label.String = '';
%         if plotBSCdB 
%             hColor.Label.String ='dB';
%         end
%     title('$\frac{b_s}{b_r}$', 'Interpreter','latex')
%     set(gca,'fontsize',fontSize)
% %%%%%%%%%%%%%%%%%%%%%%%%%% Delta g ratio in dB %%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%% Delta s %%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% units           = 1E3;
% bmodeFull       = bmode_sam;
% colorImg        = s_ratio;
% range_bmode     = [-60 0];
% range_img       = [];
% transparency    = 0.65;
% x_img           = spectralData_sam.lateral*units;
% z_img           = spectralData_sam.depth*units;
% xFull           = SAM.x*units;
% zFull           = SAM.z*units;
% [X, Z] = meshgrid(xFull, zFull);
% roi = and(X >= x_img(1), X <= x_img(end)) & ...
%       and(Z >= z_img(1), Z <= z_img(end));
% 
% 
% t = nexttile;
%     [~,hB,hColor] = imOverlayInterp(bmodeFull, colorImg, range_bmode, range_img, ...
%                         transparency, x_img, z_img, roi, xFull, zFull);
%     hold on;
%     contour(xFull, zFull, roi, 1,'w--', 'LineWidth', 2)
%     hold off;
%     xlabel('Lateral'), ylabel('Depth');
%     hColor.Label.String = '';
%     title('$\Delta n$', 'Interpreter','latex')
%     set(gca,'fontsize',fontSize)
% %%%%%%%%%%%%%%%%%%%%%%%%%% Delta s %%%%%%%%%%%%%%%%%%%%%%%%%%
% end

%% SAVE ALL BSC EST GAUSS

% Delta_b*b_ref*f^(delta_n_prior +n_ref)

freq = spectralData_sam.band; % Given choosen BW

% OPTION A
% b_sam     = b_ratio * b_ref;
% n_sam     = n_ratio + n_ref;
% bsc_sam   = b_sam *(freq_bsc.^n_sam);
% med_bsc   = median(bsc_sam(:));

% OPTION B
b_sam_est      = median(g_ratio, 'all');
n_sam_est      = median(s_ratio, 'all');
bsc_est_powlaw = b_sam_est *(freq.^n_sam_est);

bsc_results{2, iMet} = bsc_est_powlaw;

%% SAVE ALL MAPS

maps_results{2, iMet} = acs_sam; 
maps_results{3, iMet} = g_ratio_dB; 
maps_results{4, iMet} = s_ratio; 


end

%%
% keyboard

%% NEW DISTRIBUTION DELTAs** BOX PLOT
% Define method labels

delta_b_theo = exp(delta_b_prior);
delta_n_theo = delta_n_prior;
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

% if (methodsRegu); ylim([-0.19 1.81])
% else              ylim([-80 80])
% end
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
% yline(0, 'k--')
% if (methodsRegu); %ylim([-0.5 0.5])
% else              ylim([-60 20])
% end

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
% if (methodsRegu); ylim([-0.25 0.05]); %yticks([0:0.01:0.05])
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

%% METRICS BSC & PLOTS BSC TOGETHER

BSC_gt     = bsc_rpm_powlaw;
xlim_range = pars.bw + [-0.01 0.01];
ylim_range = [10^-1 10^1]; % Y-axis limits

% Define properties for customization
line_width = 2.85; % Set line width
font_size = 16; % Adjust font size

BSC_gt_dB  = 10*log10(BSC_gt);
diff_fit_dB = @(bsc_pred, bsc_gt) mean ( abs ( 10*log10(bsc_pred) - 10*log10(bsc_gt) ) );

clear m_3dof m_2dofa m_2dofb m_2dofn MetricsBSC

m_3dof  = get_metrics_homo_gt(bsc_results{2, 1}, true(size(bsc_results{2, 1})), BSC_gt, '3-DoF');
m_3dof.diff_dB = diff_fit_dB(bsc_results{2, 1}, BSC_gt);
m_3dof.param = 'BSC';

m_2dofa  = get_metrics_homo_gt(bsc_results{2, 2}, true(size(bsc_results{2, 2})), BSC_gt, '2-DoF-a');
m_2dofa.diff_dB = diff_fit_dB(bsc_results{2, 2}, BSC_gt);
m_2dofa.param = 'BSC';

m_2dofb  = get_metrics_homo_gt(bsc_results{2, 3}, true(size(bsc_results{2, 3})), BSC_gt, '2-DoF-b');
m_2dofb.diff_dB = diff_fit_dB(bsc_results{2, 3}, BSC_gt);
m_2dofb.param = 'BSC';

m_2dofn  = get_metrics_homo_gt(bsc_results{2, 4}, true(size(bsc_results{2, 4})), BSC_gt, '2-DoF-n');
m_2dofn.diff_dB = diff_fit_dB(bsc_results{2, 4}, BSC_gt);
m_2dofn.param = 'BSC';

MetricsBSC(1:4) = [m_3dof; m_2dofa; m_2dofb; m_2dofn]; 

Tbsc        = struct2table(MetricsBSC);
Tbsc.method = categorical(Tbsc.method);
Tbsc.param  = categorical(Tbsc.param);

%
bsc_results{1, 1} = sprintf('3-DoF     (NRMSE = %.2f%%) \n', 100*MetricsBSC(1).rmse_homo);
bsc_results{1, 2} = sprintf('2-DoF "a" (NRMSE = %.2f%%) \n', 100*MetricsBSC(2).rmse_homo);
bsc_results{1, 3} = sprintf('2-DoF "b" (NRMSE = %.2f%%) \n', 100*MetricsBSC(3).rmse_homo);
bsc_results{1, 4} = sprintf('2-DoF "n" (NRMSE = %.2f%%) \n', 100*MetricsBSC(4).rmse_homo);

% Convert hexadecimal colors to RGB (MATLAB requires values between 0 and 1)
colot_gt = '#000000';  % Black
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
% ylim(ylim_range);
xlim(xlim_range);
title('Backscatter Coefficient (BSC)', 'FontSize', font_size + 2);
legend('Location', 'best', 'FontSize', font_size + 2);
set(gca, 'FontSize', font_size);
hold off;

%%
% keyboard

%% METRICS MAPS a,b,n 

clear m_3dof m_2dofa m_2dofb m_2dofn MetricsParam

%%%%%%%%%%%%%%%%%%% Metricas a %%%%%%%%%%%%%%%%%%%
m_3dof  = get_metrics_homo_gt(maps_results{2, 1}, true(size(maps_results{2, 1})), alpha_sam, '3-DoF');
m_2dofa = get_metrics_homo_gt(maps_results{2, 2}, true(size(maps_results{2, 1})), NaN, '2-DoF-a');
m_2dofb = get_metrics_homo_gt(maps_results{2, 3}, true(size(maps_results{2, 1})), alpha_sam, '2-DoF-b');
m_2dofn = get_metrics_homo_gt(maps_results{2, 4}, true(size(maps_results{2, 1})), alpha_sam, '2-DoF-n');
m_3dof.param  = 'a';
m_2dofa.param = 'a';
m_2dofb.param = 'a';
m_2dofn.param = 'a';

MetricsParam(1:4) = [m_3dof; m_2dofa; m_2dofb; m_2dofn];
%%%%%%%%%%%%%%%%%%% Metricas a %%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% Metricas b %%%%%%%%%%%%%%%%%%%
m_3dof  = get_metrics_homo_gt(maps_results{3, 1}, true(size(maps_results{3, 1})), pow2db(delta_b_theo), '3-DoF');
m_2dofa = get_metrics_homo_gt(maps_results{3, 2}, true(size(maps_results{3, 2})), pow2db(delta_b_theo), '2-DoF-a');
m_2dofb = get_metrics_homo_gt(maps_results{3, 3}, true(size(maps_results{3, 3})), NaN, '2-DoF-b');
m_2dofn = get_metrics_homo_gt(maps_results{3, 4}, true(size(maps_results{3, 4})), pow2db(delta_b_theo), '2-DoF-n');
m_3dof.param  = 'b';
m_2dofa.param = 'b';
m_2dofb.param = 'b';
m_2dofn.param = 'b';

MetricsParam(5:8) = [m_3dof; m_2dofa; m_2dofb; m_2dofn];
%%%%%%%%%%%%%%%%%%% Metricas b %%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% Metricas n %%%%%%%%%%%%%%%%%%%
m_3dof  = get_metrics_homo_gt(maps_results{4, 1}, true(size(maps_results{4, 1})), delta_n_theo, '3-DoF');
m_2dofa = get_metrics_homo_gt(maps_results{4, 2}, true(size(maps_results{4, 2})), delta_n_theo, '2-DoF-a');
m_2dofb = get_metrics_homo_gt(maps_results{4, 3}, true(size(maps_results{4, 3})), delta_n_theo, '2-DoF-b');
m_2dofn = get_metrics_homo_gt(maps_results{4, 4}, true(size(maps_results{4, 4})), NaN, '2-DoF-n');
m_3dof.param  = 'n';
m_2dofa.param = 'n';
m_2dofb.param = 'n';
m_2dofn.param = 'n';

MetricsParam(9:12) = [m_3dof; m_2dofa; m_2dofb; m_2dofn];
%%%%%%%%%%%%%%%%%%% Metricas n %%%%%%%%%%%%%%%%%%%

Tparam        = struct2table(MetricsParam);
Tparam.method = categorical(Tparam.method);
Tparam.param  = categorical(Tparam.param);

%% PLOTS ACS MAPS TOGETHER

fontSize = 16;

figure,
% set(gcf,'units','normalized','outerposition',[0 0.1 1 0.8]); box on;

tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'tight');
sgtitle('\alpha', 'FontSize', fontSize+2, 'FontWeight', 'bold');

for iMet = 1:length(methods)
%%%%%%%%%%%%%%%%%%%%%%%%%% alpha_s (ACS) %%%%%%%%%%%%%%%%%%%%%%%%%%
acs_sam = alpha_ratio + alpha_ref;

units           = 1E3;
bmodeFull       = bmode_sam;
colorImg        = maps_results{2, iMet};
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
    xlabel('Lateral [mm]'), ylabel('Depth [mm]');
    hColor.Label.String = 'dB\cdotcm^{-1}\cdotMHz^{-1}';

    title([maps_results{1, iMet}, sprintf(': %.2f±%.2f', mean2d(maps_results{2, iMet}), std2d(maps_results{2, iMet}))])

    set(gca,'fontsize',fontSize)

end
%%
%% BSC  GT
% keyboard
nameBSC = strcat('BSC_', sprintf('sam%.2f_ref%.2f', alpha_sam, alpha_ref ));
nameBSC = strrep(nameBSC, '.', 'p');
% BSC = calculateBSC_RPM_ok(SAM, REF, pars); % slow only once**

BSC = calculateBSC_RPM_fast(SAM, REF, pars); % slow only once**

% bsc_rpm = BSC.BSCcurve_Uni(:,1); % mean
bsc_rpm = BSC.BSCcurve_Uni(:,2); % median
freq    = BSC.band;

% MM = [log(freq), ones(size(freq))];
% rr = log(bsc_rpm);
% pp = cgs(MM'*MM,MM'*rr,1e-16,10000); % coeffs
% b = exp(pp(1));
% n = pp(2);


% Perform linear regression  ln(bsc) = d_n . ln(f) + ln(d_b) 
coeffs   = polyfit(log(freq), log(bsc_rpm), 1); % Fit y = mx + c
d_n      = coeffs(1); % Slope = d_n
ln_db    = coeffs(2); % Intercept = ln(d_b) 
d_b      = exp(ln_db); % 

% Display results
fprintf('-----RPM PowLaw (b.(f.^n))-----\n')
fprintf('Δb           = %f\n', d_b);
fprintf('b_s/b_r [dB] = %f\n', 10*log10(d_b));
fprintf('Δn           = %f\n', d_n);
fprintf('---------\n')

bsc_rpm_powlaw = d_b*(freq.^d_n);


%%
vertical = false;

fontSize = 16;

for iMet = 1:length(methods)

estim_method = maps_results{1, iMet};

figure,

if vertical
set(gcf,'units','normalized','outerposition',[0 0.05 0.2 0.9]); box on;
tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'tight');

else
set(gcf,'units','normalized','outerposition',[0 0.05 0.6 0.4]); box on;
tiledlayout(1, 4, 'TileSpacing', 'compact', 'Padding', 'tight');
end

sgtitle(estim_method, 'FontSize', fontSize+2, 'FontWeight', 'bold');

%%%%%%%%%%%%%%%%%%%%%%%%%% Bmode %%%%%%%%%%%%%%%%%%%%%%%%%%
units           = 1E3;
tb = nexttile;
imagesc(SAM.x*units, SAM.z*units, bmode_sam, range_bmode)
axis image
% xlim([spectralData_sam.x_roi(1) spectralData_sam.x_roi(end)]*units),
% ylim([spectralData_sam.z_roi(1) spectralData_sam.z_roi(end)]*units),
xlabel('Lateral'), ylabel('Depth');
c = colorbar; c.Label.String = 'dB';
title('Bmode')
colormap(tb,gray)
set(gca,'fontsize',fontSize)

%%%%%%%%%%%%%%%%%%%%%%%%%% Bmode %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% alpha_s (ACS) %%%%%%%%%%%%%%%%%%%%%%%%%%
acs_sam = maps_results{2, iMet};

units           = 1E3;
bmodeFull       = bmode_sam;
colorImg        = bigImg(acs_sam, spectralData_sam.rf_roi);
range_bmode     = [-60 0];
range_img       = [0.1 1.2];
transparency    = 0.65;
x_img           = spectralData_sam.x_roi*units;
z_img           = spectralData_sam.z_roi*units;
xFull           = SAM.x*units;
zFull           = SAM.z*units;

t = nexttile;
    [~,hB,hColor] = imOverlaySimple(bmodeFull, colorImg, range_bmode, range_img, ...
                        transparency, x_img, z_img, xFull, zFull);   
    hold on;
    rectangle('Position', units*[pars.x_roi(1) pars.z_roi(1) pars.x_roi(2)-pars.x_roi(1) pars.z_roi(2)-pars.z_roi(1)], ...
        'EdgeColor','w', 'LineWidth', 2, 'LineStyle','--'), 
    hold off;
    xlabel('Lateral'), ylabel('Depth');
    hColor.Label.String = 'dB\cdotcm^{-1}\cdotMHz^{-1}';
    title('$\alpha_s$', 'Interpreter', 'latex')
    colormap(t,turbo)
    set(gca,'fontsize',fontSize)
%%%%%%%%%%%%%%%%%%%%%%%%%% Delta alpha (ACS) %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% Delta b (BSC) %%%%%%%%%%%%%%%%%%%%%%%%%%
g_ratio_dB = maps_results{3, iMet};

units           = 1E3;
bmodeFull       = bmode_sam;
colorImg        = bigImg(g_ratio_dB, spectralData_sam.rf_roi);
range_bmode     = [-60 0];
range_img       = [];
transparency    = 0.65;
x_img           = spectralData_sam.x_roi*units;
z_img           = spectralData_sam.z_roi*units;
xFull           = SAM.x*units;
zFull           = SAM.z*units;

    if plotBSCdB 
       colorImg = bigImg(g_ratio_dB, spectralData_sam.rf_roi);
    end

t = nexttile;
    [~,hB,hColor] = imOverlaySimple(bmodeFull, colorImg, range_bmode, range_img, ...
                        transparency, x_img, z_img, xFull, zFull);
    hold on;
    rectangle('Position', units*[pars.x_roi(1) pars.z_roi(1) pars.x_roi(2)-pars.x_roi(1) pars.z_roi(2)-pars.z_roi(1)], ...
        'EdgeColor','w', 'LineWidth', 2, 'LineStyle','--'), 
    hold off;
    xlabel('Lateral'), ylabel('Depth');
    hColor.Label.String = '';
        if plotBSCdB 
            hColor.Label.String ='dB';
        end
    % title('$\frac{b_s}{b_r}$', 'Interpreter','latex')
    title('$\Delta b$', 'Interpreter','latex')
    colormap(t,turbo)
    set(gca,'fontsize',fontSize)
%%%%%%%%%%%%%%%%%%%%%%%%%% Delta g ratio in dB %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% Delta n %%%%%%%%%%%%%%%%%%%%%%%%%%
s_ratio = maps_results{4, iMet};
units           = 1E3;
bmodeFull       = bmode_sam;
colorImg        = bigImg(s_ratio, spectralData_sam.rf_roi);
range_bmode     = [-60 0];
range_img       = [];
transparency    = 0.65;
x_img           = spectralData_sam.x_roi*units;
z_img           = spectralData_sam.z_roi*units;
xFull           = SAM.x*units;
zFull           = SAM.z*units;

t = nexttile;
    [~,hB,hColor] = imOverlaySimple(bmodeFull, colorImg, range_bmode, range_img, ...
                        transparency, x_img, z_img, xFull, zFull);
    hold on;
    rectangle('Position', units*[pars.x_roi(1) pars.z_roi(1) pars.x_roi(2)-pars.x_roi(1) pars.z_roi(2)-pars.z_roi(1)], ...
        'EdgeColor','w', 'LineWidth', 2, 'LineStyle','--'), 
    hold off;
    xlabel('Lateral'), ylabel('Depth');
    hColor.Label.String = '[a.u.]';
    title('$\Delta n$', 'Interpreter','latex')
    colormap(t,turbo)
    set(gca,'fontsize',fontSize)

colormap(tb,gray)
%%%%%%%%%%%%%%%%%%%%%%%%%% Delta s %%%%%%%%%%%%%%%%%%%%%%%%%%

end