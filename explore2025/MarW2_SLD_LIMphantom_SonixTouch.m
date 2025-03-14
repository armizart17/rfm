% SLD PHANTOMS LIM 
% Manual or by code roi selection 

% USING EXPERIMENTAL PHANTOM DATA ACQ BY OLD SONIX TOUCH
% Phantom: 261, acs 0.54 bsc 4.81E-4@3MHz
% Phantom: 544, acs 0.53 bsc 6.73E-4@3MHz
%% GENERAL UTILS
clear all, 
% clc, 
warning('off');
% close all;

addpath(genpath(pwd))

Np2dB       = 20*log10(exp(1));
dB2Np       = 1/Np2dB;
range_bmode = [-60 0];

plotroi      = false;
plotBmode    = true;
manualroi    = true;

coloracs     = false;

plotSpectrumSam = false;
plotSpectrumRef = false;

pathData = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\sonixtouch';

%% Spectral parameters
fLow = 3e6; fHigh = 9e6;
c0               = 1540; % sos [m/s]
fc               = 0.5*(fLow + fHigh); 

pars.overlap     = 0.8;
pars.window_type = 3; %  (1) Hanning, (2) Tuckey, (3) Hamming, (4) Tchebychev
pars.saran_layer = false;

pars.lambda     = c0/fc;
% pars.bw         = [2 12]; % [MHz]
pars.bw         = [3 8.4]; % [MHz]
pars.blocksize  = 10; % wavelengths
pars.ratio_zx   = 1.25;

%% DATA FILES

% folderDataSam = 'ID544V2'; numPhantomSam = 544;
folderDataSam = 'ID261V2'; numPhantomSam = 261; 
% folderDataSam = 'ID316V2'; % targets T6

% folderDataRef = 'ID261V2'; numPhantomRef = 261;
folderDataRef = 'ID544V2'; numPhantomRef = 544;


j_sam = 1.0;
j_ref = 1.0;

% SAMPLE SPECS

switch numPhantomRef
    case 261
        alpha_ref   = 0.54; % [dB/cm/MHz] % manufacturer
        % alpha_ref   = 0.48; % [dB/cm/MHz] % real tested
    case 544
        alpha_ref   = 0.53; % [dB/cm/MHz]
end

filesSam = dir(fullfile(pathData, folderDataSam,'*.mat'));
% samList = ["16-19-52","16-20-50","16-21-22", "16-21-22"];
samName = filesSam(1).name;
SAM     = load( fullfile(pathData, folderDataSam, samName) );
if isfield (SAM, 'RF')
    SAM.rf  = SAM.RF;
end


filesRef = dir(fullfile(pathData, folderDataRef,'*.mat'));
% refList = ["16-27-24","16-27-58","16-28-33", "16-29-07"];
% filesRef = filesRef(x,:); % for select specific "x" filesRef
numRefs  = length(filesRef); 
REF     = load( fullfile(pathData, folderDataRef, filesRef(1).name ) );
newrf  = nan([size(REF.rf), numRefs], 'like', REF.rf); % Use 'like' for type consistency
for i = 1:numRefs
    newrf(:,:,i) = load(fullfile(pathData,folderDataRef,filesRef(1).name ), 'rf').rf(:,:,1); % Directly extract rf, avoiding redundant variables
end

REF.rf = newrf;
REF.acs = alpha_ref;
REF.acs = alpha_ref;

refName = filesRef(1).name;
REF     = load( fullfile(pathData, folderDataRef, refName) );

% B-MODE CHECK
bmode_sam = db(hilbert(SAM.rf(:,:,1)));
bmode_sam = bmode_sam - max(bmode_sam(:));

bmode_ref = db(hilbert(REF.rf(:,:,1)));
bmode_ref = bmode_ref - max(bmode_ref(:));

%% ROI SELECTION
if ~manualroi

    pars.x_roi       = [0.5 37]*1E-3; % [m] 
    % pars.x_roi       = [29.5 38]*1E-3; % [m] right
    % pars.x_roi       = [0 8]*1E-3; % [m] left
    pars.z_roi       = [6 33]*1E-3; % [m] 
    
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

%% PLOT BMODE
if (plotroi)
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
% Plot
% [pxx,fpxx] = pwelch(rfdata_sam_roi-mean(rfdata_sam_roi),nz,nz-wz,nz,fs);
% meanSpectrum = mean(pxx,2);
% meanSpectrum(1) = 0;
% 
% figure,
% plot(fpxx/1e6,db(meanSpectrum/max(meanSpectrum))),grid on
% xlim([0, fs/2e6])
% hold on
% xline(freq_L/1e6, 'k--')
% xline(freq_H/1e6, 'k--')
% hold off
% xlabel('Frequency [MHz]')
% ylabel('Magnitude [dB]')

%% Calculate Power Spectra

% SAMPLE
% spectralData_sam = calc_powerSpectra_prox_dis(SAM, pars);
spectralData_sam = calc_powerSpectraFull_prox_dis(SAM, pars);

% REFERENCE
% spectralData_ref = calc_powerSpectra_prox_dis(REF, pars);
spectralData_ref = calc_powerSpectraFull_prox_dis(REF, pars);

%% PLOT FULL Spectra
[m, n, ~] = size(spectralData_sam.Sfull);

if (plotSpectrumSam)

nLines = 5;
lin_cen = round(n / 2); 
lat_range = max(1, lin_cen-fix(nLines/2)):min(n, lin_cen+fix(nLines/2)); 

S_2d = squeeze(mean(spectralData_sam.Sfull(:, lat_range, :), 2)); % Mean over 2nd dim (lateral)
S_2d_dB  = pow2db(S_2d ./ max(S_2d, [], 2));

figure;
set(gcf,'units','normalized','outerposition',[0 0.1 0.5 0.5]); box on;

subplot(1,2,1)
imagesc(spectralData_sam.bandFull, spectralData_sam.depth*1e3, S_2d_dB),
xline(pars.bw(1), 'w--', 'LineWidth', 2)
xline(pars.bw(2), 'w--', 'LineWidth', 2)
xlim([0 SAM.fs/2]*1e-6); % visualize only positive size
xlabel('Frequency [MHz]');
ylabel('Depth [mm]');
h2 = colorbar; 
ylabel(h2,'dB');
title('SAM Norm Power Spectrum by depth');

subplot(1,2,2)
plot(spectralData_sam.bandFull, S_2d_dB(1, :), 'DisplayName', 'Top')
hold on, grid on;
plot(spectralData_sam.bandFull, S_2d_dB(round(m/2), :), 'DisplayName', 'Half')
plot(spectralData_sam.bandFull, S_2d_dB(end, :), 'DisplayName', 'Bottom')
yline(-20, 'k--', 'DisplayName', '')
xlim([0 SAM.fs/2]*1e-6); % visualize only positive size
hold off;
xlabel('Frequency [MHz]');
ylabel('Norm Power Spectrum [dB]');
title('SAM Norm Power Spectrum');
legend('Location', 'Best');
end

if (plotSpectrumRef)

nLines = 5;
lin_cen = round(n / 2); 
lat_range = max(1, lin_cen-fix(nLines/2)):min(n, lin_cen+fix(nLines/2)); 

S_2d = squeeze(mean(spectralData_ref.Sfull(:, lat_range, :), 2)); % Mean over 2nd dim (lateral)
S_2d_dB  = pow2db(S_2d ./ max(S_2d, [], 2));

figure;
set(gcf,'units','normalized','outerposition',[0 0.1 0.5 0.5]); box on;

subplot(1,2,1)
imagesc(spectralData_ref.bandFull, spectralData_ref.depth*1e3, S_2d_dB),
xline(pars.bw(1), 'w--', 'LineWidth', 2)
xline(pars.bw(2), 'w--', 'LineWidth', 2)
xlim([0 SAM.fs/2]*1e-6); % visualize only positive size
xlabel('Frequency [MHz]');
ylabel('Depth [mm]');
h2 = colorbar; 
ylabel(h2,'dB');
title('REF Norm Power Spectrum by depth');

subplot(1,2,2)
plot(spectralData_ref.bandFull, S_2d_dB(1, :), 'DisplayName', 'Top')
hold on, grid on;
plot(spectralData_ref.bandFull, S_2d_dB(round(m/2), :), 'DisplayName', 'Half')
plot(spectralData_ref.bandFull, S_2d_dB(end, :), 'DisplayName', 'Bottom')
yline(-20, 'k--', 'DisplayName', '')
xlim([0 SAM.fs/2]*1e-6); % visualize only positive size
hold off;
xlabel('Frequency [MHz]');
ylabel('Norm Power Spectrum [dB]');
title('REF Norm Power Spectrum');
legend('Location', 'Best');
end

%% SLD STEP BY STEP

% SAMPLE Sp an Sd
Sp_sam     = spectralData_sam.Sp;
Sd_sam     = spectralData_sam.Sd;
bandFull   = spectralData_sam.bandFull;  % [MHz]
zd_zp      = spectralData_sam.zd_zp;     % [m]

% REFERENCE Sp an Sd
Sp_ref = spectralData_ref.Sp;
Sd_ref = spectralData_ref.Sd;

% COMPENSATION ATTENUATION
acs_ref     = alpha_ref; % [dB/cm/MHz]
att_ref     = acs_ref*bandFull /Np2dB; % vector [Np/cm]
att_ref_map = ones(size(Sp_sam)) .* reshape(att_ref, 1, 1, []); % 3D array as SLogRatio

SLogRatioFull = log(Sp_sam./Sd_sam) - ( log(Sp_ref./Sd_ref) - 4*zd_zp*1E2*att_ref_map ); % (rows, cols, freqChannels)

% USABLE FREQUENCY RANGE
range = bandFull > pars.bw(1) & bandFull < pars.bw(2);
SLogRatio   = SLogRatioFull(:,:,range);
band        = bandFull(range);

SLogRatio_vec = reshape(SLogRatio, [], size(SLogRatio, 3));   
mux_SLogRatio = 1./(abs(mean(SLogRatio_vec, 1, 'omitnan')) ./ std(SLogRatio_vec, 0, 1, 'omitnan') + 1E-5);
weightEstimators = rescale(mux_SLogRatio, 1, max(mux_SLogRatio));
%% REGULARIZED METHODS

%%%%%%%%%%%%%%%%%% RSLD %%%%%%%%%%%%%%%%%%
mu1 = 10^3.5;
mu2 = 10^3.5;
[m, n, p] = size(SLogRatio);
mask = ones(m,n,p);
A1 = kron( 4*zd_zp*1E2*band , speye(m*n) );
A2 = kron( ones(size(band)) , speye(m*n) );
tol = 1e-3;
[Bn,~] = AlterOpti_ADMM(A1, A2, SLogRatio(:),mu1,mu2,m,n,tol,mask(:));
acs_rsld = reshape(Bn*Np2dB,m,n);
%%%%%%%%%%%%%%%%%% RSLD %%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%% TNV Den %%%%%%%%%%%%%%%%%% 
bestTau = 0.001;
mu_wtnv = 2.7696; % best JASA
tol = 0.5e-4; % tolerance error
stableIter = 200;
maxIter = 1000;
SLD_tnv = pdo_den_wtnv(SLogRatio, mu_wtnv, bestTau, maxIter, tol, stableIter, weightEstimators);       
acs_tnv = cgs_ACS([A1 A2], SLD_tnv);
%%%%%%%%%%%%%%%%%% TNV Den %%%%%%%%%%%%%%%%%% 

%% METRICS PLOT

calc2dStats = {@(x) mean(x(:)), @(x) std(x(:)), @(x) 100 * std(x(:)) / mean(x(:))};

met_rsld = cellfun(@(f) f(acs_rsld), calc2dStats);
met_tnv  = cellfun(@(f) f(acs_tnv), calc2dStats);

mm = 1e3;
acs_range = [0.2 1.2];

if coloracs
figure,
subplot (1,2,1)
imagesc(spectralData_sam.lateral*mm, spectralData_sam.depth*mm, acs_rsld, acs_range), colorbar
axis("image");
xlabel('Lateral [mm]'), ylabel('Depth [mm]'), colormap turbo;
title({'RSLD $\alpha_s$:', ...
       [num2str(round(met_rsld(1), 3)), ' $\pm$ ', num2str(round(met_rsld(2), 3)), ', CV = ', num2str(round(met_rsld(3), 3))]}, ...
      'Interpreter', 'latex');
h2 = colorbar; 
ylabel(h2,'dB\cdotcm^{-1}\cdotMHz^{-1}');


subplot (1,2,2)
imagesc(spectralData_sam.lateral*mm, spectralData_sam.depth*mm, acs_tnv, acs_range), colorbar
axis("image");
xlabel('Lateral [mm]'), ylabel('Depth [mm]'), colormap turbo;
title({'TNV $\alpha_s$:', ...
       [num2str(round(met_tnv(1), 3)), ' $\pm$ ', num2str(round(met_tnv(2), 3)), ', CV = ', num2str(round(met_tnv(3), 3))]}, ...
      'Interpreter', 'latex');
h2 = colorbar; 
ylabel(h2,'dB\cdotcm^{-1}\cdotMHz^{-1}');

end
%% SIMPLE OVERLAY ACS MAPS + BMODE 

estim_method = 'RSLD';
fontSize     = 14;
figure,
set(gcf,'units','normalized','outerposition',[0 0.1 1 0.8]); box on;

tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
sgtitle(estim_method, 'FontSize', fontSize+2, 'FontWeight', 'bold');

%%%%%%%%%%%%%%%%%%%%%%%%%% RSLD (ACS) %%%%%%%%%%%%%%%%%%%%%%%%%%

units           = 1E3;
bmodeFull       = bmode_sam;
colorImg        = bigImg(acs_rsld, spectralData_sam.rf_roi);
range_bmode     = [-60 0];
range_img       = [0 1.1];
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

    title({'RSLD $\alpha_s$:', ...
       [num2str(round(met_rsld(1), 3)), ' $\pm$ ', num2str(round(met_rsld(2), 3))]}, ...
      'Interpreter', 'latex');
    set(gca,'fontsize',fontSize)

%%%%%%%%%%%%%%%%%%%%%%%%%% RSLD (ACS) %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% TNV (ACS) %%%%%%%%%%%%%%%%%%%%%%%%%%

units           = 1E3;
bmodeFull       = bmode_sam;
colorImg        = bigImg(acs_tnv, spectralData_sam.rf_roi);
range_bmode     = [-60 0];
range_img       = [0 1.1];
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

    title({'TNV $\alpha_s$:', ...
       [num2str(round(met_tnv(1), 3)), ' $\pm$ ', num2str(round(met_tnv(2), 3))]}, ...
      'Interpreter', 'latex');
    set(gca,'fontsize',fontSize)

%%%%%%%%%%%%%%%%%%%%%%%%%% TNV (ACS) %%%%%%%%%%%%%%%%%%%%%%%%%%


%% INTERP OVERLAY ACS MAPS + BMODE 

% estim_method = 'RSLD';
% fontSize     = 14;
% figure,
% set(gcf,'units','normalized','outerposition',[0 0.1 1 0.8]); box on;
% 
% tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
% sgtitle(estim_method, 'FontSize', fontSize+2, 'FontWeight', 'bold');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%% RSLD (ACS) %%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% units           = 1E3;
% bmodeFull       = bmode_sam;
% colorImg        = acs_rsld;
% range_bmode     = [-60 0];
% range_img       = [0 1.1];
% transparency    = 0.65;
% x_img           = spectralData_sam.lateral*units;
% z_img           = spectralData_sam.depth*units;
% xFull           = SAM.x*units;
% zFull           = SAM.z*units;
% [X, Z] = meshgrid(xFull, zFull);
% roi = and(X >= x_img(1), X <= x_img(end)) & ...
%       and(Z >= z_img(1), Z <= z_img(end));
% t = nexttile;
%     [~,hB,hColor] = imOverlayInterp(bmodeFull, colorImg, range_bmode, range_img, ...
%                         transparency, x_img, z_img, roi, xFull, zFull);   
%     hold on;
%     contour(xFull, zFull, roi, 1,'w--', 'LineWidth', 2)
%     hold off;
%     xlabel('Lateral'), ylabel('Depth');
%     hColor.Label.String = 'dB\cdotcm^{-1}\cdotMHz^{-1}';
% 
%     title({'RSLD $\alpha_s$:', ...
%        [num2str(round(met_rsld(1), 3)), ' $\pm$ ', num2str(round(met_rsld(2), 3))]}, ...
%       'Interpreter', 'latex');
%     set(gca,'fontsize',fontSize)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%% RSLD (ACS) %%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%% TNV (ACS) %%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% units           = 1E3;
% bmodeFull       = bmode_sam;
% colorImg        = acs_tnv;
% range_bmode     = [-60 0];
% range_img       = [0 1.1];
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
% 
%     title({'TNV $\alpha_s$:', ...
%        [num2str(round(met_tnv(1), 3)), ' $\pm$ ', num2str(round(met_tnv(2), 3))]}, ...
%       'Interpreter', 'latex');
%     set(gca,'fontsize',fontSize)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%% TNV (ACS) %%%%%%%%%%%%%%%%%%%%%%%%%%