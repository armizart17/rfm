% SLD PHANTOMS LIM 
% Manual or by code roi selection 

% USING EXPERIMENTAL PHANTOM DATA ACQ BY EMZ AND C. Soto
% Phantom: 261, acs 0.54 bsc 4.81E-4@3MHz
% Phantom: 544, acs 0.53 bsc 6.73E-4@3MHz
% TARGETS options: T1_F T2_F T3_F T4_F T5_F T6_F T7_F T8F
% PHANTOMS options: 544 

%% GENERAL UTILS
% clear all, 
% clc, 
warning('off');
% close all;

addpath(genpath(pwd))

Np2dB       = 20*log10(exp(1));
dB2Np       = 1/Np2dB;
range_bmode = [-60 0];

plotBmode   = false;
manualroi   = true;

pathData = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\SavedDataQUSPhantom';

%% DATA FILES

samName = 'T6_F';

folderDataSam = 'bf\targets';
folderDataRef = 'bf';

% numPhantomSam = 544;
% numPhantomRef = 261;

% numPhantomSam = 261;
numPhantomRef = 544;

j_sam = 1.0;
j_ref = 1.0;

% SAMPLE SPECS
if numPhantomRef == 544 
    % alpha_sam           = 0.54; % [dB/cm/MHz]
    alpha_ref           = 0.53; % [dB/cm/MHz]

elseif numPhantomRef == 261
    % alpha_sam           = 0.53; % [dB/cm/MHz]
    alpha_ref           = 0.54; % [dB/cm/MHz]

end

% samName = sprintf('%d_F_3', numPhantomSam);

SAM     = load( fullfile(pathData, folderDataSam, samName) );

% refName = sprintf('%d_F_2', numPhantomRef);

refEnd = ["","_2","_3"];
numRefs  = length(refEnd);
REF     = load( fullfile(pathData, folderDataRef, numPhantomRef+"_F"+refEnd(1)) );
newrf  = nan([size(REF.rf), numRefs], 'like', REF.rf); % Use 'like' for type consistency
for i = 1:numRefs
    newrf(:,:,i) = load(fullfile(pathData,folderDataRef,numPhantomRef+"_F"+refEnd(i)), 'rf').rf; % Directly extract rf, avoiding redundant variables
end
REF.rf = newrf;
REF.acs = alpha_ref;

%%

% B-MODE CHECK
bmode_sam = db(hilbert(SAM.rf));
bmode_sam = bmode_sam - max(bmode_sam(:));

bmode_ref = db(hilbert(REF.rf));
bmode_ref = bmode_ref - max(bmode_ref(:));

%% ROI SELECTION
if ~manualroi

    pars.x_roi       = [-15 15]*1E-3; % [m] 
    pars.z_roi       = [10 50]*1E-3; % [m] 
    
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

%% SLD STEP by STEP

% JT
pars.bw         = [2 12]; % [MHz]
% pars.bw         = [3 8.4]; % [MHz]
pars.blocksize  = 12; % wavelengths
pars.ratio_zx   = 1.5;
pars.window_type = 3; %  (1) Hanning, (2) Tuckey, (3) Hamming, (4) Tchebychev
pars.saran_layer = false;

% SAMPLE
spectralData_sam = calc_powerSpectra_prox_dis(SAM, pars);

Sp_sam = spectralData_sam.Sp;
Sd_sam = spectralData_sam.Sd;

% REFERENCE
num_ref = 1;
spectralData_ref = calc_powerSpectra_prox_dis(REF, pars);

band   = spectralData_sam.band; % [MHz]
zd_zp  = spectralData_sam.zd_zp; % [m]
Sp_ref = spectralData_ref.Sp;
Sd_ref = spectralData_ref.Sd;

% COMPENSATION ATTENUATION
acs_ref     = alpha_ref; % [dB/cm/MHz]
att_ref     = acs_ref*band/Np2dB; % vector [Np/cm]
att_ref_map = ones(size(Sp_sam)) .* reshape(att_ref, 1, 1, []); % 3D array as SLogRatio

SLogRatio = log(Sp_sam./Sd_sam) - ( log(Sp_ref./Sd_ref) - 4*zd_zp*1E2*att_ref_map ); % (rows, cols, freqChannels)
SLogRatio_vec = reshape(SLogRatio, [], size(SLogRatio, 3));   
mux_SLogRatio = 1./(abs(mean(SLogRatio_vec, 1, 'omitnan')) ./ std(SLogRatio_vec, 0, 1, 'omitnan') + 1E-5);
weightEstimators = rescale(mux_SLogRatio, 1, max(mux_SLogRatio));

%% PLOT FULL Spectral
% Get size of Sp
[m, n, pFull] = size(spectralData_sam.Sdfull);

nLines = 5;
lin_cen = round(n / 2); 
lat_range = max(1, lin_cen-fix(nLines/2)):min(n, lin_cen+fix(nLines/2)); 

Sd_matrix = squeeze(mean(spectralData_sam.Sdfull(:, lat_range, :), 2)); % Mean over 2nd dim (lateral)

freq_Pos   = (0:(pFull/2-1))*SAM.fs/pFull*1e-6; % freq_rangePos = size ( (0:(1/2-1/pFull))'*SAM.fs) ;
Sd_Pos     = Sd_matrix(:,1:pFull/2);

Sd_Pos_dB  = pow2db(Sd_Pos ./ max(Sd_Pos, [], 2));
%% Plot the FULL extracted spectra
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
hold off;
xlabel('Frequency [MHz]');
ylabel('NormMax [dB]');
title('Power Spectrum');
legend('Location', 'Best');

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

if plotBmode
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
%% OVERLAY ACS MAPS + BMODE 

estim_method = 'SLD';
fontSize     = 14;
figure,
set(gcf,'units','normalized','outerposition',[0 0.1 1 0.8]); box on;

tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
sgtitle(estim_method, 'FontSize', fontSize+2, 'FontWeight', 'bold');

%%%%%%%%%%%%%%%%%%%%%%%%%% RSLD (ACS) %%%%%%%%%%%%%%%%%%%%%%%%%%

units           = 1E3;
bmodeFull       = bmode_sam;
colorImg        = acs_rsld;
range_bmode     = [-60 0];
range_img       = [0 1.1];
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

    title({'RSLD $\alpha_s$:', ...
       [num2str(round(met_rsld(1), 3)), ' $\pm$ ', num2str(round(met_rsld(2), 3))]}, ...
      'Interpreter', 'latex');
    set(gca,'fontsize',fontSize)

%%%%%%%%%%%%%%%%%%%%%%%%%% RSLD (ACS) %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% TNV (ACS) %%%%%%%%%%%%%%%%%%%%%%%%%%

units           = 1E3;
bmodeFull       = bmode_sam;
colorImg        = acs_tnv;
range_bmode     = [-60 0];
range_img       = [0 1.1];
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

    title({'TNV $\alpha_s$:', ...
       [num2str(round(met_tnv(1), 3)), ' $\pm$ ', num2str(round(met_tnv(2), 3))]}, ...
      'Interpreter', 'latex');
    set(gca,'fontsize',fontSize)

%%%%%%%%%%%%%%%%%%%%%%%%%% TNV (ACS) %%%%%%%%%%%%%%%%%%%%%%%%%%