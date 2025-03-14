% SLD SIMULATION HOMOGENEOUS LIM 

%% GENERAL UTILS
% clear all, 
% clc, 
warning('off');
% close all;

addpath(genpath(pwd))

Np2dB       = 20*log10(exp(1));
dB2Np       = 1/Np2dB;
range_bmode = [-80 0];

plotBmode   = false;
manualroi   = true;
homoCase = true;

pathData = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\dataACS_kwave';

%% DATA FILES
if homoCase
alpha_sam = 0.7; % ACS 0.5 0.7 1 
j_sam = 1.1;

folderDataSam = sprintf('sim_TUFFC25\\acs%.1f_pow%.1f', alpha_sam, j_sam);
folderDataSam = strrep(folderDataSam, '.', 'p');

rf_sam_name = strcat('rfref_', sprintf('%.3f', j_sam));
rf_sam_name = strrep(rf_sam_name, '.', 'p');
rf_sam_name = strcat(rf_sam_name, '.mat');
SAM = load(fullfile(pathData, folderDataSam, rf_sam_name));

alpha_ref = 0.5; % ACS 0.5 0.7 1 
j_ref     = j_sam;

folderDataRef = sprintf('sim_TUFFC25\\acs%.1f_pow%.1f', alpha_ref, j_ref);
folderDataRef = strrep(folderDataRef, '.', 'p');

rf_ref_name = strcat('rfref_', sprintf('%.3f', j_sam));
rf_ref_name = strrep(rf_ref_name, '.', 'p');
rf_ref_name = strcat(rf_ref_name, '.mat');
REF = load(fullfile(pathData, folderDataRef, rf_ref_name));
REF.acs = alpha_ref;
else

    SAM = load(fullfile(pathData, 'sam0.mat'));
    REF = load(fullfile(pathData, 'ref0.mat')); % 0.4 dB/cm/MHz
    REF.acs = 0.4;

end

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

%% SLD STEP by STEP

pars.bw          = [3.5 8.5]; % [MHz]
pars.overlap     = 0.8;
pars.blocksize   = 20; % wavelengths
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
% weightEstimators = rescale(mux_SLogRatio, 1, max(mux_SLogRatio));

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
bestTau = 0.01;
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

    title({'TNV $\alpha_s$:', ...
       [num2str(round(met_tnv(1), 3)), ' $\pm$ ', num2str(round(met_tnv(2), 3))]}, ...
      'Interpreter', 'latex');
    set(gca,'fontsize',fontSize)

%%%%%%%%%%%%%%%%%%%%%%%%%% TNV (ACS) %%%%%%%%%%%%%%%%%%%%%%%%%%