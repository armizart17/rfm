% Reference Frequency Method Local Attenuation Layed Media by SM
% AMZ 

% clear all, clc, close all;

Np2dB       = 20*log10(exp(1));
dB2Np       = 1/Np2dB;
range_bmode = [-60 0];
range_acs   = [0.35 1.15];

mean2d = @(x) mean(x(:));
std2d = @(x) std(x(:));
cv2d = @(x) 100*std(x(:))/mean(x(:));

calc2dStats = {@(x) mean(x(:)), @(x) std(x(:)), @(x) 100 * std(x(:)) / mean(x(:))};

%% SPECTRAL PARAMETERS
pars.bw          = [3.5 8.5]; % [MHz]
pars.overlap     = 0.8;
pars.blocksize   = 12; % wavelengths
pars.z_roi       = [4 36]*1E-3; % [m] 
pars.x_roi       = [-15 15]*1E-3; % [m] 
pars.saran_layer = false;
pars.ratio_zx    = 1.25;
pars.window_type = 3; %  (1) Hanning, (2) Tuckey, (3) Hamming, (4) Tchebychev

blocksize_wv_r = 30;

%% DATAPATHS

dataDir = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\dataACS_kwave\layeredSM\24_04_04_layered';
refDir  = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\dataACS_kwave\layeredSM\24_04_25_ref';
resultsDir = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\dataACS_kwave\RFM\layeredSM';

resultsDir = strcat(resultsDir, '\ref_', num2str(blocksize_wv_r), 'wv');
tableName  = strcat('simuLay_ref', num2str(blocksize_wv_r), 'wv.xlsx');

if ~exist(resultsDir) mkdir(resultsDir); end

samFiles = dir([dataDir,'\rf*.mat']);
refFiles = dir([refDir,'\rf*.mat']);

% GT
groundTruthTop    = [0.5,0.5,0.5,0.75,0.75];
groundTruthBottom = [1.0,1.0,1.0,0.75,0.75];

%% For looping

for iAcq = 1:3

SAM = load(fullfile(dataDir,samFiles(iAcq).name));
REF = load(fullfile(refDir,refFiles(iAcq).name));

fprintf("Acquisition no. %i, id %s\n",iAcq,samFiles(iAcq).name);

%%%%%%%%%%%%%%%%%%%%%%%%%% SLD %%%%%%%%%%%%%%%%%%%%%%%%%%
% SAMPLE
spectralData_sam = calc_powerSpectra_prox_dis(SAM, pars);

Sp_sam = spectralData_sam.Sp;
Sd_sam = spectralData_sam.Sd;

% REFERENCE
spectralData_ref = calc_powerSpectra_prox_dis(REF, pars);

band   = spectralData_sam.band; % [MHz]
zd_zp  = spectralData_sam.zd_zp; % [m]
Sp_ref = spectralData_ref.Sp;
Sd_ref = spectralData_ref.Sd;

% COMPENSATION ATTENUATION
acs_ref     = REF.medium.alpha_coeff(1,1); % [dB/cm/MHz]
att_ref     = acs_ref*band/Np2dB; % vector [Np/cm]
att_ref_map = reshape(att_ref, 1, 1, []); % 3D array as SLogRatio [Np/cm]

SLogRatio = log(Sp_sam./Sd_sam) - ( log(Sp_ref./Sd_ref) - 4*zd_zp*1E2*att_ref_map ); % (rows, cols, freqChannels)
clear('Sp_sam','Sd_sam', 'Sp_ref', 'Sd_ref', 'att_ref_map', 'att_ref');

[m, n, ~] = size(SLogRatio);
A1 = kron( 4*zd_zp*1E2*band , speye(m*n) );
A2 = kron( ones(size(band)) , speye(m*n) );
[ACS_SLD, ~] = cgs_ACS([A1 A2], SLogRatio);
%%%%%%%%%%%%%%%%%%%%%%%%%% SLD %%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%% RFM %%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading experiment settings parameters
bw              = pars.bw;
overlap         = pars.overlap;
blocksize_wv    = pars.blocksize;
z_ini           = pars.z_roi(1);
z_end           = pars.z_roi(2);
x_ini           = pars.x_roi(1);
x_end           = pars.x_roi(2);
saran_layer     = pars.saran_layer;
ratio_zx        = pars.ratio_zx;

% Spectral Ratio INIT
z            = SAM.z;
x            = SAM.x;
fs           = SAM.fs;

% Default window
window_type     = 2; %  (1) Hanning, (2) Tukey, (3) Hamming, (4) Tchebychev
if isfield(pars, 'window_type')
    window_type = pars.window_type;
end

% Reading phantoms parameters
if isfield(SAM, 'rf')
    rfdata_sam   = SAM.rf;
end
if isfield(SAM, 'RF')
    rfdata_sam   = SAM.RF;
end

dx = x(2)-x(1);
dz = z(2)-z(1);
c0 = 1540; % [m/s]
lambda = c0/mean(bw)*1E-6;

% Cropping and finding sample sizes
% Limits for ACS estimation
ind_x = x_ini <= x & x <= x_end;
ind_z = z_ini <= z & z <= z_end;
x = x(ind_x);
z = z(ind_z);
rfdata_sam_roi = rfdata_sam(ind_z,ind_x);

% Lateral samples
wx = round(blocksize_wv*lambda*(1-overlap)/dx);  % Between windows  
nx = round(blocksize_wv*lambda/ dx);

x0 = 1:wx:length(x)-nx;
x_ACS = x(x0+round(nx/2));
x_ACS = x_ACS(:); %Ensure col vector
n  = length(x0);

% Axial samples
wz = round(blocksize_wv*lambda*(1-overlap)/dz * ratio_zx); % Between windows
nz = 2*round(blocksize_wv*lambda/dz /2 * ratio_zx); % Window size

z0 = 1:wz:length(z)-nz;
z_ACS = z(1, z0+round(nz/2));
z_ACS = z_ACS(:); %Ensure col vector
m  = length(z0);

% Frequency samples
% NFFT = 2^(nextpow2(nz/2)+1);
NFFT = 2^(nextpow2(nz/2)+2); %**

axis_f = (0:NFFT-1)'/NFFT * fs;   % [Hz] freq axis as default because "spectra" function
freq_L = bw(1)*1E6; % [Hz] 
freq_H = bw(2)*1E6; % [Hz]

ind_f = axis_f >= freq_L & axis_f <= freq_H ;   

band  = axis_f(ind_f)*1E-6; % [MHz]
bandFull  = axis_f*1E-6; % [MHz]

p = length(band);


fprintf('\nBandwith: %.2f - %.2f MHz\n',freq_L*1e-6,freq_H*1e-6)
fprintf('Blocksize in wavelengths: %i\n',blocksize_wv)
fprintf('Blocksize x: %.2f mm, z: %.2f mm\n',nx*dx*1e3,nz*dz*1e3)
fprintf('Blocksize in pixels nx: %i, nz: %i\n',nx,nz);
fprintf('Region of interest columns: %i, rows: %i\n\n',m,n);

% Saran Layer compensation 
t_saran = 1;
if (saran_layer)
    t_saran = transmission_saran(band);
end

% Windows for spectrum

windowing = window_choice(nz, window_type); %@ nz/2 before
windowing = windowing*ones(1,nx);

nSamples = size(rfdata_sam_roi,3);

% € 
SNR = zeros(m,n,nSamples);

% Spectrum Sample
Sp_k = zeros(m,n,NFFT);

samRef = rfdata_sam_roi;

% Spectrum
for jj=1:n
    for ii=1:m
        xw = x0(jj) ;   % x window
        zw = z0(ii);

        sub_block = samRef(zw:zw+nz-1,xw:xw+nx-1);

        [tempSp,~] = spectra(sub_block,windowing,0,nz,NFFT); %@ nz/2 before

        Sp_k(ii,jj,:) = tempSp();
     
    end
end

% RFM
RSp_k = zeros(m,n,NFFT);
RSp_k(:,:, 2:end) = Sp_k(:,:, 2:end) ./ Sp_k(:,:, 1:end-1);
% For the first slice, keep the ratio the same as the first slice
RSp_k(:,:, 1) = RSp_k(:,:, 2); % Assuming the first slice ratios are 1 as there is no "i-1"

%%% Reference Depth
% blocksize_wv_r = blocksize_wv *2;
blocksize_wv_r = blocksize_wv_r;

wz_r = round(blocksize_wv_r*lambda*(1-overlap)/dz * ratio_zx); % Between windows
nz_r = 2*round(blocksize_wv_r*lambda/dz /2 * ratio_zx); % Window size

z0_r = 1:wz_r:length(z)-nz_r;
z_ACS_r = z(z0_r+ nz/2);
z_ACS_r = z_ACS_r(:); %Ensure col vector
m_r  = length(z0_r);

% Frequency samples
% NFFT = 2^(nextpow2(nz/2)+1);
NFFT = 2^(nextpow2(nz/2)+2); %**

axis_f = (0:NFFT-1)'/NFFT * fs;   % [Hz] freq axis as default because "spectra" function
freq_L = bw(1)*1E6; % [Hz] 
freq_H = bw(2)*1E6; % [Hz]

ind_f = axis_f >= freq_L & axis_f <= freq_H ;   

band  = axis_f(ind_f)*1E-6; % [MHz]
bandFull  = axis_f*1E-6; % [MHz]

fprintf('\nBandwith: %.2f - %.2f MHz\n',freq_L*1e-6,freq_H*1e-6)
fprintf('Blocksize in wavelengths: %i\n',blocksize_wv)
fprintf('Blocksize x: %.2f mm, z: %.2f mm\n',nx*dx*1e3,nz*dz*1e3)
fprintf('Blocksize in pixels nx: %i, nz: %i\n',nx,nz);
fprintf('Region of interest columns: %i, rows: %i\n\n',m,n);

% Windows for spectrum
windowing = window_choice(nz, window_type); %@ nz/2 before
windowing = windowing*ones(1,nx);

% Spectrum Sample
Sp_r = zeros(m_r,n,NFFT);

samRef = rfdata_sam_roi;

% figure, 
% set(gcf, 'Position', [0 0 1 1], 'Units', 'normalized');
font = 12;
for jj=1:n
    for ii=1:m_r
        xw = x0(jj) ;   % x window
        zw = z0_r(ii);

        sub_block = samRef(zw:zw+nz-1,xw:xw+nx-1);

        [tempSp,~] = spectra(sub_block,windowing,0,nz,NFFT); %@ nz/2 before

        Sp_r(ii,jj,:) = tempSp();
     
    end
end

% RFM
RSp_r = zeros(m_r,n,NFFT);
RSp_r(:,:, 2:end) = Sp_r(:,:, 2:end) ./ Sp_r(:,:, 1:end-1);
% For the first slice, keep the ratio the same as the first slice
RSp_r(:,:, 1) = RSp_r(:,:, 2); % Assuming the first slice ratios are 1 as there is no "i-1"

%%% ATTEMPT RFM A_local UFR

% UFR strategy
freqL = 3.5; freqH = 8.5;
range = bandFull >= freqL & bandFull <= freqH;

RSp_k_ufr = RSp_k(:,:,range);
RSp_r_ufr = RSp_r(:,:,range);

band_ufr = bandFull(range);
p_ufr = length(band_ufr);

% Convert depth values to cm
z_ACS_cm = z_ACS * 1e2;      % Convert from meters to cm
z_ACS_r_cm = z_ACS_r * 1e2;  % Convert reference depths to cm

a_local_ufr = zeros(m, n); % Preallocate local attenuation matrix (depth x lateral)
% WAY LOOP
tic;
for jj = 1:n  % Loop over lateral positions (x_j)
    for ii = 1:m  % Loop over depth positions (z_k)
        
        y_vec = []; % Initialize y vector for this location
        X_mat = []; % Initialize X matrix for this location
        
        for r = 1:m_r  % Loop over reference depths
            for i = 2:p_ufr  % Loop over frequency bins
                % Compute y = log(RSnorm) at this depth & lateral position
                y = log(RSp_k_ufr(ii, jj, i)) - log(RSp_r_ufr(r, jj, i));
                
                % Define X = -4 * (fi - fi-1) * (zk - zr)
                X = -4 * (band_ufr(i) - band_ufr(i-1)) * (z_ACS_cm(ii) - z_ACS_r_cm(r));

                % Store values for least squares regression
                y_vec = [y_vec; y(:)];
                X_mat = [X_mat; X];
            end
        end

        % Solve for local attenuation a(z_k, BSWIFTx_j) using least squares
        if ~isempty(y_vec)
            a_local_ufr(ii, jj) = ( (X_mat' * X_mat) \ (X_mat' * y_vec) )*Np2dB ;
        end
    end
end
t = toc;
fprintf('Loop way Elapsed time %.2f \n', t);

ACS_RFM = a_local_ufr;
%%%%%%%%%%%%%%%%%%%%%%%%%% RFM %%%%%%%%%%%%%%%%%%%%%%%%%%

%% Setting Up

z_ACS = spectralData_sam.depth;
x_ACS = spectralData_sam.lateral;
z     = spectralData_sam.z_roi;
x     = spectralData_sam.x_roi;
Bmoderoi = my_RF2Bmode(spectralData_sam.rf_roi);

if iAcq == 4
    z = z-0.02;
end
% Creating reference
[X,Z]   = meshgrid(x_ACS, z_ACS);
[Xq,Zq] = meshgrid(x, z);

attIdeal        = ones(size(Z));
attIdeal(Z<=2E-2)  = groundTruthTop(iAcq);
attIdeal(Z>2E-2)   = groundTruthBottom(iAcq);

% 0.1 cm interface
top     = Zq < 1.9E-2; 
bottom  = Zq > 2.1E-2;
%% SLD METRICS

axialSLD = mean(ACS_SLD,2);

AttInterp       = interp2(X,Z,ACS_SLD,Xq,Zq);
m_sld.meanTop       = mean(AttInterp(top),"omitnan");
m_sld.stdTop        = std(AttInterp(top),"omitnan");
m_sld.meanBottom    = mean(AttInterp(bottom),"omitnan");
m_sld.stdBottom     = std(AttInterp(bottom),"omitnan");
m_sld.maeTop        = mean( abs(AttInterp(top) - groundTruthTop(iAcq))/groundTruthTop(iAcq),"omitnan");
m_sld.maeBottom     = mean( abs(AttInterp(bottom) - groundTruthBottom(iAcq))/groundTruthBottom(iAcq),"omitnan");
m_sld.nrmseTop      = sqrt(mean( (AttInterp(top) - groundTruthTop(iAcq)).^2,"omitnan")) / m_sld.meanTop;
m_sld.nrmseBottom   = sqrt(mean( (AttInterp(bottom) - groundTruthBottom(iAcq)).^2,"omitnan")) /m_sld.meanBottom;
m_sld.cnr           = abs(m_sld.meanBottom - m_sld.meanTop)/sqrt(m_sld.stdTop^2 + m_sld.stdBottom^2);
m_sld.method        = 'SLD';
m_sld.sample        = iAcq;
MetricsSLD(iAcq) = m_sld;

% clear m_sld

%% RFM METRICS

axialRFM        = mean(ACS_RFM,2);

AttInterp       = interp2(X,Z,ACS_RFM,Xq,Zq);
m_rfm.meanTop       = mean(AttInterp(top),"omitnan");
m_rfm.stdTop        = std(AttInterp(top),"omitnan");
m_rfm.meanBottom    = mean(AttInterp(bottom),"omitnan");
m_rfm.stdBottom     = std(AttInterp(bottom),"omitnan");
m_rfm.maeTop        = mean( abs(AttInterp(top) - groundTruthTop(iAcq))/groundTruthTop(iAcq),"omitnan");
m_rfm.maeBottom     = mean( abs(AttInterp(bottom) - groundTruthBottom(iAcq))/groundTruthBottom(iAcq),"omitnan");
m_rfm.nrmseTop      = sqrt(mean( (AttInterp(top) - groundTruthTop(iAcq)).^2,"omitnan")) / m_rfm.meanTop;
m_rfm.nrmseBottom   = sqrt(mean( (AttInterp(bottom) - groundTruthBottom(iAcq)).^2,"omitnan")) /m_rfm.meanBottom;
m_rfm.cnr           = abs(m_rfm.meanBottom - m_rfm.meanTop)/sqrt(m_rfm.stdTop^2 + m_rfm.stdBottom^2);
m_rfm.method = 'RFM';
m_rfm.sample = iAcq;
MetricsRFM(iAcq) = m_rfm;

% clear m_rfm

%% METRICS VEMZ

attIdeal_big = bigImg(attIdeal, Bmoderoi);

ACS_SLD_big = bigImg(ACS_SLD, Bmoderoi);

top2         = mask_rect(x,z,pars.x_roi(1),pars.x_roi(2),pars.z_roi(1),1.9E-2, ACS_SLD_big);
bottom2      = mask_rect(x,z,pars.x_roi(1),pars.x_roi(2),2.1E-2,pars.z_roi(2), ACS_SLD_big);

ACS_RFM_big = bigImg(ACS_RFM, Bmoderoi);


m_sld2 = get_metrics_gt(ACS_SLD_big, top2, bottom2, 'SLD', groundTruthTop(iAcq), groundTruthBottom(iAcq));
m_rfm2 = get_metrics_gt(ACS_RFM_big, top2, bottom2, 'RFM', groundTruthTop(iAcq), groundTruthBottom(iAcq));

%% Plotting BIG IMAGE VEMZ
%% Plotting BIG IMAGE VEMZ
cm = 1e2;
fontSize = 14;

% figure('Units','centimeters', 'Position',[5 5 30 10]);
figure,
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1550, 500]); % [x, y, width, height]
tiledlayout(1,4, "Padding","tight", 'TileSpacing','compact');

% B-mode Image
t1 = nexttile;
imagesc(x*cm, z*cm, Bmoderoi, range_bmode)
axis image
colormap(t1, gray)
title('B-mode')
ylabel('Depth [cm]')
xlabel('Lateral [cm]')
set(gca,'fontsize',fontSize)

% Ideal Image
t2 = nexttile;
imagesc(x*cm, z*cm, attIdeal, range_acs)
xlabel('Lateral [cm]')
colormap(t2, jet)
axis image
title('Ideal')
set(gca,'fontsize',fontSize)

% SLD Image with Formatted Title
t3 = nexttile; 
imagesc(x*cm, z*cm, ACS_SLD_big, range_acs)
colormap(t3, jet)
axis image
xlabel('Lateral [cm]')
set(gca,'fontsize',fontSize)

% Generate title with LaTeX formatting
title(t3, {['SLD:'], ...
    sprintf('Top = %.2f $\\pm$ %.2f', m_sld.meanTop, m_sld.stdTop), ...
    sprintf('Bot = %.2f $\\pm$ %.2f', m_sld.meanBottom, m_sld.stdBottom)}, ...
    'Interpreter', 'latex')

% RFM Image with Formatted Title
t4 = nexttile; 
imagesc(x*cm, z*cm, ACS_RFM_big, range_acs)
colormap(t4, jet)
axis image
xlabel('Lateral [cm]')
set(gca,'fontsize',fontSize)

% Generate title with LaTeX formatting
title(t4, {['RFM:'], ...
    sprintf('Top = %.2f $\\pm$ %.2f', m_rfm.meanTop, m_rfm.stdTop), ...
    sprintf('Bot = %.2f $\\pm$ %.2f', m_rfm.meanBottom, m_rfm.stdBottom)}, ...
    'Interpreter', 'latex')

% Add Colorbar
c = colorbar;
c.Label.String = 'ACS [dB\cdotcm^{-1}\cdotMHz^{-1}]';



%% Plotting
% cm = 1e2;
% fontSize = 14;
% 
% % figure('Units','centimeters', 'Position',[5 5 30 10]);
% figure,
% set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1300, 375]); % [x, y, width, height]
% tiledlayout(1,4, "Padding","tight", 'TileSpacing','compact');
% 
% t1 = nexttile;
% imagesc(x*cm,z*cm,Bmoderoi,range_bmode)
% axis image
% % xlim([x_ACS(1) x_ACS(end)]*cm),
% % ylim([z_ACS(1) z_ACS(end)]*cm),
% colormap(t1,gray)
% % c = colorbar(t1, 'westoutside');
% % c.Label.String = '[dB]';
% title('B-mode')
% ylabel('Depth [cm]')
% xlabel('Lateral [cm]')
% set(gca,'fontsize',fontSize)
% % hold on 
% % xline(2.05, 'w--', 'LineWidth',1.5)
% % hold off
% 
% t2 = nexttile;
% imagesc(x_ACS*cm,z_ACS*cm,attIdeal,range_acs)
% xlabel('Lateral [cm]'), % ylabel('Axial [cm]')
% colormap(t2,turbo)
% axis image
% title('Ideal')
% set(gca,'fontsize',fontSize)
% 
% t1 = nexttile; 
% imagesc(x_ACS*cm,z_ACS*cm,ACS_SLD, range_acs)
% colormap(t1,turbo)
% axis image
% title('SLD')
% % ylabel('Axial [cm]')
% xlabel('Lateral [cm]')
% set(gca,'fontsize',fontSize)
% % hold on 
% % yline(2, 'w--', 'LineWidth',1.5)
% % hold off
% 
% t4 = nexttile; 
% imagesc(x_ACS*cm,z_ACS*cm,ACS_RFM, range_acs)
% colormap(t4,turbo)
% axis image
% title('RFM')
% c = colorbar;
% c.Label.String = 'ACS [dB\cdotcm^{-1}\cdotMHz^{-1}]';
% % ylabel('Axial [cm]')
% xlabel('Lateral [cm]')
% set(gca,'fontsize',fontSize)
% % hold on 
% % yline(2, 'w--', 'LineWidth',1.5)
% % hold off
% 
% % fontsize(gcf,8,'points')

%% Lateral and axial profiles
lineColors = [0.635 0.078 0.184; 0.466 0.674 0.188; 0.301 0.745 0.933];
lw = 2;

% figure('Units','centimeters', 'Position',[5 5 14 8])
figure, 
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 700, 500]); % [x, y, width, height]
% tiledlayout(1,2, 'TileSpacing','compact', 'Padding','compact')
% nexttile,
plot(z_ACS*cm, axialSLD, '-.', 'LineWidth',lw, 'Color', 'r' ), hold on
plot(z_ACS*cm, axialRFM, '-.','LineWidth',lw, 'Color', 'b' ), 
plot(z_ACS*cm,mean(attIdeal,2), '--', 'LineWidth',lw, 'Color', 'k')
hold off
grid on
ylim([-1 2.5])
xlim(cm*[z_ACS(1) z_ACS(end)])
title('\bf Axial Profile')
%title('Axial profile')
xlabel('Depth [cm]')
ylabel('ACS [dB\cdotcm^{-1}\cdotMHz^{-1}]')
% if iAcq==1
    legend({'SLD','RFM'}, 'Location','Best') ;
% end
set(gca,'fontsize',fontSize+2)
%%
save_all_figures_to_directory(resultsDir,['Fig_simuLay',num2str(iAcq)]);
close all
pause(0.05);

end

%%
results1 = struct2table(MetricsSLD);
results2 = struct2table(MetricsRFM);

T = [results1; results2];
writetable(T, fullfile(resultsDir,tableName), 'WriteRowNames',true);

%%
