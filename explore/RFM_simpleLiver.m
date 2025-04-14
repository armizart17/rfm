% Script for C.Soto by EMZ
% Description : RFM Method to achieve ACS maps ACQ C52V 

% I already beamformed the data (check script acqLIM_curve.m)
% Just in case I made a copy of data in NAS2
% Q:\CAvendanoData\bf_M04_D02 or /mnt/nfs2/CAvendanoData/bf_M04_D02

% Modify iAcq for each case or loop all from the folder
% Inspect ROI manually or change it

% Modify Spectral Parameters: 
% - pars.blocksize      : 10 wl - 16 wl (sam depth of roi)
% - pars.blocksize_wv_r : 14 - 20 wl (ref depths from same roi)
% - pars.blocklines     : 8 - 10 lines 

% There's usually high variability, so only take median/mean ACS value to
% characterize tissue

% WHEN USING COLORBAR IN CURVE PLOT ACS MOVES (WATCH OUT)
%%
clear all;

manualroi   = true;

Np2dB       = 20*log10(exp(1));
dB2Np       = 1/Np2dB;
range_bmode = [-60 0];
range_acs   = [0 1.5];
fact_transp = 0.7;
fontSize    = 15;

% UTILS SIMPLE
mean2d      = @(x) mean(x(:));
std2d       = @(x) std(x(:));
cv2d        = @(x) 100*std(x(:))/mean(x(:));
median2d    = @(x) median(x(:));

calc2dStats = {@(x) mean(x(:)), @(x) std(x(:)), @(x) 100*std(x(:)) / mean(x(:)), @(x) median(x(:))};

%% SPECTRAL PARAMETERS
pars.bw          = [1.2 3.9]; % [MHz]
pars.overlap     = 0.8;
pars.blocksize   = 10; % wavelengths
pars.blocklines  = 8;
pars.saran_layer = false;
pars.ratio_zx    = 1;
pars.window_type = 3; %  (1) Hanning, (2) Tuckey, (3) Hamming, (4) Tchebychev

pars.blocksize_wv_r   = 14;

%% LOAD DATA

% Change *
iAcq = 1; 

% PATH PC (change if needed it)
pathData = 'Q:\CAvendanoData\bf_M04_D02';

% To save results (optional)
% figsDir = 'Q:\CAvendanoData\rfmV1'; 
% if ~exist(figsDir) mkdir(figsDir); end

acqDir   = dir(fullfile(pathData,'*.mat'));

samName  = acqDir(iAcq).name;
fileName = fullfile(pathData, samName);

SAM = load(fileName);

bmode_sam   = SAM.bMode;
SAM.x       = SAM.xr;
SAM.z       = SAM.zr;

caption = strrep(samName(1:end-4), '_', ' ');

%% ROI SELECTION
if manualroi 
figure('Units','centimeters', 'Position',[5 5 15 15]),
imagesc(SAM.x, SAM.z*1E3,bmode_sam,range_bmode);
colormap gray; clim(range_bmode);
hb2=colorbar; ylabel(hb2,'dB')
xlabel('Lateral [°]'), ylabel('Depth [mm]'); 
title(caption)

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

pars.x_roi     = [rect(1), rect(1)+rect(3)];      % [º]
pars.z_roi     = [rect(2), rect(2)+rect(4)]*1E-3; % [m]
end

% PRINT ROI
% fprintf('==========(ROI)===========\n')
fprintf('Depth [mm] = %s;\n', mat2str(round(1e3*pars.z_roi,2)));
fprintf('Lateral [°] = %s;\n', mat2str(round(pars.x_roi,2)));

%% RFM V1 

% Reading experiment settings parameters
bw              = pars.bw;
overlap         = pars.overlap;
blocksize_wv    = pars.blocksize;
blocklines      = pars.blocklines;
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

% Lateral samples (LINEAR)
% wx = round(blocksize_wv*lambda*(1-overlap)/dx);  % Between windows  
% nx = round(blocksize_wv*lambda/ dx);

% Lateral samples (CURVE)
wx = round(blocklines*(1-overlap));  % Between windows
nx = blocklines;                 % Window size

x0 = 1:wx:length(x)-nx;
x_ACS = x(x0+round(nx/2));
x_ACS = x_ACS(:); %Ensure col vector
n  = length(x0);


%%%%%%%%%%%%%%%%%%%%% Sample Depth %%%%%%%%%%%%%%%%%%%%%

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

band     = axis_f(ind_f)*1E-6; % [MHz]
bandFull = axis_f*1E-6; % [MHz]

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
% clear Sp_k
RSp_k(:,:, 1) = RSp_k(:,:, 2); % Assuming the first slice ratios are 1 as there is no "i-1"
RSp_k = log(RSp_k); % @@

%%%%%%%%%%%%%%%%%%%%% Reference Depth %%%%%%%%%%%%%%%%%%%%%

blocksize_wv_r = pars.blocksize_wv_r;

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
% clear Sp_r
RSp_r(:,:, 1) = RSp_r(:,:, 2); % Assuming the first slice ratios are 1 as there is no "i-1"
RSp_r = log(RSp_r); % @@

%% UFR Strategy

bw_ufr = pars.bw;
freqL  = bw_ufr(1); freqH = bw_ufr(2);
range  = bandFull >= freqL & bandFull <= freqH;

band_ufr    = bandFull(range);
p_ufr       = length(band_ufr);

RSp_k_ufr   = RSp_k(:,:,range);
RSp_r_ufr   = RSp_r(:,:,range);

%% NORMALIZED FPR

% Convert depth values to cm
z_ACS_cm = z_ACS * 1e2;      % Convert from meters to cm
z_ACS_r_cm = z_ACS_r * 1e2;  % Convert reference depths to cm

% Delta MHz 
df_MHz = band_ufr(2) - band_ufr(1);

% Preallocate cell arrays for storing x_temp and y_temp
x_temp_all = cell(m, n);
y_temp_all = cell(m, n);

tic;
for jj = 1:n  % Loop over lateral positions (x_j)
    for ii = 1:m  % Loop over depth positions (z_k)

        % Temporary storage for this location
        y_temp = nan(m_r, p_ufr);  % (Reference depth, Frequency) ** m_r
        x_temp = nan(m_r, p_ufr);  % (Reference depth, Frequency) ** m_r

        for r = 1:m_r  % Loop over reference depths
            % if (ii==1 && r==1)
            y_col = squeeze( ( (RSp_k_ufr(ii, jj, :)) - (RSp_r_ufr(r, jj, :)) ) /(4*df_MHz) *Np2dB ); % p_ufr x 1
     
            X = z_ACS_cm(ii) - z_ACS_r_cm(r);

            y_temp(r, :) = y_col; %**
            x_temp(r, :) = X; % **

        end
        x_temp_all{ii, jj} = x_temp;
        y_temp_all{ii, jj} = y_temp;

    end
end
t = toc;
fprintf('Loop way Elapsed time %.2f \n', t);
%%%%%%%%%%%%%%%%%%%% FAST WAY %%%%%%%%%%%%%%%%%%%%

% NOW ESTIMATION CGS RFM
a_rfm = zeros(m, n); 

for jj = 1:n  % Loop over lateral positions (x_j)
    for ii = 1:m  % Loop over depth positions (z_k)

        x_temp = x_temp_all{ii, jj};
        y_temp = y_temp_all{ii, jj};

        X_vec = x_temp(:);
        y_vec = y_temp(:);

        %%%%%%%%%%%%%%%% Option 1 %%%%%%%%%%%%%%%%
        % a_rfm(ii, jj) = -(X_vec' * X_vec) \ (X_vec' * y_vec);

        %%%%%%%%%%%%%%%% Option 2 %%%%%%%%%%%%%%%%
        A2_full    = kron(speye(p_ufr), x_temp(:, 1));
        A2t_full   = A2_full';

        slope_full = -(cgs(A2t_full*A2_full, A2t_full*y_vec));

        a_rfm(ii, jj) = mean(slope_full);
        % a_rfm(ii, jj) = median(slope_full);

    end
end

% Simple metrics
[m_a, s_a, cv_a, med_a]    = deal(calc2dStats{1}(a_rfm), calc2dStats{2}(a_rfm), calc2dStats{3}(a_rfm), calc2dStats{4}(a_rfm));

%% VISUALIZATION (RECT PLOT)

xFull = SAM.x;
zFull = SAM.z;
[X,Z] = meshgrid(xFull, zFull);

roi = X >= x_ACS(1) & X <= x_ACS(end) & Z >= z_ACS(1) & Z <= z_ACS(end);

figure('Units', 'pixels', 'Position', [100, 100, 1200, 600]); % [x, y, width, height]

tiledlayout(1, 2, 'TileSpacing','compact', 'Padding','compact')

t1 = nexttile();
imagesc(xFull,zFull*1e3, bmode_sam, range_bmode); % axis image;
title(caption)

% title([caption, ' bk', mat2str(round(pars.blocksize,2)), ...
%     'br', mat2str(round(blocksize_wv_r,2)), ...
%     'nL', mat2str(round(pars.blocklines,2))])
% ylim(yLimits)
hold on
contour(xFull,zFull,roi,1,'w--')
hold off
xlabel('Lateral [°]')
ylabel('Axial [mm]')
hBm = colorbar('Ticks',-60:20:0);
hBm.Label.String = 'dB';
% hBm.Location = 'westoutside';
set(gca,'fontsize',fontSize)

nexttile,
[~,hB,hColor] = imOverlayInterp(bmode_sam,a_rfm,range_bmode,range_acs,fact_transp,...
    x_ACS,z_ACS*1e3,roi,xFull,zFull*1e3);
title(sprintf('RFM: %.3f ± %.3f, CV=%.2f%%', m_a, s_a, cv_a));
% subtitle(['ACS = ',num2str(m_a,2),' dB/cm/MHz'])
axis normal
% ylim(yLimits)
hold on
contour(xFull,zFull*1e3,roi,1,'w--')
hold off
% axis off
hColor.Label.String = 'dB\cdotcm\cdotMHz^{-1}';
xlabel('Lateral [°]')

colormap(t1,'gray')
set(gca,'fontsize',fontSize)

%% VISUALIZATION (CURVE PLOT)

xPolar  = SAM.xp;
zPolar  = SAM.zp;
z0Polar = SAM.z0p; 
r0      = SAM.r0;
% r0 = 0.0401;

fontSize = 26;
roi = X >= x_ACS(3) & X <= x_ACS(end) & Z >= z_ACS(1) & Z <= z_ACS(end);

[TH_acs,R_acs] = meshgrid(-x_ACS*pi/180 + pi/2,z_ACS + r0);
[xPolarACS, zPolarACS] = pol2cart(TH_acs,R_acs);
zPolarACS = zPolarACS + z0Polar;


figure('Units','pixels', 'Position', [100, 100, 1200, 600]);

[ax1,~] = imOverlayPolar(bmode_sam,a_rfm,range_bmode,range_acs,0.7, ...
    xPolar,zPolar,xPolarACS,zPolarACS);

title(ax1, sprintf('RFM: %.3f ± %.3f, CV=%.2f%%', m_a, s_a, cv_a));
xlabel(ax1, 'Lateral [cm]');
ylabel(ax1, 'Axial [cm]');
set(ax1, 'FontSize', fontSize);
% ylim([-1 17])
% xlim([-12 12])
hold on;
% contour(xPolar*1e2, zPolar*1e2, roi, 1, 'w--');
% hb2=colorbar; ylabel(hb2,'dB\cdotcm\cdotMHz^{-1}', 'FontSize', fontSize)

% % Add label box "FL"
% annotation('textbox', [0.28, 0.785, 0.05, 0.10], 'String', 'FL', ...
%     'EdgeColor', 'w', 'BackgroundColor', 'k', 'Color', 'w','FontSize', fontSize+6, ...
%     'FontWeight', 'bold', 'HorizontalAlignment', 'center');
hold off;

% hb2=colorbar; ylabel(hb2,'dB\cdotcm\cdotMHz^{-1}')

set(ax1,'fontsize',fontSize)