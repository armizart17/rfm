% Reference Frequency Method v2 March
% AMZ 

clear all, 
% clc, close all;

Np2dB       = 20*log10(exp(1));
dB2Np       = 1/Np2dB;
range_bmode = [-60 0];
manualroi    = false;

dataVerasonics = true;
dataSonix      = false;
%% DATA VERASONICS

if dataVerasonics
pathData = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\SavedDataQUSPhantom\bf';
folderDataSam = '544'; numPhantomSam = '544'; alpha_sam = 0.53;
% folderDataSam = '261'; numPhantomSam = '261'; alpha_sam = 0.54;
samName = numPhantomSam + "_F_3";
end
%% DATA SONIX
if dataSonix
pathData = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\sonixtouch';
% folderDataSam = 'ID261V2'; numPhantomSam = 261; 
folderDataSam = 'ID544V2'; numPhantomSam = 544; alpha_sam = 0.53;

filesSam = dir(fullfile(pathData, folderDataSam,'*.mat'));
% samList = ["16-19-52","16-20-50","16-21-22", "16-21-22"];
samName = filesSam(1).name;
end
%%

SAM = load (fullfile(pathData, folderDataSam, samName));
SAM.rf = SAM.rf(:,:,1);

gt_acs  = alpha_sam;

% B-MODE CHECK
bmode_sam = db(hilbert(SAM.rf));
bmode_sam = bmode_sam - max(bmode_sam(:));


%% SPECTRAL PARAMETERS
pars.bw          = [3.5 8.5]; % [MHz]
pars.overlap     = 0.8;
pars.blocksize   = 7; % wavelengths
pars.saran_layer = false;
pars.ratio_zx    = 1.25;
pars.window_type = 3; %  (1) Hanning, (2) Tuckey, (3) Hamming, (4) Tchebychev

%% ROI SELECTION
if ~manualroi

    % pars.z_roi       = [5 35]*1E-3; % [m] 
    pars.z_roi       = [25 50]*1E-3; % [m] 
    pars.x_roi       = [-17 17]*1E-3; % [m] 
        
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


%% RFM V1 
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

figure,
imagesc(SAM.x*1E3, SAM.z*1E3, bmode_sam, range_bmode), axis("image");
rectangle('Position', 1E3*[pars.x_roi(1) pars.z_roi(1) pars.x_roi(2)-pars.x_roi(1) pars.z_roi(2)-pars.z_roi(1)], ...
        'EdgeColor','r', 'LineWidth', 2, 'LineStyle','--'), hold off;
xlabel('Lateral [mm]'), ylabel('Depth [mm]');
colormap gray; clim(range_bmode);
title('SAM')
colormap('gray')

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
x_ACS = x(1,x0+round(nx/2));
n  = length(x0);

% Axial samples
wz = round(blocksize_wv*lambda*(1-overlap)/dz * ratio_zx); % Between windows
nz = 2*round(blocksize_wv*lambda/dz /2 * ratio_zx); % Window size

z0p = 1:wz:length(z)-nz;
z0d = z0p + nz/2;
z_ACS = z(z0p+ nz/2);
m  = length(z0p);

% Frequency samples
% NFFT = 2^(nextpow2(nz/2)+1);
NFFT = 2^(nextpow2(nz/2)+2); %**

% axis_f = (0:(NFFT/2-1))'*fs/NFFT;  % [Hz] (so 0-fs/2),  it should be
axis_f = (0:NFFT-1)'/NFFT * fs;   % [Hz] freq axis as default because "spectra" function
freq_L = bw(1)*1E6; % [Hz] 
freq_H = bw(2)*1E6; % [Hz]

ind_f = axis_f >= freq_L & axis_f <= freq_H ;   

band  = axis_f(ind_f)*1E-6; % [MHz]
p = length(band);

zd_zp = (nz/2)*dz;   % zd_zp = 2*\delta_z = nz/2 =  %  [m]

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
Sp_ref = zeros(m,n,p,nSamples);
Sd_ref = zeros(m,n,p,nSamples);

% € 
SNR = zeros(m,n,nSamples);

% z0 = z0p
z0 = 1:wz:length(z)-nz;
z_ACS = z(z0+nz/2);
m  = length(z0);

% Spectrum Sample
Sp = zeros(m,n,p);
Sd = zeros(m,n,p);


samRef = SAM.rf;
samRef = samRef(ind_z,ind_x); % Cropping

    
f_MHz = axis_f((1:NFFT/2+1))*1e-6;
% figure, 
% set(gcf, 'Position', [0 0 1 1], 'Units', 'normalized');
font = 12;
for jj=1:n
    for ii=1:m
        xw = x0(jj) ;   % x window
        zw = z0(ii);

        sub_block = samRef(zw:zw+nz-1,xw:xw+nx-1);

        [tempSp,~] = spectra(sub_block,windowing,0,nz,NFFT); %@ nz/2 before

        Sp(ii,jj,:) = tempSp(ind_f);
     
    end
end

% RFM
RSp = zeros(m,n,p);
RSp(:, :, 2:end) = Sp(:, :, 2:end) ./ Sp(:, :, 1:end-1);
% For the last slice, keep the ratio the same as the previous slice
RSp(:, :, 1) = RSp(:, :, 2); % Assuming the first slice ratios are 1 as there is no "i-1"

%% SAVE FPR terms WITH FIT UPDATE Feb W3

% Make full ACS all freq
z_ACS_cm = z_ACS*1E2;
df_MHz = min(diff(band));

% Initialize arrays to store FPR, linear fits, and slopes
FPR_all = zeros(length(z_ACS_cm), n, p); % Store original FPR
FPR_fit_all = zeros(length(z_ACS_cm), n, p); % Store linear fit (lin_y)
slopes_all = zeros(n, p); % Store slope values

windowSize = 9;

for ff = 1:p 
    for jj = 1:n
        % Compute FPR
        FPR = log(squeeze(squeeze(RSp(:,jj,ff)))) / (4*df_MHz)*Np2dB;

        % FPR = medfilt1(FPR, windowSize); % median filter
        
        % Perform linear fit
        [lin_slope , lin_intercept, lin_y, R_fit] = fit_linear(z_ACS_cm, FPR, 2); 
        
        % Store FPR, linear fit, and slope
        FPR_all(:, jj, ff) = FPR;
        FPR_fit_all(:, jj, ff) = lin_y;
        slopes_all(jj, ff) = lin_slope;
        
        % Store ACS value
        ACS_RFM(ff, jj) = -lin_slope; % -8.68*m/(4*delta_f) [dB/MHz/cm]
    end
end

%% WITH SLOPE PRINT

% Find freqs MHz
freq1 = 3.75; freq2 = 6; freq3 = 8.25;
idx_f1 = find(band >= freq1, 1, 'first');
idx_f2 = find(band >= freq2, 1, 'first');
idx_f3 = find(band >= freq3, 1, 'first');

% Define the range of central horizontal positions to average
pos_range = 1:n; % Selected range of horizontal positions

% Compute the mean FPR and the mean of the linear fit over the selected range
FPR_avg_f1 = mean(FPR_all(:, pos_range, idx_f1), 2); 
FPR_avg_f2 = mean(FPR_all(:, pos_range, idx_f2), 2); 
FPR_avg_f3 = mean(FPR_all(:, pos_range, idx_f3), 2); 

FPR_fit_avg_f1 = mean(FPR_fit_all(:, pos_range, idx_f1), 2);
FPR_fit_avg_f2 = mean(FPR_fit_all(:, pos_range, idx_f2), 2);
FPR_fit_avg_f3 = mean(FPR_fit_all(:, pos_range, idx_f3), 2);

% Compute the average slope over the selected positions
slope_avg_f1 = mean(slopes_all(pos_range, idx_f1));
slope_avg_f2 = mean(slopes_all(pos_range, idx_f2));
slope_avg_f3 = mean(slopes_all(pos_range, idx_f3));

% Plot the results
figure,
sgtitle('Frequency Power Ratio')
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 900, 700]); % [x, y, width, height] in pixels

% Plot raw FPR data
plot(z_ACS_cm, FPR_avg_f1, 'r-', 'DisplayName', sprintf('S_{%.2fMHz}/S_{%.2fMHz}', band(idx_f1), band(idx_f1-1)), 'LineWidth', 1.5), hold on
plot(z_ACS_cm, FPR_avg_f2, 'b-', 'DisplayName', sprintf('S_{%.2fMHz}/S_{%.2fMHz}', band(idx_f2), band(idx_f2-1)), 'LineWidth', 1.5), hold on
plot(z_ACS_cm, FPR_avg_f3, 'k-', 'DisplayName', sprintf('S_{%.2fMHz}/S_{%.2fMHz}', band(idx_f3), band(idx_f3-1)), 'LineWidth', 1.5), hold on

% Plot linear fits with slope values in the legend
plot(z_ACS_cm, FPR_fit_avg_f1, 'r--', 'DisplayName', sprintf('Fit: S_{%.2fMHz}/S_{%.2fMHz} (Slope: %.3f)', band(idx_f1), band(idx_f1-1), slope_avg_f1), 'LineWidth', 2), hold on
plot(z_ACS_cm, FPR_fit_avg_f2, 'b--', 'DisplayName', sprintf('Fit: S_{%.2fMHz}/S_{%.2fMHz} (Slope: %.3f)', band(idx_f2), band(idx_f2-1), slope_avg_f2), 'LineWidth', 2), hold on
plot(z_ACS_cm, FPR_fit_avg_f3, 'k--', 'DisplayName', sprintf('Fit: S_{%.2fMHz}/S_{%.2fMHz} (Slope: %.3f)', band(idx_f3), band(idx_f3-1), slope_avg_f3), 'LineWidth', 2), hold on

grid on
xlabel('Depth [cm]')
ylabel('FPR [dB/MHz]')
ylabel('FPR [dB\cdotMHz^{-1}]')
legend('Location', 'Best')
set(gca, 'FontSize', 14)

%% BOX PLOT RESULTS RFM

% Compute mean and standard deviation
mu      = mean(ACS_RFM(:));
sigma   = std(ACS_RFM(:));
cv      = sigma/mu;

% Generate the box plot
figure;
boxplot(ACS_RFM(:)); 
grid on;
yline(gt_acs, 'k--', 'LineWidth', 1.5); % Dashed line for GT ACS
ylim([-5 5]);

% Set title with mean and standard deviation 
title(sprintf('RFM ACS: %.3f \\pm %.3f, %%CV=%.2f', mu, sigma, 100*cv));

% Axis labels
xlabel('ACS');
ylabel('Values');

% Add a text box showing GT ACS
text(0.05, 0.95, sprintf('GT ACS = %.2f', gt_acs), ...
    'Units', 'normalized', 'FontSize', 12, ...
    'Color', 'k', 'BackgroundColor', 'w', ...
    'EdgeColor', 'k', 'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'top');

