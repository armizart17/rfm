% Reference Frequency Method v2 March
% AMZ 

% clear all, clc, close all;

Np2dB       = 20*log10(exp(1));
dB2Np       = 1/Np2dB;
range_bmode = [-60 0];

%% LOAD SAM

% DATA NEW AMZ
pathData = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\dataACS_kwave';

%%%%%%%%%%%%%%% OLD FEBRUARY ** %%%%%%%%%%%%%%% 
% alpha_sam = 1; % ACS 0.5 0.7 1  
% alpha_ref = 0.7; % ACS 0.5 0.7 1 
% j_sam = 1.1;
% 
% folderDataSam = sprintf('sim_TUFFC25\\acs%.1f_pow%.1f', alpha_sam, j_sam);
% folderDataSam = strrep(folderDataSam, '.', 'p');
% rf_sam_name = strcat('rfref_', sprintf('%.3f', j_sam));
% rf_sam_name = strrep(rf_sam_name, '.', 'p');
% SAM = load(fullfile(pathData, folderDataSam, rf_sam_name + ".mat"));
% 
% folderDataRef = sprintf('sim_TUFFC25\\acs%.1f_pow%.1f', alpha_ref, j_sam);
% folderDataRef = strrep(folderDataRef, '.', 'p');
% rf_ref_name = strcat('rfref_', sprintf('%.3f', j_sam));
% rf_ref_name = strrep(rf_ref_name, '.', 'p');
% REF = load(fullfile(pathData, folderDataRef, rf_ref_name + ".mat"));
%%%%%%%%%%%%%%% OLD FEBRUARY ** %%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%% NEW MARCH %%%%%%%%%%%%%%%
alpha_sam = 0.7; % ACS 0.4 0.5 0.6 0.7 1
alpha_ref = 0.4; % ACS 0.4 0.5 0.6 0.7 1
folderDataSam = 'a_pow1p';
rf_sam_name = sprintf('rf1_a_%.2g', alpha_sam);
% rf_sam_name = strrep(rf_sam_name, '.', 'p');
SAM = load(fullfile(pathData, folderDataSam, rf_sam_name + ".mat"));

folderDataRef = 'a_pow1p';
rf_ref_name = sprintf('rf1_a_%.2g', alpha_ref);
% rf_ref_name = strrep(rf_ref_name, '.', 'p');
REF = load(fullfile(pathData, folderDataRef, rf_ref_name + ".mat"));
%%%%%%%%%%%%%%% NEW MARCH %%%%%%%%%%%%%%%

gt_acs  = alpha_sam;

% B-MODE CHECK
bmode_sam = db(hilbert(SAM.rf));
bmode_sam = bmode_sam - max(bmode_sam(:));

% figure,
% % subplot(121), 
% imagesc(SAM.x*1E3, SAM.z*1E3, bmode_sam, range_bmode), axis("image");
% xlabel('Lateral [mm]'), ylabel('Depth [mm]');
% cb = colorbar;
% cb.Label.String = 'dB'; % Add the label "dB"
% title('SAM')
% colormap('gray')

%% SPECTRAL PARAMETERS
pars.bw          = [3 9]; % [MHz]
pars.overlap     = 0.8;
pars.blocksize   = 10; % wavelengths
pars.z_roi       = [5 45]*1E-3; % [m] 
pars.x_roi       = [-17 17]*1E-3; % [m] 
pars.saran_layer = false;
pars.ratio_zx    = 1.25;
pars.window_type = 3; %  (1) Hanning, (2) Tuckey, (3) Hamming, (4) Tchebychev

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
if isfield(SAM, 'RF')
    rfdata_sam   = SAM.RF;
end

% figure,
% imagesc(SAM.x*1E3, SAM.z*1E3, bmode_sam), axis("image");
% rectangle('Position', 1E3*[pars.x_roi(1) pars.z_roi(1) pars.x_roi(2)-pars.x_roi(1) pars.z_roi(2)-pars.z_roi(1)], ...
%         'EdgeColor','r', 'LineWidth', 2, 'LineStyle','--'), hold off;
% xlabel('Lateral [mm]'), ylabel('Depth [mm]');
% title('SAM')
% colormap('gray')

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
% wx = round(blocksize(1)*(1-overlap)/dx);  % Between windows OLD
% nx = round(blocksize(1)/dx);  % Window size OLD
wx = round(blocksize_wv*lambda*(1-overlap)/dx);  % Between windows  
nx = round(blocksize_wv*lambda/ dx);

x0 = 1:wx:length(x)-nx;
x_ACS = x(1,x0+round(nx/2));
n  = length(x0);

% Axial samples
% wz = round(blocksize(2)*(1-overlap)/dz); % Between windows OLD
% nz = 2*round(blocksize(2)/dz /2); % Window size OLD
% ratio_zx = 1.5; % ^
% ratio_zx = 1; % @
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

for ff = 1:p 
    for jj = 1:n
        % Compute FPR
        FPR = log(squeeze(squeeze(RSp(:,jj,ff)))) / (4*df_MHz)*Np2dB;
        
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
freq1 = 3.5; freq2 = 6; freq3 = 8.5;
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

% % Plot the results
% figure,
% sgtitle('Frequency Power Ratio')
% set(gcf, 'Units', 'pixels', 'Position', [100, 100, 900, 700]); % [x, y, width, height] in pixels

% % Plot raw FPR data
% plot(z_ACS_cm, FPR_avg_f1, 'r-', 'DisplayName', sprintf('S_{%.2fMHz}/S_{%.2fMHz}', band(idx_f1), band(idx_f1-1)), 'LineWidth', 1.5), hold on
% plot(z_ACS_cm, FPR_avg_f2, 'b-', 'DisplayName', sprintf('S_{%.2fMHz}/S_{%.2fMHz}', band(idx_f2), band(idx_f2-1)), 'LineWidth', 1.5), hold on
% plot(z_ACS_cm, FPR_avg_f3, 'k-', 'DisplayName', sprintf('S_{%.2fMHz}/S_{%.2fMHz}', band(idx_f3), band(idx_f3-1)), 'LineWidth', 1.5), hold on
% 
% % Plot linear fits with slope values in the legend
% plot(z_ACS_cm, FPR_fit_avg_f1, 'r--', 'DisplayName', sprintf('Fit: S_{%.2fMHz}/S_{%.2fMHz} (Slope: %.3f)', band(idx_f1), band(idx_f1-1), slope_avg_f1), 'LineWidth', 2), hold on
% plot(z_ACS_cm, FPR_fit_avg_f2, 'b--', 'DisplayName', sprintf('Fit: S_{%.2fMHz}/S_{%.2fMHz} (Slope: %.3f)', band(idx_f2), band(idx_f2-1), slope_avg_f2), 'LineWidth', 2), hold on
% plot(z_ACS_cm, FPR_fit_avg_f3, 'k--', 'DisplayName', sprintf('Fit: S_{%.2fMHz}/S_{%.2fMHz} (Slope: %.3f)', band(idx_f3), band(idx_f3-1), slope_avg_f3), 'LineWidth', 2), hold on
% 
% grid on
% xlabel('Depth [cm]')
% ylabel('FPR [dB/MHz]')
% ylabel('FPR [dB\cdotMHz^{-1}]')
% legend('Location', 'Best')
% set(gca, 'FontSize', 14)


%% RESULTS

%% HISTOGRAM

% figure;
% histogram(ACS_RFM, 20, 'Normalization', 'probability'); % 20 bins
% xlabel('Value');
% ylabel('Frequency');
% title('Histogram of ACS Elements');
% grid on;

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


%% SLD METHOD

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
acs_ref     = alpha_ref; % [dB/cm/MHz]
att_ref     = acs_ref*band/Np2dB; % vector [Np/cm]
att_ref_map = ones(size(Sp_sam)) .* reshape(att_ref, 1, 1, []); % 3D array as SLogRatio

SLogRatio = log(Sp_sam./Sd_sam) - ( log(Sp_ref./Sd_ref) - 4*zd_zp*1E2*att_ref_map ); % (rows, cols, freqChannels)
clear('Sp_sam','Sd_sam', 'Sp_ref', 'Sd_ref', 'att_ref_map', 'att_ref');

%% SLD
[m, n, p] = size(SLogRatio);
A1 = kron( 4*zd_zp*1E2*band , speye(m*n) );
A2 = kron( ones(size(band)) , speye(m*n) );
[ACS_SLD, ~] = cgs_ACS([A1 A2], SLogRatio);

% Compute mean and standard deviation
mu      = mean(ACS_SLD(:));
sigma   = std(ACS_SLD(:));
cv      = sigma/mu;

% Generate the box plot
figure;
boxplot(ACS_SLD(:)); 
grid on;
yline(gt_acs, 'k--', 'LineWidth', 1.5); % Dashed line for GT ACS
% ylim([-5 5]);

% Set title with mean and standard deviation 
title(sprintf('SLD ACS: %.3f \\pm %.3f, %%CV=%.2f', mu, sigma, 100*cv));

% Axis labels
xlabel('ACS');
ylabel('Values');

% Add a text box showing GT ACS
text(0.05, 0.95, sprintf('GT ACS = %.2f', gt_acs), ...
    'Units', 'normalized', 'FontSize', 12, ...
    'Color', 'k', 'BackgroundColor', 'w', ...
    'EdgeColor', 'k', 'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'top');

%%
keyboard
%% OPTIMIZATION RSLD

% Optimization constants
tol = 1e-3;
mu1 = 10^4.5; 
mu2 = 10^4.5;

[m, n, p] = size(SLogRatio);
mask = ones(m,n,p);
A1 = kron( 4*zd_zp*1E2*band , speye(m*n) );
A2 = kron( ones(size(band)) , speye(m*n) );

[Bn,~] = AlterOpti_ADMM(A1,A2,SLogRatio(:),mu1,mu2,m,n,tol,mask(:));
ACS_RSLD = reshape(Bn*Np2dB,m,n);

% Compute mean and standard deviation
mu      = mean(ACS_RSLD(:));
sigma   = std(ACS_RSLD(:));
cv      = sigma/mu;

% % Generate the box plot
% figure;
% boxplot(ACS_RSLD(:)); 
% grid on;
% yline(gt_acs, 'k--', 'LineWidth', 1.5); % Dashed line for GT ACS
% % ylim([-5 5]);
% 
% % Set title with mean and standard deviation 
% title(sprintf('RSLD ACS: %.3f \\pm %.3f, %%CV=%.2f', mu, sigma, 100*cv));
% 
% % Axis labels
% xlabel('ACS');
% ylabel('Values');
% 
% % Add a text box showing GT ACS
% text(0.05, 0.95, sprintf('GT ACS = %.2f', gt_acs), ...
%     'Units', 'normalized', 'FontSize', 12, ...
%     'Color', 'k', 'BackgroundColor', 'w', ...
%     'EdgeColor', 'k', 'HorizontalAlignment', 'left', ...
%     'VerticalAlignment', 'top');


%% BOX PLOTS ALL

% Find the maximum number of rows among the arrays
maxLength  = max([length(ACS_RFM(:)), length(ACS_SLD(:)), length(ACS_RSLD(:))]);

% Pad each array with NaNs to match the maximum number of rows
ACS_RFM_padded  = [ACS_RFM(:); NaN(maxLength - length(ACS_RFM(:)), 1)];
ACS_SLD_padded  = [ACS_SLD(:); NaN(maxLength - length(ACS_SLD(:)), 1)];
ACS_RSLD_padded = [ACS_RSLD(:); NaN(maxLength - length(ACS_RSLD(:)), 1)];

% Concatenate horizontally
acs_data = [ACS_RFM_padded, ACS_SLD_padded, ACS_RSLD_padded];

labels = {'RFM', 'SLD', 'RSLD'};

% Create the box plot
figure;
boxplot(acs_data, labels);
yline(gt_acs, 'k--', 'LineWidth', 1.5); % Dashed line for GT ACS
xlabel('ACS Method');
ylabel('\alpha [dB\cdotcm^{-1}\cdotMHz^{-1}]');
title('Box Plot of ACS Maps');
% ylim([-2 2])
grid on;

%%
keyboard
%% DENOISING

% WINDOW FILTER 1D (AVERAGE AND MEDFILT)
windowSize = 9; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;

% TV
mu_TV = 3;
iter_TV = 1000;

% SAVE FPR terms WITH FIT UPDATE Feb W3

[m, n, p] = size(RSp);
% Make full ACS all freq
z_ACS_cm = z_ACS*1E2;
df_MHz = min(diff(band));

% Initialize arrays to store FPR, linear fits, and slopes
FPR_all = zeros(length(z_ACS_cm), n, p); % Store original FPR
FPR_fit_all = zeros(length(z_ACS_cm), n, p); % Store linear fit (lin_y)

slopes_all = zeros(p, n); % Store slope values
FPR_snr_all = zeros(p, n); % Store standard deviation of FPR 
ACS_RFM = zeros(p, n);

for ff = 1:p 
    for jj = 1:n
        % Compute FPR
        FPR = log(squeeze(squeeze(RSp(:,jj,ff)))) / (4*df_MHz)*Np2dB;

        % Denoising     
        % FPR   = filter(b, a, FPR); % average filter
        FPR = medfilt1(FPR, windowSize); % median filter
        % FPR = tv_denoise(FPR, mu_TV, iter_TV); % TV-1D
        
        % Perform linear fit
        [lin_slope , lin_intercept, lin_y, R_fit] = fit_linear(z_ACS_cm, FPR, 2); 
        
        % Store FPR, linear fit, and slope
        FPR_all(:, jj, ff) = FPR;
        FPR_fit_all(:, jj, ff) = lin_y;


        slopes_all(ff, jj) = lin_slope;
        
        % Store ACS value
        ACS_RFM(ff, jj) = -lin_slope; % -8.68*m/(4*delta_f) [dB/MHz/cm]

        % std
        FPR_snr_all(ff, jj) = mean(FPR)/std(FPR);
    end
end

% WITH SLOPE PRINT

% Find freqs MHz
freq1 = 3.5; freq2 = 6; freq3 = 8.5;
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
slope_avg_f1 = mean(slopes_all(idx_f1, pos_range), 2);
slope_avg_f2 = mean(slopes_all(idx_f2, pos_range), 2);
slope_avg_f3 = mean(slopes_all(idx_f3, pos_range), 2);

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
% ylabel('FPR [dB\cdotMHz^{-1}]')
legend('Location', 'Best')
ylim([-10 12]);
set(gca, 'FontSize', 14)

% BOX PLOT RESULTS RFM

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

%%

y = 0.1*randn(1, 100) + sin(linspace(0, 4*pi, 100)); % Generate noisy signal
lambda = 0.2;  % Regularization parameter
x_denoised = tv_denoise(y, lambda, 100);

figure;
plot(y, 'r'); hold on;
plot(x_denoised, 'b', 'LineWidth', 2);
legend('Noisy Signal', 'Denoised Signal');
title('Total Variation Denoising - 1D');

function x_denoised = tv_denoise(y, lambda, iter)
    % Implements Chambolle's algorithm for Total Variation Denoising (1D)
    % 
    % Inputs:
    %   y      - Noisy input signal (1D array)
    %   lambda - Regularization parameter (higher = smoother)
    %   iter   - Number of iterations (default: 50)
    %
    % Output:
    %   x_denoised - Denoised signal
    y = y(:);
    if nargin < 3
        iter = 50; % Default number of iterations
    end
    
    % Initialize variables
    x = y; 
    p = zeros(length(y), 1); % p should be one element smaller than x
    tau = 0.25;

    for k = 1:iter
        % Compute divergence (div_p) and update x
        div_p = [diff(p); 0];  
        x_new = y - lambda * div_p; 
        
        % Update p, making sure dimensions match
        div_x_nw = [diff(x_new); 0];
        p = p + tau * div_x_nw;
        p = p ./ max(1, abs(p)); 
        
        % Update x
        x = x_new;
    end
    
    x_denoised = x;
end

%%
pos_range = 1:n;
pos_range = floor(n/2);
% pos_range = 
% Find freqs MHz
freq1 = 3.5; freq2 = 6; freq3 = 8.5;
idx_f1 = find(band >= freq1, 1, 'first');
idx_f2 = find(band >= freq2, 1, 'first');
idx_f3 = find(band >= freq3, 1, 'first');

FPR_avg_f3 = mean(FPR_all(:, pos_range, idx_f3), 2); 

figure, 
plot(z_ACS_cm, FPR_avg_f3, 'k-', 'DisplayName', sprintf('S_{%.2fMHz}/S_{%.2fMHz}', band(idx_f3), band(idx_f3-1)), 'LineWidth', 1.5), hold on
grid on
xlabel('Depth [cm]')
ylabel('FPR [dB/MHz]')
% ylabel('FPR [dB\cdotMHz^{-1}]')
legend('Location', 'Best')
ylim([-10 5]);
set(gca, 'FontSize', 14)

figure,
list_mu_TV = logspace(log10(0.001), log10(7), 15);
iter_TV = 100;

for uu = 1:length(list_mu_TV)
    mu_TV = list_mu_TV(uu);
    FPR_avg_f3_den = tv_denoise(FPR, mu_TV, iter_TV); % TV-1D

    plot(z_ACS_cm, FPR_avg_f3_den, '-', 'DisplayName', sprintf('u = %.3g', mu_TV), 'LineWidth', 1.5), hold on
    hold on
    ylim([-19 5]);
    pause(0.5)
end
grid on
xlabel('Depth [cm]')
ylabel('FPR [dB/MHz]')
% ylabel('FPR [dB\cdotMHz^{-1}]')
legend('Location', 'Best')

set(gca, 'FontSize', 14)

%%
% Define log-spaced values for mu_TV
iter_TV = 100;
mu_TV_list = logspace(log10(0.05), log10(7), 15); 

% Initialize arrays to store ACS results for each mu_TV
ACS_RFM_all = cell(length(mu_TV_list), 1);

% Loop over different mu_TV values
for k = 1:length(mu_TV_list)
    mu_TV = mu_TV_list(k);
    
    % Initialize arrays to store FPR, linear fits, and slopes
    FPR_all = zeros(length(z_ACS_cm), n, p);
    FPR_fit_all = zeros(length(z_ACS_cm), n, p);
    slopes_all = zeros(n, p);

    for ff = 1:p 
        for jj = 1:n
            % Compute FPR
            FPR = log(squeeze(RSp(:,jj,ff))) / (4*df_MHz) * Np2dB;

            % Denoising using current mu_TV
            % FPR = tv_denoise(FPR, mu_TV, iter_TV);
            FPR = tv_denoise(FPR, mu_TV*fact_freq(ff), iter_TV); % TV-1D

            % Perform linear fit
            [lin_slope, lin_intercept, lin_y, R_fit] = fit_linear(z_ACS_cm, FPR, 2); 

            % Store results
            FPR_all(:, jj, ff) = FPR;
            FPR_fit_all(:, jj, ff) = lin_y;
            slopes_all(jj, ff) = lin_slope;
            
            % Store ACS value
            ACS_RFM(ff, jj) = -lin_slope; 
        end
    end
    
    % Store ACS results
    ACS_RFM_all{k} = ACS_RFM(:);
end

%% BOX PLOTS FOR EACH mu_TV
figure;
for k = 1:length(mu_TV_list)
    % Extract ACS values for current mu_TV
    acs_values = ACS_RFM_all{k};
    
    % Compute mean, standard deviation, and coefficient of variation
    mu = mean(acs_values);
    sigma = std(acs_values);
    cv = (sigma / mu) * 100;

    % Create subplot for each mu_TV
    subplot(2, 5, k); % Arrange in 2 rows, 5 columns
    boxplot(acs_values);
    grid on;
    yline(gt_acs, 'k--', 'LineWidth', 1.5); % Dashed line for GT ACS
    ylim([-5 5]);

    % Set title with mu_TV value
    title(sprintf('\\mu_{TV} = %.3f', mu_TV_list(k)));
    
    % Axis labels
    xlabel('ACS');
    ylabel('Values');

    % Add text annotation with mean, std, and CV
    text(0.7, 0.1, sprintf('%.3f \\pm %.3f, CV=%.2f%%', mu, sigma, cv), ...
        'Units', 'normalized', 'FontSize', 10, ...
        'Color', 'k', 'BackgroundColor', 'w', ...
        'EdgeColor', 'k', 'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'bottom');
end

% Global title for the entire figure
sgtitle('Effect of Different \mu_{TV} on ACS Estimation');


%% PLOT SNR vs SHIFT FREQ

figure, 
subplot(121)
imagesc(x_ACS*1e3, band, FPR_snr_all)
% axis("image")
xlabel('Lateral [mm]'), ylabel('Freq [MHz]')
title('\bfSNR FPR');

pos_range = 1:n;
% Compute the average slope over the selected positions
snr_FPR_pos_range = mean(FPR_snr_all(:, pos_range), 2);

subplot(122)
plot(band, snr_FPR_pos_range)
grid on;
xlabel('Freq [MHz]')
ylabel('SNR')
title('\bfSNR FPR');

%%
pos_range = 1:n;
% Compute the average slope over the selected positions
snr_FPR_pos_range = mean(FPR_snr_all(:, pos_range), 2);
snr_FPR_pos_range = rescale(snr_FPR_pos_range, 0.1, 1);
figure, 
plot(band, snr_FPR_pos_range)
grid on;
xlabel('bfFreq [MHz]')
ylabel('\bfSNR FPR ')
title('\bfSNR FPR')

figure, 
plot(band, pow2db(snr_FPR_pos_range))
grid on;
xlabel('bfFreq [MHz]')
ylabel('\bfSNR [dB]')
title('\bfNorm SNR FPR')

figure, 
semilogy(band, snr_FPR_pos_range)
grid on;
xlabel('bfFreq [MHz]')
ylabel('\bfSNR FPR [dB]')
title('\bfSNR FPR')

%% DENOISING

% WINDOW FILTER 1D (AVERAGE AND MEDFILT)
windowSize = 9; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;

% TV
mu_TV = 10;
fact_freq = snr_FPR_pos_range;
iter_TV = 1000;

% SAVE FPR terms WITH FIT UPDATE Feb W3

[m, n, p] = size(RSp);
% Make full ACS all freq
z_ACS_cm = z_ACS*1E2;
df_MHz = min(diff(band));

% Initialize arrays to store FPR, linear fits, and slopes
FPR_all = zeros(length(z_ACS_cm), n, p); % Store original FPR
FPR_fit_all = zeros(length(z_ACS_cm), n, p); % Store linear fit (lin_y)

slopes_all = zeros(p, n); % Store slope values
FPR_snr_all = zeros(p, n); % Store standard deviation of FPR 
ACS_RFM = zeros(p, n);

for ff = 1:p 
    for jj = 1:n
        % Compute FPR
        FPR = log(squeeze(squeeze(RSp(:,jj,ff)))) / (4*df_MHz)*Np2dB;

        % Denoising     
        % FPR   = filter(b, a, FPR); % average filter
        % FPR = medfilt1(FPR, windowSize); % median filter
        % FPR = tv_denoise(FPR, mu_TV, iter_TV); % TV-1D
        FPR = tv_denoise(FPR, mu_TV*fact_freq(ff), iter_TV); % TV-1D
        
        % Perform linear fit
        [lin_slope , lin_intercept, lin_y, R_fit] = fit_linear(z_ACS_cm, FPR, 2); 
        
        % Store FPR, linear fit, and slope
        FPR_all(:, jj, ff) = FPR;
        FPR_fit_all(:, jj, ff) = lin_y;

        slopes_all(ff, jj) = lin_slope;
        
        % Store ACS value
        ACS_RFM(ff, jj) = -lin_slope; % -8.68*m/(4*delta_f) [dB/MHz/cm]

        % std
        FPR_snr_all(ff, jj) = mean(FPR)/std(FPR);
    end
end

% WITH SLOPE PRINT

% Find freqs MHz
freq1 = 3.5; freq2 = 6; freq3 = 8.5;
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
slope_avg_f1 = mean(slopes_all(idx_f1, pos_range), 2);
slope_avg_f2 = mean(slopes_all(idx_f2, pos_range), 2);
slope_avg_f3 = mean(slopes_all(idx_f3, pos_range), 2);

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
% ylabel('FPR [dB\cdotMHz^{-1}]')
legend('Location', 'Best')
ylim([-10 12]);
set(gca, 'FontSize', 14)

% BOX PLOT RESULTS RFM

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