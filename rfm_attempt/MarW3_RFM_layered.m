% Reference Frequency Method v2 March
% AMZ 

% clear all, clc, close all;

Np2dB       = 20*log10(exp(1));
dB2Np       = 1/Np2dB;
range_bmode = [-60 0];

mean2d = @(x) mean(x(:));
std2d = @(x) std(x(:));
cv2d = @(x) 100*std(x(:))/mean(x(:));

calc2dStats = {@(x) mean(x(:)), @(x) std(x(:)), @(x) 100 * std(x(:)) / mean(x(:))};

groundTruthTop    = [0.5,0.5,0.5,0.75,0.75];
groundTruthBottom = [1.0,1.0,1.0,0.75,0.75];

%% LOAD SAM

% DATA NEW AMZ
pathData = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\dataACS_kwave\layeredSM';

%%%%%%%%%%%%%%% NEW MARCH %%%%%%%%%%%%%%%
alpha_sam = 0.7; % ACS 0.4 0.5 0.6 0.7 1

folderDataSam = '24_04_04_layered';
% rf_sam_name = sprintf('rf1_a_%.2g', alpha_sam);
rf_sam_name = 'rf_layered1';
% rf_sam_name = strrep(rf_sam_name, '.', 'p');
% rf_sam_name = 'sam0';
SAM = load(fullfile(pathData, folderDataSam, rf_sam_name + ".mat"));

alpha_ref = 0.6; % ACS 0.6
folderDataRef = '24_04_25_ref';
rf_ref_name = 'rf_ref1';
% rf_ref_name = strrep(rf_ref_name, '.', 'p');
% rf_ref_name = 'ref0';
REF = load(fullfile(pathData, folderDataRef, rf_ref_name + ".mat"));
%%%%%%%%%%%%%%% NEW MARCH %%%%%%%%%%%%%%%

gt_acs  = alpha_sam;

% B-MODE CHECK
bmode_sam = db(hilbert(SAM.rf));
bmode_sam = bmode_sam - max(bmode_sam(:));

figure,
% subplot(121), 
imagesc(SAM.x*1E3, SAM.z*1E3, bmode_sam, range_bmode), axis("image");
xlabel('Lateral [mm]'), ylabel('Depth [mm]');
cb = colorbar;
cb.Label.String = 'dB'; % Add the label "dB"
title('SAM')
colormap('gray')

%% SPECTRAL PARAMETERS
pars.bw          = [3.5 8.5]; % [MHz]
pars.overlap     = 0.8;
pars.blocksize   = 12; % wavelengths
pars.z_roi       = [25 35]*1E-3; % [m] 
pars.x_roi       = [-14 14]*1E-3; % [m] 
pars.saran_layer = false;
pars.ratio_zx    = 1.25;
pars.window_type = 3; %  (1) Hanning, (2) Tuckey, (3) Hamming, (4) Tchebychev


blocksize_wv_r = 24;
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

%% Reference Depth

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

%% ATTEMPT RFM A_local UFR

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

        % Solve for local attenuation a(z_k, x_j) using least squares
        if ~isempty(y_vec)
            a_local_ufr(ii, jj) = ( (X_mat' * X_mat) \ (X_mat' * y_vec) )*Np2dB ;
        end
    end
end
t = toc;
fprintf('Loop way Elapsed time %.2f \n', t);

%% WAY MATRIX
% Precompute frequency and depth differences
% delta_f = diff(band_ufr, 1, 1); % (p_ufr-1, 1)
% delta_z = z_ACS_cm - z_ACS_r_cm'; % (m × m_r)
% 
% delta_f = reshape(delta_f, [1, 1, p_ufr - 1]);  % (1 × 1 × (p_ufr - 1))
% 
% % Define X for all depth and reference depth pairs
% X_mat = -4 * reshape(delta_f, [1, 1, p_ufr - 1]) .* delta_z; % (m × m_r × (p_ufr - 1))
% 
% 
% % Preallocate attenuation map
% a_local_ufr = zeros(m, n);  
% 
% tic;
% for jj = 1:n  % Loop over lateral positions (x_j)
% 
%     % Compute y for all (m × m_r × (p_ufr - 1)) at once
%     Y_mat = log(RSp_k_ufr(:, jj, 2:end)) - log(RSp_r_ufr(:, jj, 2:end)); % (m × m_r × (p-1))
% 
%     % Solve least squares for all depths at once
%     for ii = 1:m  % Loop over depth positions
%         y_vec = squeeze(Y_mat(ii, :, :));  % Extract row for depth ii (m_r × (p-1))
%         X_vec = squeeze(X_mat(ii, :, :));  % Extract row for depth ii (m_r × (p-1))
% 
%         if ~isempty(y_vec)
%             % Solve for local attenuation at (z_k, x_j)
%             a_local_ufr(ii, jj) = (X_vec(:)' * X_vec(:)) \ (X_vec(:)' * y_vec(:)) * 8.686; % Convert Np/cm/MHz to dB/cm/MHz
%         end
%     end
% 
% end
% t = toc;
% fprintf('Optimized Matrix Elapsed time %.2f seconds\n', t);
% 

% Apply constraints to stabilize attenuation estimates
% a_min = 0.1; % Lower bound for attenuation
% a_max = 2.0; % Upper bound for attenuation
% a_local = max(min(a_local, a_max), a_min); % Constrain within bounds


%% CLASSICAL SLD  METHOD

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
att_ref_map = reshape(att_ref, 1, 1, []); % 3D array as SLogRatio

SLogRatio = log(Sp_sam./Sd_sam) - ( log(Sp_ref./Sd_ref) - 4*zd_zp*1E2*att_ref_map ); % (rows, cols, freqChannels)
clear('Sp_sam','Sd_sam', 'Sp_ref', 'Sd_ref', 'att_ref_map', 'att_ref');

[m, n, ~] = size(SLogRatio);
A1 = kron( 4*zd_zp*1E2*band , speye(m*n) );
A2 = kron( ones(size(band)) , speye(m*n) );
[ACS_SLD, ~] = cgs_ACS([A1 A2], SLogRatio);

%% FIGURES COMPARISON RFM VS SLD

caxis_acs = [0 1.1];
fontSize = 14;

% ACS_SLD = ACS_RSLD;
figure, 
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1000, 600]); % [x, y, width, height]

% SLD
[m_a, s_a, cv_a] = deal(calc2dStats{1}(ACS_SLD), calc2dStats{2}(ACS_SLD), calc2dStats{3}(ACS_SLD));

subplot(121)
imagesc(x_ACS * 1e3, z_ACS * 1e3, ACS_SLD, caxis_acs); % Convert to mm
axis("image")
colorbar; colormap("turbo")
xlabel('Lateral [mm]');
ylabel('Depth [mm]');
hb2=colorbar; ylabel(hb2,'dB\cdotcm^{-1}\cdotMHz^{-1}')
% title('Local Attenuation Coefficient');
title(sprintf('SLD Local AC (GT= %.2f)\n%.3f $\\pm$ %.3f,  \\%%CV = %.2f', ...
               alpha_sam, m_a, s_a, cv_a), ...
      'Interpreter', 'latex');
set(gca,'fontsize',fontSize)

% RFM
[m_a, s_a, cv_a] = deal(calc2dStats{1}(a_local_ufr), calc2dStats{2}(a_local_ufr), calc2dStats{3}(a_local_ufr));

subplot(122)
imagesc(x_ACS * 1e3, z_ACS * 1e3, a_local_ufr, caxis_acs); % Convert to mm
axis("image")
colorbar; colormap("turbo")
xlabel('Lateral [mm]');
ylabel('Depth [mm]');
hb2=colorbar; ylabel(hb2,'dB\cdotcm^{-1}\cdotMHz^{-1}')
% title('Local Attenuation Coefficient');
title(sprintf('RFM Local AC (GT= %.2f)\n%.3f $\\pm$ %.3f,  \\%%CV = %.2f', ...
               alpha_sam, m_a, s_a, cv_a), ...
      'Interpreter', 'latex');
set(gca,'fontsize',fontSize)



%%
keyboard