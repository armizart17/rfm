% Febrery Reference Frequency Method
% AMZ 

clear all, clc, close all;

Np2dB       = 20*log10(exp(1));
dB2Np       = 1/Np2dB;

%%
%% LOAD SAM
pathData = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\dataACS_kwave';

SAM = load(fullfile(pathData, 'ref0.mat'));

% B-MODE CHECK
bmode_sam = db(hilbert(SAM.rf));
bmode_sam = bmode_sam - max(bmode_sam(:));

figure,

% subplot(121), 
imagesc(SAM.x*1E3, SAM.z*1E3, bmode_sam), axis("image");
xlabel('Lateral [mm]'), ylabel('Depth [mm]');
title('SAM')
colormap('gray')

%%
pars.bw          = [3 9]; % [MHz]
pars.overlap     = 0.8;
pars.blocksize   = 20; % wavelengths
pars.z_roi       = [5 45]*1E-3; % [m] 
pars.x_roi       = [-17 17]*1E-3; % [m] 
pars.saran_layer = false;

pars.window_type = 3; %  (1) Hanning, (2) Tuckey, (3) Hamming, (4) Tchebychev

% Reading experiment settings parameters
bw              = pars.bw;
overlap         = pars.overlap;
blocksize_wv    = pars.blocksize;
z_ini           = pars.z_roi(1);
z_end           = pars.z_roi(2);
x_ini           = pars.x_roi(1);
x_end           = pars.x_roi(2);
saran_layer     = pars.saran_layer;

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

figure,
imagesc(SAM.x*1E3, SAM.z*1E3, bmode_sam), axis("image");
rectangle('Position', 1E3*[pars.x_roi(1) pars.z_roi(1) pars.x_roi(2)-pars.x_roi(1) pars.z_roi(2)-pars.z_roi(1)], ...
        'EdgeColor','r', 'LineWidth', 2, 'LineStyle','--'), hold off;
xlabel('Lateral [mm]'), ylabel('Depth [mm]');
title('SAM')
colormap('gray')

%% Spectral Ratio INIT
z            = SAM.z;
x            = SAM.x;
fs           = SAM.fs;

dx = x(2)-x(1);
dz = z(2)-z(1);
c0 = 1540; % [m/s]
lambda = c0/mean(bw)*1E-6;

%% Cropping and finding sample sizes
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
ratio_zx = 1; % @
wz = round(blocksize_wv*lambda*(1-overlap)/dz * ratio_zx); % Between windows
nz = 2*round(blocksize_wv*lambda/dz /2 * ratio_zx); % Window size

z0p = 1:wz:length(z)-nz;
z0d = z0p + nz/2;
z_ACS = z(z0p+ nz/2);
m  = length(z0p);

 
% Frequency samples
NFFT = 2^(nextpow2(nz/2)+1);
NFFT = 4096;

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

%% RFM

RSp = zeros(m,n,p);
RSp(:, :, 2:end) = Sp(:, :, 2:end) ./ Sp(:, :, 1:end-1);
% For the last slice, keep the ratio the same as the previous slice
RSp(:, :, 1) = RSp(:, :, 2); % Assuming the first slice ratios are 1 as there is no "i-1"


% Choose half lateral  
x_half = ceil(n/2);

% Find freqs MHz
freq1 = 3.5; freq2 = 6; freq3 = 8.5;
idx_f1 = find(band >= freq1, 1, 'first');
idx_f2 = find(band >= freq2, 1, 'first');
idx_f3 = find(band >= freq3, 1, 'first');

figure, 
sgtitle('Spectrum 1b')
plot(z_ACS*1E2, log(RSp(:, x_half, idx_f1)), 'DisplayName', sprintf('S_{%.2fMHz}/S_{%.2fMHz}', band(idx_f1), band(idx_f1-1)) ), hold on, grid on;
plot(z_ACS*1E2, log(RSp(:, x_half, idx_f2)), 'DisplayName', sprintf('S_{%.2fMHz}/S_{%.2fMHz}', band(idx_f2), band(idx_f2-1)) ), hold on
plot(z_ACS*1E2, log(RSp(:, x_half, idx_f3)), 'DisplayName', sprintf('S_{%.2fMHz}/S_{%.2fMHz}', band(idx_f3), band(idx_f3-1)) ), hold on

xlabel('Depth [cm]')
ylabel('Spectrum ln [dB]')
legend('Location','Best')
set(gca, 'FontSize', 14)

%% TWO DEPTHS
z0p = 1:wz:length(z)-nz;
z0d = z0p + nz/2;
Sp = zeros(m,n,p);
Sd = zeros(m,n,p);

windowing = window_choice(nz/2, window_type);
windowing = windowing*ones(1,nx);

for jj=1:n
    for ii=1:m
        xw = x0(jj) ;   % x window
        zp = z0p(ii);
        zd = z0d(ii);

        sub_block_p = samRef(zp:zp+nz/2-1,xw:xw+nx-1);
        sub_block_d = samRef(zd:zd+nz/2-1,xw:xw+nx-1);

        [tempSp,~] = spectra(sub_block_p,windowing,0,nz/2,NFFT);
        [tempSd,~] = spectra(sub_block_d,windowing,0,nz/2,NFFT);

        Sp(ii,jj,:) = tempSp(ind_f) ./ t_saran;
        Sd(ii,jj,:) = tempSd(ind_f) ./ t_saran;

    end
end

RSp = zeros(m,n,p);
RSd = zeros(m,n,p);

% Compute the ratio along the third dimension
RSp(:, :, 2:end) = Sp(:, :, 2:end) ./ Sp(:, :, 1:end-1);
% For the last slice, keep the ratio the same as the previous slice
RSp(:, :, 1) = 1; % Assuming the first slice ratios are 1 as there is no "i-1"

% Compute the ratio along the third dimension
RSd(:, :, 2:end) = Sd(:, :, 2:end) ./ Sd(:, :, 1:end-1);
% For the last slice, keep the ratio the same as the previous slice
RSd(:, :, 1) = 1; % Assuming the first slice ratios are 1 as there is no "i-1"

RSnorm = RSd ./ RSp;

%% SIMPLER WAY
df_MHz = min(diff(band));
delta_z_cm = zd_zp*1e-1; % [cm]

factorDiv = -4*df_MHz*delta_z_cm;
array3D = log(RSnorm) / factorDiv *Np2dB; %[dB/MHz/cm]
% HISTOGRAM

figure;
histogram(array3D(:,:,90), 20, 'Normalization', 'probability'); % 20 bins
xlabel('Value');
ylabel('Frequency');
title('Histogram of ACS Elements');
grid on;

figure, 
plot(squeeze(array3D(16,16,:)))
%% MATRIX WAY
b = log(RSd) - log(RSp);
b_vec = b(:);

df_MHz = min(diff(band));

zd_zp = (nz/2)*dz;   % zd_zp = 2*\delta_z = nz/2 =  %  [m]
delta_z_cm = zd_zp*1e-1; % [cm]

% System of eq
A1 = kron( -4*delta_z_cm*df_MHz*ones(p) , speye(m*n) );
A2 = kron( ones(p) , speye(m*n) );

A = [A1 A2];

x_opt = cgs(A'*A, A'*b(:));

Bn = x_opt(1:m*n*p);
Cn = x_opt(m*n*p+1:end);

BR_RFM = reshape(Bn*Np2dB,m,n,p);
CR_RFM = reshape(Cn,m,n,p);

% acs_map = reshape(BR_RFM*Np2dB,m,n);
% 
% figure, imagesc(acs_map), colorbar, colormap("jet");

%%
acs_map = mean(BR_RFM, 3);
figure, 
imagesc(x_ACS, z_ACS, acs_map), colorbar
