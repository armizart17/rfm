% Average 2D-filters
 clear; %close all; clc;
% 
mask = ones(9,9)/81;

figure; imagesc(mask);
colorbar; title('Mask')

Rmask = conv2(mask, mask);

figure; imagesc(Rmask); 
colorbar; title('Autocorrelation')

figue; surf(Rmask)
title('Xcorr')


NFFT = 2^10;
SRmask = fftshift(abs(fft2(Rmask,NFFT,NFFT)));

figure; imagesc(10*log10(SRmask)); 
colorbar; title('FFT(Autocorrelation)')

figure; surf(10*log10(SRmask))
title('Xcorr')


% spatial frequency 
dx = 0.075e-3 /2;
fs_space = 1540/dx;
fs_space / 2
freq = ((-NFFT/2):(NFFT/2-1))/NFFT*fs_space*1e-6;
figure; plot(10*log10(SRmask(:,512)))
figure; plot(freq, 10*log10(SRmask(:,512)))

return

% Gaussian 2D-filters
load("maskgauss.mat")

% mask6: gaussian, size=5, sgm=1
% mask7: gaussian, size=9, sgm=1
% mask8: gaussian, size=5, sgm=2
% mask9: gaussian, size=9, sgm=2

% mask=mask6;
% 
% figure; imagesc(mask);
% colorbar; title('Mask')
% 
% Rmask = conv2(mask, mask);
% 
% figure; imagesc(Rmask); 
% colorbar; title('Autocorrelation')
% 
% NFFT = 2^10;
% Smask = fftshift(abs(fft2(Rmask,NFFT,NFFT)));
% 
% figure; imagesc(10*log10(Smask)); 
% colorbar; title('FFT(Autocorrelation)')
% 

mask = mask7;

figure; imagesc(mask);
colorbar; title('Mask')

Rmask = conv2(mask, mask);

figure; imagesc(Rmask); 
colorbar; title('Autocorrelation')

NFFT = 2^10;
SRmask = fftshift(abs(fft2(Rmask,NFFT,NFFT)));

figure; imagesc(20*log10(SRmask)); 
colorbar; title('FFT(Autocorrelation)')

% spatial frequency 
dx = 0.3e-3 /16; % sigma is one dx.
ft_space = (1540/2)/(dx); % c0/2 / dx
ft_space / 2
freq = ((-NFFT/2):(NFFT/2-1))/NFFT*ft_space*1e-6;
figure; plot(10*log10(SRmask(:,NFFT/2)))
figure; plot(freq, 10*log10(SRmask(:,NFFT/2)),'k')



log_bsc = log(SRmask(:,NFFT/2));
fL = 3;
fH = 7.5;
f0 = mean([fL fH]);

[~,idx_fL] = min(abs(freq-fL));
[~,idx_fH] = min(abs(freq-fH));

new_freq = freq(idx_fL:idx_fH);
new_log_bsc = log_bsc(idx_fL:idx_fH);
slope = polyfit(new_freq,new_log_bsc,1)
n = f0 * slope(1); % Equal to n/f0
n

%%

mask=mask8;

figure; imagesc(mask);
colorbar; title('Mask')

Rmask = conv2(mask, mask);

figure; imagesc(Rmask); 
colorbar; title('Autocorrelation')

NFFT = 2^10;
Smask = fftshift(abs(fft2(Rmask,NFFT,NFFT)));

figure; imagesc(10*log10(Smask)); 
colorbar; title('FFT(Autocorrelation)')

mask=mask9;

figure; imagesc(mask);
colorbar; title('Mask')

Rmask = conv2(mask, mask);

figure; imagesc(Rmask); 
colorbar; title('Autocorrelation')

%NFFT = 2^10;
SRmask = fftshift(abs(fft2(Rmask,NFFT,NFFT)));

figure; imagesc(20*log10(SRmask)); 
colorbar; title('FFT(Autocorrelation)')

figure; plot(10*log10(SRmask(:,512)))

% Circular 2D-filters

load('maskcircle.mat')
mask=mask11;

figure; imagesc(mask);
colorbar; title('Mask')

Rmask = conv2(mask, mask);

figure; imagesc(Rmask); 
colorbar; title('Autocorrelation')

NFFT = 2^10;
SRmask = fftshift(abs(fft2(Rmask,NFFT,NFFT)));

figure; imagesc(10*log10(SRmask)); 
colorbar; title('FFT(Autocorrelation)')

% spatial frequency 
dx = 0.075e-3 /2;
ft_space = (1540/2)/(dx);
ft_space / 2
freq = ((-NFFT/2):(NFFT/2-1))/NFFT*ft_space*1e-6;
figure; plot(10*log10(SRmask(:,NFFT/2)))
figure; plot(freq, 10*log10(SRmask(:,NFFT/2)))

%%
% clear; 
% close all; clc;

pathData = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\NonUniformBSC2D';

normZero = @(x) x-max(x(:));
rf2Bmode = @(x) 20*log10(abs(hilbert(x)));

%filenames = {'reference','sample'}
%filenames = {'rf1','rf2'}
filenames = {'rf7','rf9'}

% filename_ref = 'rf7';  % case1_rfref.mlx | map_sd0p005_bsc1  
% filename_sam = 'rf9';  % case3_rfsam.mlx | mapmask7_sd0p005_bsc2 | mask7: gaussian, size=9, sgm=1

%filenames = {'rf_fs1_mask0','rf_fs1_mask7'}
% Only when the attenuation of sample and reference is the same
ratio  = 2;
for kk = 1:length(filenames)

media = filenames{kk}

f0 = 6.66;   % MHz
size_wl = 10;

width = 10; % Width to smooth laterally for each BoA lines
overlap_pct = 0.5;    % 80%
overlap = round((1-overlap_pct)*width);

col_last = width;
n = 0;

load(fullfile(pathData,media)); 
sam = rf;
% clear sam
%fs = 1/ (1.461038961038961e-08);
%fs = 1/ (7.305194805194804e-09);

NFFT = 2^12;

fNyquist = fs/2;
%h1 = fir1(200,[2e6/(fNyquist) 9e6/(fNyquist)]');

%ratio = 2

env_sam = 20*log10(abs(hilbert(sam(200:end,:))));
bmode_sam = env_sam - max(max(env_sam));
figure; imagesc(bmode_sam); colormap gray;
title("B-mode")


figure; imagesc(x*1e3, z*1e3, normZero(rf2Bmode(rf)),[-60 0])
title('B-mode'); colormap gray
xlabel('Lateral distance (mm)')
ylabel('Depth (mm)')
axis image


%roi = rf(2500:3000,:);
roi = rf(6000/ratio:7000/ratio,:);
% for ff = 1:128
%      roi(:,ff) = conv(roi(:,ff),h1,'same');
% end

env_sam = 20*log10(abs(hilbert(roi)));
bmode_sam = env_sam - max(max(env_sam));
figure; imagesc(bmode_sam); colormap gray;
title('ROI')

% PS = nan(NFFT,128);
% for ii = 54:74
for ii = 1:128
    PS(:,ii) = power(abs(fftshift(fft(roi(:,ii),NFFT))), 2);
end

freq = (-NFFT/2:NFFT/2-1)*(fs*1e-6)/NFFT; % [MHz]
PS_0 = mean(PS(:,1:128),2);
PS = PS_0./max(PS_0);
%PS = 20*log10(PS);
figure; plot(freq,10*log10(PS)); xlim([0 25]); 
ylim([-60 0]), grid on;
xlabel('MHz'); ylabel('Power [dB]'), title('10log(PS)')


% a = 9*(0.0375e-3)/2;  
% c0 = 1540;
% wl = c0./(freq*1e6);
% k = 2*pi./wl; 
% ka =k*a;
%figure; %plot(freq, ka); 

%plot(freq,  10*log10(( 2*besselj(1,2*ka)./(2*ka)).^2) ); xlim([0 25]); 
%hold on; plot(freq, (abs(sin(ka)./ka)).^2);
%hold off;

if kk == 1
PS1 = PS;
PS1_0 = PS_0;
%figure; plot(freq,20*log10(fftshift(abs(fft(h1,NFFT)))))
%xlim([0 20])
end
figure; plot(freq,10*log10(PS./PS1)); 
xlim([0 12])
ylim([-15 0])
xlabel('MHz'); ylabel('Power [dB]'),title('(PS norm)'); %hold on;

figure; plot(freq,10*log10(PS_0./PS1_0)); 
xlim([0 25])
% ylim([-15 0])
xlabel('MHz'); ylabel('Power [dB]'), title('(PS_0)')

dx = 0.3e-3 /16; % technically dy dx of kgrid element_pitch / element_width
ft_space = (1540/2)/(dx);
ft_space / 2
freq_s = ((-NFFT/2):(NFFT/2-1))/NFFT*ft_space*1e-6;
ff = exp(-2 * (sqrt(2) * dx*1e3)^2 * (2*pi/1.54)^2 * freq_s.^2 ); %ok
%ff = exp(-2 * ( dx*1e3)^2 * (2*pi/1.54)^2 * freq_s.^2 ); %ok

figure; plot(freq_s,10*log10(ff),'r'); xlim([1 25])
%ylim([10^-2 10^1])
xlabel('MHz'); title('Form Factor')

end


log_bsc = log(PS_0./PS1_0);
fL = 3.5;
fH = 7.0; % -6 dB
f0 = mean([fL fH]);

[~,idx_fL] = min(abs(freq-fL));
[~,idx_fH] = min(abs(freq-fH));

new_freq = freq(idx_fL:idx_fH);
new_log_bsc = log_bsc(idx_fL:idx_fH);
slope = polyfit(new_freq,new_log_bsc,1)
n = f0 * slope(1); % Equal to n/f0
n


log_bsc = log(ff);
%fL = 3.6;
%fH = 8.7;
f0 = mean([fL fH]);

[~,idx_fL] = min(abs(freq_s-fL));
[~,idx_fH] = min(abs(freq_s-fH));

new_freq = freq_s(idx_fL:idx_fH);
new_log_bsc = log_bsc(idx_fL:idx_fH);
slope = polyfit(new_freq,new_log_bsc,1)
n = f0 * slope(1); % Equal to n/f0
n

%% VERSION EMZ OK

normZero = @(x) x-max(x(:));
rf2Bmode = @(x) 20*log10(abs(hilbert(x)));

% Define filenames
% Reference filename
filename_ref = 'rf7';  % case1_rfref.mlx | map_sd0p005_bsc1  
filename_sam = 'rf9';  % case3_rfsam.mlx | mapmask7_sd0p005_bsc2 | mask7: gaussian, size=9, sgm=1

% Only when the attenuation of sample and reference is the same
ratio  = 2; 

% **Second Power Spectrum Computation (REFERENCE)**
media = filename_ref;
disp(['Processing Reference: ', media])

f0 = 6.66;   % MHz
size_wl = 10;

width = 10; % Width to smooth laterally for each BoA lines
overlap_pct = 0.5;    % 80%
overlap = round((1-overlap_pct)*width);

% RE
load(media); 
sam = rf;

NFFT = 2^12;
fNyquist = fs/2;

% Compute Envelope & B-mode
env_sam = 20*log10(abs(hilbert(sam(200:end,:))));
bmode_sam = env_sam - max(max(env_sam));
figure; imagesc(bmode_sam); colormap gray;
title(['B-mode ', media]);

figure; imagesc(x*1e3, z*1e3, normZero(rf2Bmode(rf)), [-60 0])
title(['B-mode ', media]); colormap gray
xlabel('Lateral distance (mm)')
ylabel('Depth (mm)')
axis image

% Select ROI
roi = rf(6000/ratio:7000/ratio,:);

% Compute Power Spectrum
% PS1 = nan(NFFT, 128);
% for ii = 1:128
%     PS1(:,ii) = power(abs(fftshift(fft(roi(:,ii), NFFT))), 2);
% end
PS1 = abs(fftshift(fft(roi, NFFT, 1))).^2;

freq = (-NFFT/2:NFFT/2-1)*(fs*1e-6)/NFFT; % [MHz]
PS1_0 = mean(PS1, 2);
PS1 = PS1_0 ./ max(PS1_0);

% Plot Power Spectrum
figure; plot(freq, 10*log10(PS1)); xlim([0 25]); 
ylim([-60 0]), grid on;
xlabel('MHz'); ylabel('Power [dB]');
title(['Power Spectrum e', media]);


% **Second Power Spectrum Computation (SAMPLE)**
media = filename_sam;
disp(['Processing Sample: ', media])

load(media); 
sam = rf;

% Compute Envelope & B-mode
env_sam = 20*log10(abs(hilbert(sam(200:end,:))));
bmode_sam = env_sam - max(max(env_sam));
figure; imagesc(bmode_sam); colormap gray;
title(['B-mode ', media]);

figure; imagesc(x*1e3, z*1e3, normZero(rf2Bmode(rf)), [-60 0])
title(['B-mode ', media]); colormap gray
xlabel('Lateral distance (mm)')
ylabel('Depth (mm)')
axis image

% Select ROI
roi = rf(6000/ratio:7000/ratio,:);

% Compute Power Spectrum
% PS2 = nan(NFFT, 128);
% for ii = 1:128
%     PS2(:,ii) = power(abs(fftshift(fft(roi(:,ii), NFFT))), 2);
% end
PS2 = abs(fftshift(fft(roi, NFFT, 1))).^2;

PS2_0 = mean(PS2(:,:),2);
PS2 = PS2_0 ./ max(PS2_0);

% Plot Power Spectrum
figure; plot(freq, 10*log10(PS2)); xlim([0 25]); 
ylim([-60 0]), grid on;
xlabel('MHz'); ylabel('Power e[dB]');
title(['Power Spectrum ', media]);

% Ratio Between Power Spectra
figure; plot(freq, 10*log10(PS2./PS1)); 
xlim([0 12])
ylim([-15 0])
xlabel('MHz'); ylabel('Power [dB]');
title('PS2 / PS1 e');

figure; 
hold on;
plot(freq, 10*log10(PS1_0), 'DisplayName', 'Reference'); 
plot(freq, 10*log10(PS2_0), 'DisplayName', 'Sample'); 
plot(freq, 10*log10(PS2_0./PS1_0), 'DisplayName', 'Ratio'); 
xlim([0 12])
ylim([-15 0])
xlabel('MHz'); ylabel('Power [dB]');
legend('Location','Best')
title('PS2_0 / PS1_0 e');

%%  MASK
% **Form Factor EACH MASK
dx = 0.3e-3 /16; % technically dy dx of kgrid element_pitch / element_width

% MASK7 
sigma = 1; kernel= 9;
mask = fspecial('gaussian', 9, sigma);

title_str = sprintf('Mask Gaussian %dx%d, \\sigma=%d', kernel, kernel, sigma);

figure; 
imagesc(mask);
colorbar; title(['Mask: ', title_str]);

% Compute autocorrelation
Rmask = conv2(mask, mask, 'same');

figure; imagesc(Rmask);
colorbar; title(['Autocorrelation: ', title_str]);

% Compute 2D FFT
SRmask = fftshift(abs(fft2(Rmask, NFFT, NFFT)));

figure; imagesc(20*log10(SRmask));
colorbar; title(['FFT(Autocorrelation): ', title_str]);

% Compute spatial frequency
ft_space = (1540/2) / dx;
freq_s = ((-NFFT/2):(NFFT/2-1)) / NFFT * ft_space * 1e-6;
ff_k = SRmask(:, NFFT/2); clear SR_mask

figure; plot(freq_s, 10*log10(SRmask(:, NFFT/2)), 'k');
xlabel('Frequency (MHz)');
ylabel('Power (dB)');
title(['Frequency Response: ', title_str]);

%% **Form Factor Computation Theoretical**
freq_s = ((-NFFT/2):(NFFT/2-1)) / NFFT * ft_space * 1e-6;
ff_theo = exp(-2 * (sqrt(2) * dx*1e3)^2 * (2*pi/1.54)^2 * freq.^2 ); %ok

figure; plot(freq_s,10*log10(ff_theo),'r'); xlim([1 25])
xlabel('MHz'); title('Form Factor')

% Compare
figure, 
hold on
plot(freq_s,10*log10(ff_k),'b', 'DisplayName','Kernel'); xlim([1 25])
plot(freq_s,10*log10(ff_theo),'r--', 'DisplayName','Theo'); xlim([1 25])
aa = ff_k' ./ ff_theo;
plot(freq_s,aa,'g--', 'DisplayName','Diff'); xlim([1 25])
legend('Location','Best')
xlabel('MHz'); 
title('Form Factor')

%% FINAL PLOT FACTOR
figure; 
hold on;
plot(freq, 10*log10(PS1_0), 'DisplayName', 'Reference'); 
plot(freq, 10*log10(PS2_0), 'DisplayName', 'Sample'); 
plot(freq, 10*log10(PS2_0./PS1_0), 'DisplayName', 'Ratio'); 

plot(freq_s,10*log10(ff_k),'b', 'DisplayName','Kernel'); xlim([1 25])
plot(freq_s,10*log10(ff_theo),'r--', 'DisplayName','Theo'); xlim([1 25])

xlim([0 12])
ylim([-15 0])
xlabel('MHz'); ylabel('Power [dB]');
legend('Location','Best')
title('PS2_0 / PS1_0 e');

%% TEST MASK
% Define the triangular BSC in frequency domain
f = linspace(-10, 10, 1024); % Frequency vector
fc = 5; % Cutoff frequency
BSC = max(0, 1 - abs(f / fc)); % Triangular frequency response

% Inverse Fourier Transform to compute the kernel
kernel = ifftshift(ifft(sqrt(BSC))); % Kernel in spatial domain

% Plot results
figure;
subplot(2,1,1); plot(f, BSC); title('BSC (Triangular)'); xlabel('Frequency (MHz)');
subplot(2,1,2); plot(real(kernel)); title('Kernel in Spatial Domain'); xlabel('Position');


%%
% Define Frequency Domain BSC
NFFT = 64; % Kernel size in Fourier domain
[f_x, f_y] = meshgrid(-NFFT/2:NFFT/2-1, -NFFT/2:NFFT/2-1);
f_c = 0.1; % Cutoff frequency
% BSC = max(0, 1 - abs(f_x / f_c)) .* max(0, 1 - abs(f_y / f_c)); % Triangular BSC

BSC = exp(-2 * (sqrt(2) * dx*1e3)^2 * (2*pi/2.54)^2 * f_x.^2 )

% Compute Kernel (Inverse FFT)
K_custom = fftshift(ifft2(ifftshift(sqrt(BSC))));
K_custom = abs(K_custom); % Take magnitude
K_custom = K_custom / sum(K_custom(:)); % Normalize

% Plot kernel
figure; imagesc(K_custom); colormap hot; colorbar;
title('Custom Kernel from Desired BSC');

%%
clear; close all; clc;

dx = 0.3e-3 /16; % technically dy dx of kgrid element_pitch / element_width

% Parameters

lambda = 2.54e-3; % Characteristic wavelength in meters
NFFT = 64; % Number of frequency points

% Frequency axis
ft_space = (1540/2) / dx;
f_x = ((-NFFT/2):(NFFT/2-1)) / NFFT * ft_space * 1e-6; % Frequency in cycles/m
% f_x = linspace(-1/(2*dx), 1/(2*dx), NFFT); % Frequency in cycles/m

% Define the theoretical BSC in 1D
BSC_1D = exp(-2 * (sqrt(2) * dx * 1e3)^2 * (2*pi / (lambda * 1e3))^2 * f_x.^2);

% Compute 1D kernel via inverse FFT
K_1D = fftshift(ifft(ifftshift(sqrt(BSC_1D))));
K_1D = abs(K_1D); % Magnitude of the kernel
K_1D = K_1D / max(K_1D); % Normalize

% Estimate sigma1D
x = (-NFFT/2:NFFT/2-1) * dx; % Spatial domain
sigma_1D = sqrt(sum(x.^2 .* K_1D) / sum(K_1D)); % Weighted standard deviation

% Plot results
figure;
subplot(2, 1, 1);
plot(f_x * 1e-6, BSC_1D, 'r'); grid on;
title('Theoretical BSC (1D)'); xlabel('f_x (MHz)'); ylabel('Magnitude');

subplot(2, 1, 2);
plot(x * 1e6, K_1D, 'b'); grid on;
title('1D Kernel in Spatial Domain'); xlabel('x (µm)'); ylabel('Amplitude');
disp(['Estimated sigma1D (meters): ', num2str(sigma_1D)]);

%%
% Parameters
kernel_size = 9; % Size of the 2D kernel (e.g., 11x11 pixels)
sigma_2D = sigma_1D; % Use the same sigma from 1D kernel

% Create 2D Gaussian kernel
[x, y] = meshgrid(-floor(kernel_size/2):floor(kernel_size/2));
K_2D = exp(-(x.^2 + y.^2) / (2 * (sigma_2D / dx)^2)); % Gaussian function
K_2D = K_2D / sum(K_2D(:)); % Normalize

% Plot 2D Gaussian kernel
figure;
subplot(1, 2, 1);
imagesc(K_2D); colormap hot; colorbar;
title('2D Gaussian Kernel');

subplot(1, 2, 2);
surf(K_2D); shading interp;
title('2D Gaussian Kernel (3D View)');
xlabel('x (pixels)'); ylabel('y (pixels)'); zlabel('Amplitude');

%%
% clear; close all; clc;

% Parameters
dx = 1.8750e-5; % Spatial resolution in meters
lambda = 2.54e-3; % Characteristic wavelength in meters
NFFT = 2^10; % Number of frequency points

% Frequency axis
ft_space = (1540/2) / dx;
% f_x = ((-NFFT/2):(NFFT/2-1)) / NFFT * ft_space * 1e-6; % Frequency in cycles/m
f_x = linspace(-1/(2*dx), 1/(2*dx), NFFT); % Frequency axis in cycles/m


% Theoretical BSC
% BSC_1D = exp(-2 * (sqrt(2) * dx * 1e3)^2 * (2*pi / (lambda * 1e3))^2 * f_x.^2);
BSC_1D = ff_k';

% Compute 1D kernel
K_1D = fftshift(ifft(ifftshift(sqrt(BSC_1D))));
K_1D = abs(K_1D); % Take magnitude
K_1D = K_1D / max(K_1D); % Normalize

% Spatial axis
x = (-NFFT/2:NFFT/2-1) * dx; % Spatial domain
% x = ((-NFFT/2):(NFFT/2-1)) / NFFT * ft_space * 1e-6; % Frequency in cycles/m

% Compute sigma1D
sigma_1D = sqrt(sum(x.^2 .* K_1D) / sum(K_1D));
disp(['Estimated sigma1D: ', num2str(sigma_1D)]);

% Plot results
figure;
subplot(2, 1, 1);
plot(f_x * 1e-6, BSC_1D, 'r'); grid on;
title('Theoretical BSC (1D)'); xlabel('f_x (MHz)'); ylabel('Magnitude');

subplot(2, 1, 2);
plot(x * 1e6, K_1D, 'b'); grid on;
title('1D Kernel in Spatial Domain'); xlabel('x (\mum)'); ylabel('Amplitude');

%% NEWW INVRSE WAY


% clear; close all; clc;

% Given parameters
kernel_size = 9; % Desired kernel size
NFFT = 64; % FFT size (must match the original one)
dx = 0.3e-3 / 16; % Spatial resolution (m)

% Create a 2D frequency grid
[f_x, f_y] = meshgrid(linspace(-1/(2*dx), 1/(2*dx), NFFT));
f_r = sqrt(f_x.^2 + f_y.^2); % Radial frequency grid

% Assume `ff_k` is a given 1D frequency response
% Replace this with any empirical `ff_k`
freq_1D = linspace(0, max(f_r(:)), length(ff_k)); % 1D frequency axis
ff_k_r = interp1(freq_1D, ff_k, f_r, 'linear', 0); % Expand into 2D

% Ensure symmetry in the frequency domain
SRmask = ff_k_r; 

% Compute 2D kernel via inverse FFT
kernel_full = abs(ifft2(ifftshift(sqrt(SRmask)))); % Inverse FFT
kernel_full = fftshift(kernel_full); % Center the kernel properly
% kernel_full = kernel_full / max(kernel_full(:)); % Normalize

% Find true center
center_idx = ceil(NFFT / 2); % Correct center index

% Crop the kernel correctly
half_size = floor(kernel_size / 2);
kernel = kernel_full(center_idx - half_size +1: center_idx + half_size +1, ...
                     center_idx - half_size +1: center_idx + half_size +1);

% Plot results
figure;
subplot(1, 2, 1);
imagesc(SRmask); colorbar; 
title('Given Frequency Response (SRmask)');

subplot(1, 2, 2);
imagesc(kernel); colormap hot; colorbar;
title(['Recovered Kernel (Size: ', num2str(kernel_size), 'x', num2str(kernel_size), ')']);



% first obtain a kernel1D and then perhaps you can create a 2D one