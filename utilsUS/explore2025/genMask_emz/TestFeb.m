% CALCULATE POWER SPRECTRUM SAM Y REFERENCE
% ALSO NORMALIZED POWER SPECTRUM

%% **Form Factor Computation Theoretical**
freq_s = ((-NFFT/2):(NFFT/2-1)) / NFFT * ft_space * 1e-6;
ff_theo = exp(-2 * (sqrt(2) * dx*1e3)^2 * (2*pi/1.54)^2 * freq.^2 ); %ok

figure; plot(freq_s,10*log10(ff_theo),'r'); xlim([1 25])
xlabel('MHz'); title('Form Factor')

% Compare
figure, 
hold on
% plot(freq_s,10*log10(ff_k),'b', 'DisplayName','Kernel'); xlim([1 25])
plot(freq_s,10*log10(ff_theo),'r--', 'DisplayName','Theo'); xlim([1 25])
% aa = ff_k' ./ ff_theo;
% plot(freq_s,aa,'g--', 'DisplayName','Diff'); xlim([1 25])
legend('Location','Best')
xlabel('MHz'); 
title('Form Factor')

% VERSION EMZ
%% GAUSSIAN SIMULATIONS ANDRES COILA
pathData = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\NonUniformBSC2D';


normZero = @(x) x-max(x(:));
rf2Bmode = @(x) 20*log10(abs(hilbert(x)));

% Define filenames
% Reference filename
filename_ref = 'rf7';  % case1_rfref.mlx | map_sd0p005_bsc1  
filename_sam = 'rf9';  % case3_rfsam.mlx | mapmask7_sd0p005_bsc2 | mask7: gaussian, size=9, sgm=1

% Power Spectrum Computation (REFERENCE)**
media = filename_ref;
disp(['Processing Reference: ', media])

f0 = 6.66;   % MHz
size_wl = 10;

width = 10; % Width to smooth laterally for each BoA lines
overlap_pct = 0.5;    % 80%
overlap = round((1-overlap_pct)*width);

% REF
load(fullfile(pathData, media)); 
sam = rf;

NFFT = 2^12;
fNyquist = fs/2;

% Compute Envelope & B-mode
env_sam = 20*log10(abs(hilbert(sam(200:end,:))));
bmode_sam = env_sam - max(max(env_sam));
% figure; imagesc(bmode_sam); colormap gray;
% title(['B-mode ', media]);

figure; imagesc(x*1e3, z*1e3, normZero(rf2Bmode(rf)), [-60 0])
title(['B-mode REF ', media]); colormap gray
xlabel('Lateral distance (mm)')
ylabel('Depth (mm)')
axis image

ratio = 1;
% Select ROI
roi = rf(6000/ratio:8000/ratio,:);

% Compute Power Spectrum
% PS1 = nan(NFFT, 128);
% for ii = 1:128
%     PS1(:,ii) = power(abs(fftshift(fft(roi(:,ii), NFFT))), 2);
% end
PSref = abs(fftshift(fft(roi, NFFT, 1))).^2;

freq = (-NFFT/2:NFFT/2-1)*(fs*1e-6)/NFFT; % [MHz]
PSref = mean(PSref, 2);
PSref_0 = PSref ./ max(PSref);

% Plot Power Spectrum
figure; plot(freq, 10*log10(PSref_0)); xlim([0 25]); 
ylim([-60 0]), grid on;
xlabel('MHz'); ylabel('Power [dB]');
title(['Power Spectrum Ref ', media]);


% Power Spectrum Computation (SAMPLE)**
media = filename_sam;
disp(['Processing Sample: ', media])

load(fullfile(pathData,media)); 
sam = rf;

% Compute Envelope & B-mode
env_sam = 20*log10(abs(hilbert(sam(200:end,:))));
bmode_sam = env_sam - max(max(env_sam));
% figure; imagesc(bmode_sam); colormap gray;
% title(['B-mode ', media]);

figure; imagesc(x*1e3, z*1e3, normZero(rf2Bmode(rf)), [-60 0])
title(['B-mode SAM ', media]); colormap gray
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
PSsam = abs(fftshift(fft(roi, NFFT, 1))).^2;

PSsam = mean(PSsam, 2);
PSsam_0 = PSsam ./ max(PSsam);

% Plot Power Spectrum
figure; plot(freq, 10*log10(PSsam_0)); xlim([0 25]); 
ylim([-60 0]), grid on;
xlabel('MHz'); ylabel('Power e[dB]');
title(['Power Spectrum Sam ', media]);

% Plot PWELCH
blocksize_wv = 20;
lambda = 1540/6.6*1E-6;
overlap = 0.9;
ratio_zx = 1;
dz = z(2) - z(1);
nz = round(blocksize_wv*lambda/dz /2 * ratio_zx); % Window size
wz = round(blocksize_wv*lambda*(1-overlap)/dz * ratio_zx); % Between windows
    
[pxx,fpxx] = pwelch(roi-mean(roi),nz,nz-wz,nz,fs);
meanSpectrum = mean(pxx,2);
meanSpectrum(1) = 0;

figure,
plot(fpxx/1e6, db(meanSpectrum/max(meanSpectrum))),grid on
xlim([0, fs/2e6])
hold on
xline(freq_L/1e6, 'k--')
xline(freq_H/1e6, 'k--')
hold off
xlabel('Frequency [MHz]')
ylabel('Magnitude [dB]')


%% RATIO POWER SPECTRUM

PSratio = PSsam ./ PSref;
PSratio_0 = PSratio ./ max(PSratio);

% Ratio Between Power Spectra
figure; 
plot(freq, 10*log10(PSratio), 'b--'); 
xlim([0 12])
% ylim([-15 0])
xlabel('MHz'); ylabel('Power [dB]');
title('PS_{sam} / PS_{ref}');

figure; 
plot(freq, 10*log10(PSratio_0), 'b--'); 
xlim([0 12])
ylim([-15 0])
xlabel('MHz'); ylabel('Power [dB]');
title('Norm PS_{sam} / PS_{ref}');


% figure; 
% hold on;
% plot(freq, 10*log10(PSref_0), '--', 'DisplayName', 'Reference'); 
% plot(freq, 10*log10(PSsam_0), '--', 'DisplayName', 'Sample'); 
% plot(freq, 10*log10(PSsam_0./PSref_0), '--', 'DisplayName', 'Ratio'); 
% xlim([0 12])
% ylim([-15 0])
% xlabel('MHz'); ylabel('Power [dB]');
% legend('Location','Best')
% title('PSsam_0 / PSref_0');

%%

%% GAUSSIAN SIMULATIONS EMZ

list_sigma = [0.5 1 1.25 1.5 2];
colors = lines(length(list_sigma));  % Use MATLAB's built-in distinct colors

for ss = 1:length(list_sigma)
sigmaValue = list_sigma(ss);
pathData = ['C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\' ...
    'dataACS_kwave\simFeb2025'];

normZero = @(x) x-max(x(:));
rf2Bmode = @(x) 20*log10(abs(hilbert(x)));

% Define filenames
alpha_sam = 0.5; 
j_sam = 1.1;

alpha_ref = 0.5;
j_ref = j_sam;

%%%%%%%%%%%%%%%%%%%%%%%%%% REFERENCE %%%%%%%%%%%%%%%%%%%%%%%%%%%

folderDataRef = ['Gauss',num2str(sigmaValue), ...
                '_a',num2str(alpha_ref),'_pow', num2str(j_ref)];
folderDataRef = strrep(folderDataRef, '.', 'p');

rf_ref_name = strcat('rfref_', sprintf('%.3f', j_ref));
rf_ref_name = strrep(rf_ref_name, '.', 'p');
rf_ref_name = strcat(rf_ref_name, '.mat');

load(fullfile(pathData, folderDataRef, rf_ref_name)); 

% Compute Envelope & B-mode
% figure;
% imagesc(x*1e3, z*1e3, normZero(rf2Bmode(rf)), [-60 0])
% % imagesc(normZero(rf2Bmode(rf)), [-60 0])
% title(['B-mode REF ', rf_ref_name]); colormap gray
% xlabel('Lateral distance (mm)')
% ylabel('Depth (mm)')
% % axis image

% REF
NFFT = 2^12;
fNyquist = fs/2;

% Select ROI
% roi = rf;
% roi = rf(6000/ratio:7000/ratio,:);
% roi = rf(2000:6000,:);

x_ini = -15e-3; x_end = 15e-3;
z_ini = 10e-3;  z_end = 40e-3;
ind_x = x_ini <= x & x <= x_end;
ind_z = z_ini <= z & z <= z_end;
roi = rf(ind_z,ind_x);


% Compute Power Spectrum
% PS1 = nan(NFFT, 128);
% for ii = 1:128
%     PS1(:,ii) = power(abs(fftshift(fft(roi(:,ii), NFFT))), 2);
% end
PSref = abs(fftshift(fft(roi, NFFT, 1))).^2;

freq = (-NFFT/2:NFFT/2-1)*(fs*1e-6)/NFFT; % [MHz]
PSref = mean(PSref, 2);
PSref_0 = PSref ./ max(PSref);

% Plot Power Spectrum
% figure; plot(freq, 10*log10(PSref)); xlim([0 25]); 
% ylim([-60 0]), grid on;
% xlabel('MHz'); ylabel('Power [dB]');
% title(['Power Spectrum Ref ', rf_ref_name]);

%%%%%%%%%%%%%%%%%%%%%%%%%%% REFERENCE %%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% SAMPLE %%%%%%%%%%%%%%%%%%%%%%%%%%%
folderDataSam = ['Gauss',num2str(sigmaValue), ...
                '_a',num2str(alpha_sam),'_pow', num2str(j_sam)];
folderDataSam = strrep(folderDataSam, '.', 'p');

rf_sam_name = strcat('rfsam_', sprintf('%.3f', j_sam));
rf_sam_name = strrep(rf_sam_name, '.', 'p');
rf_sam_name = strcat(rf_sam_name, '.mat');

load(fullfile(pathData, folderDataSam, rf_sam_name));

% figure; imagesc(x*1e3, z*1e3, normZero(rf2Bmode(rf)), [-60 0])
% title(['B-mode SAM ', rf_sam_name]); colormap gray
% xlabel('Lateral distance (mm)')
% ylabel('Depth (mm)')
% axis image

% Select ROI
% roi = rf;
% roi = rf(6000/ratio:7000/ratio,:);
% roi = rf(2000:6000,:);

x_ini = -15e-3; x_end = 15e-3;
z_ini = 10e-3;  z_end = 40e-3;
ind_x = x_ini <= x & x <= x_end;
ind_z = z_ini <= z & z <= z_end;
roi = rf(ind_z,ind_x);

% Compute Power Spectrum
% PS2 = nan(NFFT, 128);
% for ii = 1:128
%     PS2(:,ii) = power(abs(fftshift(fft(roi(:,ii), NFFT))), 2);
% end
PSsam = abs(fftshift(fft(roi, NFFT, 1))).^2;

PSsam = mean(PSsam, 2);
PSsam_0 = PSsam ./ max(PSsam);

% Plot Power Spectrum
% figure; plot(freq, 10*log10(PSsam)); xlim([0 25]); 
% ylim([-60 0]), grid on;
% xlabel('MHz'); ylabel('Power e[dB]');
% title(['Power Spectrum Sam ', rf_sam_name]);

%%%%%%%%%%%%%%%%%%%%%%%%%%% SAMPLE %%%%%%%%%%%%%%%%%%%%%%%%%%%


% RATIOS
PSratio = PSsam ./ PSref;
PSratio_0 = PSratio ./ max(PSratio);

sigmaValue2 = 2*sigmaValue;
% MASK 
NFFT = 512;
kernel=9;
mask = fspecial('gaussian', kernel, sigmaValue);

SRmask_v1 = abs(fftshift(fft2(mask, NFFT, NFFT))).^2;
ps_gauss2d_v1 = SRmask_v1(:, NFFT/2);

Rmask = conv2(mask, mask, 'same');
SRmask_v2 = abs(fftshift(fft2(Rmask, NFFT, NFFT)));
ps_gauss2d_v2 = SRmask_v2(:, NFFT/2);

% THEORETICAL
c0 = 1540;
dx = 0.3e-3 /16;
ft_space = c0/2/dx;
NFFT = 512;

freq_s = ((-NFFT/2):(NFFT/2-1))/NFFT*ft_space; %[Hz]
freq_MHz = freq_s*1e-6;
f_theo = exp(-(sigmaValue2)^2 *(2*pi/c0 *freq_s *dx).^2 ); %ok



currentColor = colors(ss, :); 

figure(177)

% --- Plot Power Spectrum Ratio ---
plot(freq, 10*log10(PSratio_0), '.-', 'Color', currentColor, 'LineWidth', 1.5, ...
        'DisplayName', ['PS Ratio, \sigma=', num2str(sigmaValue)], 'MarkerSize', 20);
hold on;

% --- Plot Mask Spectrum ---
plot(freq_MHz, 10*log10(ps_gauss2d_v1), '--', 'Color', currentColor, 'LineWidth', 1.5, ...
        'DisplayName', ['Maskv1, \sigma=', num2str(sigmaValue)]);
hold on;
% --- Plot Theoretical Spectrum ---
plot(freq_MHz, 10*log10(f_theo), '+-', 'Color', currentColor, 'LineWidth', 1.5, ...
        'DisplayName', ['Theo, \sigma=', num2str(sigmaValue)]);
hold on;
% Update legend and labels
legend('Location', 'Best');
xlabel('Frequency [MHz]');
ylabel('Power [dB]');
title('Power Spectra Comparison for Different \sigma Values');
grid on;

xlim([-20 20]);
xlim([0 20]);
ylim([-100 0]);

end
%%
c0 = 1540;
dx = 0.3e-3 /16;
ft_space = c0/2/dx;
NFFT = 512;

freq_s = ((-NFFT/2):(NFFT/2-1))/NFFT*ft_space; %[Hz]
freq_MHz = freq_s*1e-6;
list_sigma = [0.5 1 1.25 1.5 2];


for ss = 1:length(list_sigma)

    sigma = list_sigma(ss);
    f_theo = exp(-sigma^2 *(2*pi/c0 *freq_s *dx).^2 ); %ok
    figure(11)
    plot(freq_MHz, 10*log10(f_theo), '-.', 'LineWidth', 1.5, 'DisplayName', ['Theo \sigma=', num2str(sigma)]);
    hold on, grid on,
    title('Theoretical Gauss ')
    legend('Location', 'Best')
    xlabel('MHz'); ylabel('Power [dB]');
    grid on;
end
%%

%% RATIO POWER SPECTRUM

PSratio = PSsam ./ PSref;
PSratio_0 = PSratio ./ max(PSratio);

% Ratio Between Power Spectra
figure; 
plot(freq, 10*log10(PSratio), 'b--'); 
xlim([0 12])
% ylim([-15 0])
xlabel('MHz'); ylabel('Power [dB]');
title('PS_{sam} / PS_{ref}');

figure; 
plot(freq, 10*log10(PSratio_0), 'b--'); 
xlim([0 12])
% ylim([-15 0])
xlabel('MHz'); ylabel('Power [dB]');
title('Norm PS_{sam} / PS_{ref}');


figure; 
hold on;
plot(freq, 10*log10(PSref_0), '--', 'DisplayName', 'Reference'); 
plot(freq, 10*log10(PSsam_0), '--', 'DisplayName', 'Sample'); 
plot(freq, 10*log10(PSsam_0./PSref_0), '--', 'DisplayName', 'Ratio'); 
xlim([0 12])
ylim([-15 0])
xlabel('MHz'); ylabel('Power [dB]');
legend('Location','Best')
title('PSsam_0 / PSref_0');

%%

kernel = 9;
mask = fspecial('gaussian', kernel, sigmaValue);

% Use full convolution and ensure proper centering
Rmask = conv2(mask, mask, 'same');  % Full autocorrelation

% Apply 2D FFT on the centered autocorrelation
SRmask_v2 = abs(fftshift(fft2(Rmask, NFFT, NFFT)));
ps_gauss2d_v2 = SRmask_v2(:, NFFT/2);

% Compute spatial frequency
dx = 0.3e-3 /16; % technically dy dx of kgrid element_pitch / element_width

c0 = 1540;
ft_space = c0/2/dx;
freq_s = ((-NFFT/2):(NFFT/2-1))/NFFT * ft_space * 1e-6; % [MHz]

figure, 
plot(freq_s, pow2db(ps_gauss2d_v2), 'r', 'DisplayName','conv')
grid on;



