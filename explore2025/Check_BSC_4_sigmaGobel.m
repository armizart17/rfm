%%% TEST FAST GOBEL SIMULATIONS DIFFERENT SIGMAS 
% OF GOBEL KERNEL

%% SPECTRAL METHOD PARAMETERS
pars.P = 4096; % NFFT only for calculate BSC_RPM_ok 10wl
pars.bw          = [3 9]; % [MHz] % new
pars.overlap     = 0.8;
pars.blocksize   = 12; % wavelengths

% new
pars.z_roi       = [5 45]*1E-3; % all
pars.x_roi       = [-18 18]*1E-3;

pars.window_type = 3; %  (1) Hanning, (2) Tuckey, (3) Hamming, (4) Tchebychev
pars.saran_layer = false;


%% LOAD DATA CHANGE IN SIGMA

list_sigma = [1.25 1.5 1.75];
f0 = 150;

% First row for headers, second for data
bsc_rpm_sigma = cell(2, length(list_sigma)); 

% Store headers "sigma"
bsc_rpm_sigma(1, :) = num2cell(list_sigma);

for  iSigma = 1: length(list_sigma)
    
sigmaValue = list_sigma(iSigma);

% DATA NEW AMZ
pathData = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\dataACS_kwave\simGabor\gaborV1';

alpha_sam = 0.5; 
j_sam = 1.1;

alpha_ref = 0.5;
j_ref = j_sam;

folderDataSam = sprintf('Gabor_s%.3g_f%.3g', sigmaValue, f0);
folderDataSam = strrep(folderDataSam, '.', 'p');

rf_sam_name = strcat('rfsam_', sprintf('%.3f', j_sam));
rf_sam_name = strrep(rf_sam_name, '.', 'p');
rf_sam_name = strcat(rf_sam_name, '.mat');

SAM = load(fullfile(pathData, folderDataSam, rf_sam_name));
SAM.alpha_power = j_sam;
SAM.acs = alpha_sam; % [dB/cm/MHz] 

folderDataRef = sprintf('Gabor_s%.3g_f%.3g', sigmaValue, f0);
folderDataRef = strrep(folderDataRef, '.', 'p');

rf_ref_name = strcat('rfref_', sprintf('%.3f', j_ref));
rf_ref_name = strrep(rf_ref_name, '.', 'p');
rf_ref_name = strcat(rf_ref_name, '.mat');

REF = load(fullfile(pathData, folderDataRef, rf_ref_name)); 
REF.alpha_power = j_ref; 
REF.acs = alpha_ref; % [dB/cm/MHz]

% BSC RPM METHOD
BSC = calculateBSC_RPM_ok(SAM, REF, pars); % slow only once**
% BSC = calculateBSC_RPM_fast(SAM, REF, pars); % fast TBD**
%
bsc_rpm = BSC;
% bsc_rpm = BSC.BSCcurve_Uni(:,2); % median

bsc_rpm_sigma{2, iSigma} = bsc_rpm;

end

%%
keyboard
%% PLOT ALL BSC RPM SIGMAS
% load('C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\dataACS_kwave\simFeb2025\bscOK_gobelSimu_sigmaEffect_f150.mat')
freq = BSC.band(:);

line_width = 2;
font_size = 14;
figure, 
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 700, 700]); % [x, y, width, height] in pixels

for iSigma = 1: length(list_sigma)

    nameLineBSC = ['\sigma = ', num2str(bsc_rpm_sigma{1, iSigma})];
    
    % semilogy(freq, bsc_rpm_sigma{2, iSigma}, '-',  'LineWidth', line_width, 'DisplayName', nameLineBSC);
    semilogy(freq, bsc_rpm_sigma{2, iSigma}.BSCcurve_Uni(:,2), '-',  'LineWidth', line_width, 'DisplayName', nameLineBSC);
    
    hold on, grid on;

end

xlabel('Frequency [MHz]', 'FontSize', font_size);
ylabel('BSC [cm^{-1}\cdot sr^{-1}]', 'FontSize', font_size);
title('Backscatter Coefficient (BSC)', 'FontSize', font_size + 2);
legend('Location', 'best', 'FontSize', font_size);
set(gca, 'FontSize', font_size);

%% PROCESS A PRIOR FITTING
% load('C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\dataACS_kwave\simFeb2025\bscOK_gaussSimu_sigmaEffect.mat')
for iSigma = 1: size(bsc_rpm_sigma, 2)

bsc_rpm = bsc_rpm_sigma{2, iSigma}.BSCcurve_Uni(:,2); % median
freq = bsc_rpm_sigma{2, iSigma}.band(:);
sigmaValue = bsc_rpm_sigma{1, iSigma};

fprintf( '----- (Gauss Kernel σ = %.2f) -----\n', sigmaValue);

% CODE
% GAUSSIAN MODEL (g.exp(-s.f^2))

% Linearize
% Perform linear regression  bsc = - d_s . f^2 + ln(d_g) 
coeffs   = polyfit(-freq.^2, log(bsc_rpm), 1); % Fit y = mx + c
d_s      = coeffs(1); % Slope = d_s -0.0319 (mean), -0.0317 (median)
ln_dg    = coeffs(2); % Intercept = ln(d_g) 
d_g      = exp(ln_dg); % 1.0917 (mean), 0.9079(median)

% Display results
fprintf('-----RPM Gauss (g.exp(-s.f^2))-----\n')
fprintf('Δg           = %f\n', d_g);
fprintf('g_s/g_r [dB] = %f\n', 10*log10(d_g));
fprintf('Δs           = %f\n', d_s);
fprintf('---------\n')

% Reconstruct BSC
bsc_gauss_reconstruct = d_g*exp(-d_s*freq.^2);

% POWER LAW MODEL (b.(f^n))
% Perform linear regression  bsc = d_n . log(f) + ln(d_b) 
coeffs_pl   = polyfit(log(freq), log(bsc_rpm), 1); % Fit y = m.x + c
d_n_pl      = coeffs_pl(1); % Slope = d_n  (mean),  (median)
ln_d_b_pl   = coeffs_pl(2); % Intercept = ln(d_b) 
d_b_pl      = exp(ln_d_b_pl); % 1.0917 (mean), 0.9079(median)

% Display results
fprintf('-----RPM POWER LAW (b.(f^n))-----\n')
fprintf('Δb          = %f\n', d_b_pl);
fprintf('b_s/b_r [dB] = %f\n', 10*log10(d_b_pl));
fprintf('Δn           = %f\n', d_n_pl);
fprintf('---------\n')

bsc_powlaw_reconstruct = d_b_pl * freq.^d_n_pl;

%%%%%%%% "ka" no g %%%%%%%%
a2s     = @(a,c0) 0.827 *(2*pi/c0)^2 * a^2;
s2a     = @(s,c0) sqrt(s/0.827)*c0/(2*pi);
s2ka    = @(s,fc) sqrt(s/0.827)*fc;

fc = 6.6; % [MHz]
c0 = 1540; % [m/s]

a_nog = s2a(d_s, c0);
ka_nog = s2ka(d_s, fc);
bsc_nog = exp(-0.827* (2*pi*a_nog/c0)^2 * freq.^2 );

fprintf('-----ka (no g) -----\n')
fprintf('a  [µm]      = %.2f\n', a_nog );
fprintf('ka @ %.2fMHz = %.2f\n', fc, ka_nog );

%%%%%%%% "ka" with g %%%%%%%%
yy = log(bsc_gauss_reconstruct);
% yy = log(d_g) - d_s * freq.^2;

AA = -0.827 *(2*pi/c0)^2 * freq.^2;

x_opt = cgs(AA'*AA, AA'*yy); % AA \ yy

a_g = sqrt(x_opt);
ka_g = (2*pi*fc/c0)*a_g; 
bsc_g = exp(-0.827* (2*pi*a_g/c0)^2 * freq.^2 );

fprintf('-----ka (with g) -----\n')
fprintf('a  [µm]      = %.2f\n', a_g );
fprintf('ka @ %.2fMHz = %.2f\n', fc, ka_g );

% a method var minimization search
a_minVar = bsc_rpm_sigma{2, iSigma}.ESD_Uni(1)/2*1e6; % [um]
fprintf('-----a (minVar) -----\n')
fprintf('a  [µm]      = %.2f\n', a_minVar );

% ALL PLOT
line_width = 2;
font_size = 20;
figure, 
% set(gcf, 'Units', 'pixels', 'Position', [100, 100, 700, 700]); % [x, y, width, height] in pixels

semilogy(freq, bsc_rpm, 'k-', 'DisplayName', 'GT', 'LineWidth', line_width);
hold on; grid on
semilogy(freq, bsc_gauss_reconstruct, 'g--', 'DisplayName', 'Gauss BSC', 'LineWidth', line_width) ;
semilogy(freq, bsc_powlaw_reconstruct, 'b--', 'DisplayName', 'PowLaw BSC', 'LineWidth', line_width) ;

xlabel('Frequency [MHz]', 'FontSize', font_size);
ylabel('BSC [cm^{-1}\cdot sr^{-1}]', 'FontSize', font_size);
title(['BSC \sigma=', num2str(sigmaValue)])
legend('Location', 'southwest', 'FontSize', font_size);
hold off;
set(gca,'fontsize',16)

end

%% CHANGE IN F0 MARCH 2025

% LOAD DATA CHANGE IN SIGMA

list_f0MHz = [10 12.5 15];
f0 = 150;

% First row for headers, second for data
bsc_rpm_f0 = cell(2, length(list_f0MHz)); 

% Store headers "sigma"
bsc_rpm_f0(1, :) = num2cell(list_f0MHz);

% DATA NEW AMZ
pathData = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\dataACS_kwave\simGabor\gaborV2';

for  iCount = 1: length(list_f0MHz)
    
f0 = list_f0MHz(iCount);
sigmaValue = 1.25;

alpha_sam = 0.5; 
j_sam = 1.1;

alpha_ref = 0.5;
j_ref = j_sam;

folderDataSam = sprintf('Gabor_s%.3g_f%.3g', sigmaValue, f0);
folderDataSam = strrep(folderDataSam, '.', 'p');

rf_sam_name = strcat('rfsam_', sprintf('%.3f', j_sam));
rf_sam_name = strrep(rf_sam_name, '.', 'p');
rf_sam_name = strcat(rf_sam_name, '.mat');

SAM = load(fullfile(pathData, folderDataSam, rf_sam_name));
SAM.alpha_power = j_sam;
% fprintf('Density (kg/m3) %.g ± %g \n', mean(SAM.medium.density(:)), std(SAM.medium.density(:)) );

SAM.acs = alpha_sam; % [dB/cm/MHz] 

folderDataRef = sprintf('Gabor_s%.3g_f%.3g', sigmaValue, f0);
folderDataRef = strrep(folderDataRef, '.', 'p');

rf_ref_name = strcat('rfref_', sprintf('%.3f', j_ref));
rf_ref_name = strrep(rf_ref_name, '.', 'p');
rf_ref_name = strcat(rf_ref_name, '.mat');

REF = load(fullfile(pathData, folderDataRef, rf_ref_name)); 
REF.alpha_power = j_ref; 
REF.acs = alpha_ref; % [dB/cm/MHz]

%% BSC RPM METHOD
% BSC = calculateBSC_RPM_ok(SAM, REF, pars); % slow only once**
% BSC = calculateBSC_RPM_fast(SAM, REF, pars); % fast TBD**
%
% bsc_rpm = BSC;
% bsc_rpm = BSC.BSCcurve_Uni(:,2); % median

% bsc_rpm_f0{2, iCount} = bsc_rpm;

end
%%
keyboard
%% PLOT ALL BSC RPM SIGMAS
% load('C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\dataACS_kwave\simFeb2025\bscOK_gobelSimu_sigmaEffect_f150.mat')
freq = BSC.band(:);

line_width = 2;
font_size = 14;
figure, 
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 700, 700]); % [x, y, width, height] in pixels

for iCount = 1: length(list_f0MHz)

    nameLineBSC = ['f_o = ', num2str(bsc_rpm_f0{1, iCount}/2)];
    
    % semilogy(freq, bsc_rpm_sigma{2, iCount}, '-',  'LineWidth', line_width, 'DisplayName', nameLineBSC);
    semilogy(freq, bsc_rpm_f0{2, iCount}.BSCcurve_Uni(:,2), '-',  'LineWidth', line_width, 'DisplayName', nameLineBSC);
    
    hold on, grid on;

end

xlabel('Frequency [MHz]', 'FontSize', font_size);
ylabel('BSC [cm^{-1}\cdot sr^{-1}]', 'FontSize', font_size);
title('Backscatter Coefficient (BSC)', 'FontSize', font_size + 2);
legend('Location', 'best', 'FontSize', font_size);
set(gca, 'FontSize', font_size);

%% PROCESS A PRIOR FITTING
% load('C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\dataACS_kwave\simFeb2025\bscOK_gaussSimu_sigmaEffect.mat')
for iCount = 1: size(bsc_rpm_f0, 2)

bsc_rpm = bsc_rpm_f0{2, iCount}.BSCcurve_Uni(:,2); % median
freq = bsc_rpm_f0{2, iCount}.band(:);
sigmaValue = bsc_rpm_f0{1, iCount};

fprintf( '----- (Gauss Kernel σ = %.2f) -----\n', sigmaValue);

% CODE
% GAUSSIAN MODEL (g.exp(-s.f^2))

% Linearize
% Perform linear regression  bsc = - d_s . f^2 + ln(d_g) 
coeffs   = polyfit(-freq.^2, log(bsc_rpm), 1); % Fit y = mx + c
d_s      = coeffs(1); % Slope = d_s -0.0319 (mean), -0.0317 (median)
ln_dg    = coeffs(2); % Intercept = ln(d_g) 
d_g      = exp(ln_dg); % 1.0917 (mean), 0.9079(median)

% Display results
fprintf('-----RPM Gauss (g.exp(-s.f^2))-----\n')
fprintf('Δg           = %f\n', d_g);
fprintf('g_s/g_r [dB] = %f\n', 10*log10(d_g));
fprintf('Δs           = %f\n', d_s);
fprintf('---------\n')

% Reconstruct BSC
bsc_gauss_reconstruct = d_g*exp(-d_s*freq.^2);

% POWER LAW MODEL (b.(f^n))
% Perform linear regression  bsc = d_n . log(f) + ln(d_b) 
coeffs_pl   = polyfit(log(freq), log(bsc_rpm), 1); % Fit y = m.x + c
d_n_pl      = coeffs_pl(1); % Slope = d_n  (mean),  (median)
ln_d_b_pl   = coeffs_pl(2); % Intercept = ln(d_b) 
d_b_pl      = exp(ln_d_b_pl); % 1.0917 (mean), 0.9079(median)

% Display results
fprintf('-----RPM POWER LAW (b.(f^n))-----\n')
fprintf('Δb          = %f\n', d_b_pl);
fprintf('b_s/b_r [dB] = %f\n', 10*log10(d_b_pl));
fprintf('Δn           = %f\n', d_n_pl);
fprintf('---------\n')

bsc_powlaw_reconstruct = d_b_pl * freq.^d_n_pl;

%%%%%%%% "ka" no g %%%%%%%%
a2s     = @(a,c0) 0.827 *(2*pi/c0)^2 * a^2;
s2a     = @(s,c0) sqrt(s/0.827)*c0/(2*pi);
s2ka    = @(s,fc) sqrt(s/0.827)*fc;

fc = 6.6; % [MHz]
c0 = 1540; % [m/s]

a_nog = s2a(d_s, c0);
ka_nog = s2ka(d_s, fc);
bsc_nog = exp(-0.827* (2*pi*a_nog/c0)^2 * freq.^2 );

fprintf('-----ka (no g) -----\n')
fprintf('a  [µm]      = %.2f\n', a_nog );
fprintf('ka @ %.2fMHz = %.2f\n', fc, ka_nog );

%%%%%%%% "ka" with g %%%%%%%%
yy = log(bsc_gauss_reconstruct);
% yy = log(d_g) - d_s * freq.^2;

AA = -0.827 *(2*pi/c0)^2 * freq.^2;

x_opt = cgs(AA'*AA, AA'*yy); % AA \ yy

a_g = sqrt(x_opt);
ka_g = (2*pi*fc/c0)*a_g; 
bsc_g = exp(-0.827* (2*pi*a_g/c0)^2 * freq.^2 );

fprintf('-----ka (with g) -----\n')
fprintf('a  [µm]      = %.2f\n', a_g );
fprintf('ka @ %.2fMHz = %.2f\n', fc, ka_g );

% a method var minimization search
a_minVar = bsc_rpm_f0{2, iCount}.ESD_Uni(1)/2*1e6; % [um]
fprintf('-----a (minVar) -----\n')
fprintf('a  [µm]      = %.2f\n', a_minVar );

% ALL PLOT
line_width = 2;
font_size = 20;
figure, 
% set(gcf, 'Units', 'pixels', 'Position', [100, 100, 700, 700]); % [x, y, width, height] in pixels

semilogy(freq, bsc_rpm, 'k-', 'DisplayName', 'GT', 'LineWidth', line_width);
hold on; grid on
semilogy(freq, bsc_gauss_reconstruct, 'g--', 'DisplayName', 'Gauss BSC', 'LineWidth', line_width) ;
semilogy(freq, bsc_powlaw_reconstruct, 'b--', 'DisplayName', 'PowLaw BSC', 'LineWidth', line_width) ;

xlabel('Frequency [MHz]', 'FontSize', font_size);
ylabel('BSC [cm^{-1}\cdot sr^{-1}]', 'FontSize', font_size);
title(['BSC \sigma=', num2str(sigmaValue)])
legend('Location', 'southwest', 'FontSize', font_size);
hold off;
set(gca,'fontsize',16)

end