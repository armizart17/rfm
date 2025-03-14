% Simple code to check BSC of RPM and fit Gaussian and PowLaw
% In order to see parameters of each 
% Option 1: Polyfit
% MY GT is bsc_rpm = BSC.BSCcurve_Uni(:,2); % median


%% GAUSS SIMULATION 1 sam 0.7, ref 0.5

% -----RPM Gauss (g.exp(-s.f^2))-----
% d_g          = 0.940154
% g_s/g_r [dB] = -0.268010
% Δs           = 0.022079
% ---------
% -----RPM POWER LAW (b.(f^n))-----
% d_b          = 4.506053
% b_s/b_r [dB] = 6.537963
% Δn           = -1.383433

%% GAUSS SIMULATION 1.1 sam 0.5, ref 0.5

% -----RPM Gauss (g.exp(-s.f^2))-----
% d_g          = 0.953247
% g_s/g_r [dB] = -0.207946
% Δs           = 0.022739
% ---------
% -----RPM POWER LAW (b.(f^n))-----
% d_b          = 4.786599
% b_s/b_r [dB] = 6.800270
% Δn           = -1.424616

%% GAUSS SIMULATION 2 sam 0.5 ref 0.5 by EMZ
% \sigmaGaussKernel = 1
% -----RPM Gauss (g.exp(-s.f^2))-----
% Δg           = 0.877044
% g_s/g_r [dB] = -0.569784
% Δs           = 0.091409
% ---------
% -----RPM POWER LAW (b.(f^n))-----
% Δb           = 598.333871
% b_s/b_r [dB] = 27.769436
% Δn           = -5.749143
%% BSC RPM GT
BSC = calculateBSC_RPM_ok(SAM, REF, pars); % fast TBD**
% BSC = calculateBSC_RPM_fast(SAM, REF, pars); % fast TBD**
bsc_rpm = BSC.BSCcurve_Uni(:,2); % median
freq = BSC.band(:);

%% CODE
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

c0 = 1540; fc = 6;
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
a_minVar = BSC.ESD_Uni(1)/2*1e6; % [um]
fprintf('-----a (minVar) -----\n')
fprintf('a  [µm]      = %.2f\n', a_minVar );

%% ALL PLOT

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
title('BSC')
legend('Location', 'southwest', 'FontSize', font_size);
hold off;
set(gca,'fontsize',16)
%% 