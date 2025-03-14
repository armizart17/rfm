% Test Modulation 2D
% Parameters
sigma = 1.25;          % Standard deviation of the Gaussian
f0 = 50;            % Modulation frequency
size = 9;         % Size of the grid
theta = deg2rad(90);          % Orientation angle in radians

% Create a 2D grid
x = linspace(-size/2, size/2, size);
y = linspace(-size/2, size/2, size);
[X, Y] = meshgrid(x, y);

% 2D Gaussian
gaussian = exp(-(X.^2 + Y.^2) / (2 * sigma^2)) / (2 * pi * sigma^2);

% Cosine modulation
modulation = 2*cos(2 * pi * f0 * (X * cos(theta) + Y * sin(theta)));

% Gaussian-modulated cosine (Gabor filter)
gabor_filter = gaussian .* modulation;
% gabor_filter = gaussian;

% Plotting
figure;

% Gaussian Filter
subplot(1, 3, 1);
imagesc(x, y, gaussian);
axis image;
colorbar;
title('2D Gaussian Filter');

% Cosine Modulation
subplot(1, 3, 2);
imagesc(x, y, modulation);
axis image;
colorbar;
title('Cosine Modulation');

% Gaussian Modulated Cosine (Gabor)
subplot(1, 3, 3);
imagesc(x, y, gabor_filter);
axis image;
colorbar;
title('2D Gabor Filter');

%

colors = lines(length(f0));  % Use MATLAB's built-in distinct colors
 
currentColor = colors;
% MASK 

mask = gabor_filter;
NFFT = 512;

SRmask_v1 = abs(fftshift(fft2(mask, NFFT, NFFT))).^2;
ps_gauss2d_v1 = SRmask_v1(:, NFFT/2);

Rmask = conv2(mask, mask, 'same');
SRmask_v2 = abs(fftshift(fft2(Rmask, NFFT, NFFT)));
ps_gauss2d_v2 = SRmask_v2(:, NFFT/2);

figure, 
subplot(121),
imagesc(SRmask_v1), title('V1')

subplot(122), 
imagesc(SRmask_v2), title('V2')

% THEORETICAL
c0 = 1540;
dx = 0.3e-3 /16;
ft_space = c0/2/dx;
NFFT = 512;

freq_s = ((-NFFT/2):(NFFT/2-1))/NFFT*ft_space; %[Hz]
freq_MHz = freq_s*1e-6;

figure(77)
% --- Plot Mask Spectrum v1---
plot(freq_MHz, 10*log10(ps_gauss2d_v1), '--', 'Color', currentColor, 'LineWidth', 1.5, ...
        'DisplayName', ['Maskv1, f_0=', num2str(f0)]);
hold on;

% --- Plot Mask Spectrum v2 ---
plot(freq_MHz, 10*log10(ps_gauss2d_v1), '--', 'Color', currentColor, 'LineWidth', 1.5, ...
        'DisplayName', ['Maskv2 f_0=', num2str(f0)]);
hold on;
legend('Location', 'Best');
xlabel('Frequency [MHz]');
ylabel('Power [dB]');
title('Power Spectra Comparison for Different \sigma Values');
grid on;

figure(78)
% --- Plot Mask Spectrum v1---
plot(freq_MHz, (ps_gauss2d_v1), '--', 'Color', currentColor, 'LineWidth', 1.5, ...
        'DisplayName', ['Maskv1, f_0=', num2str(f0)]);
hold on;

% --- Plot Mask Spectrum v2 ---
plot(freq_MHz, (ps_gauss2d_v1), '--', 'Color', currentColor, 'LineWidth', 1.5, ...
        'DisplayName', ['Maskv2 f_0=', num2str(f0)]);
hold on;
legend('Location', 'Best');
xlabel('Frequency [MHz]');
ylabel('Power [dB]');
title('Power Spectra Comparison for Different \sigma Values');
grid on;

%%
% Parameters
sigma = 1.25;              % Standard deviation of the Gaussian

dx = 0.3e-3 /16;
c0 = 1540;
f_shift = 10e6;
% f_shift = 6.1600e+09;
f0 = 2*f_shift*dx/c0;
% f0 = 15;                % Modulation frequency
theta1 = 0;          % First direction (45 degrees)
theta2 = 90;        % Second direction (-45 degrees)
size = 9;             % Size of the grid

% Create a 2D grid
x = linspace(-size/2, size/2, size);
y = linspace(-size/2, size/2, size);
[X, Y] = meshgrid(x, y);

% 2D Gaussian
gaussian = exp(-(X.^2 + Y.^2) / (2 * sigma^2)) / (2 * pi * sigma^2);

% First modulation (theta1 degrees)
modulation1 = cos(2 * pi * f0 * (X * cos(deg2rad(theta1)) + Y * sin(deg2rad(theta1))));

% Second modulation (theta2 degrees)
modulation2 = cos(2 * pi * f0 * (X * cos(deg2rad(theta2)) + Y * sin(deg2rad(theta2))));

% Add the two modulations
combined_modulation = modulation1 + modulation2;

% Apply Gaussian envelope
gabor_combined = gaussian .* combined_modulation;

% Plotting the result
figure;

% Combined Gabor Filter
imagesc(x, y, gabor_combined);
axis image;
colorbar;
title(sprintf('2D Gabor Filter (%d° and %d° Modulation)', theta1, theta2));


mask = gabor_combined;
NFFT = 512;

SRmask_v1 = abs(fftshift(fft2(mask, NFFT, NFFT))).^2;
ps_gauss2d_v1 = SRmask_v1(:, NFFT/2);

Rmask = conv2(mask, mask, 'same');
SRmask_v2 = abs(fftshift(fft2(Rmask, NFFT, NFFT)));
ps_gauss2d_v2 = SRmask_v2(:, NFFT/2);

figure, 
sgtitle('2D Power Spectrum')
subplot(121),
imagesc( (SRmask_v1/max(SRmask_v1(:))) ),
title('Direct FFT^2')

subplot(122), 
imagesc( (SRmask_v2/max(SRmask_v2(:))) ),
title('Autocrr + FFT')

% THEORETICAL
c0 = 1540;
dx = 0.3e-3 /16;
ft_space = c0/2/dx;
NFFT = 512;

freq_s = ((-NFFT/2):(NFFT/2-1))/NFFT*ft_space; %[Hz]
freq_MHz = freq_s*1e-6;

figure,
colors = lines(length(f0));  % Use MATLAB's built-in distinct colors
 
currentColor = colors;
% --- Plot Mask Spectrum v1---
plot(freq_MHz, 10*log10(ps_gauss2d_v1), '--', 'Color', currentColor, 'LineWidth', 1.5, ...
        'DisplayName', ['Maskv1, f_0=', num2str(f0)]);
hold on;

% --- Plot Mask Spectrum v2 ---
plot(freq_MHz, 10*log10(ps_gauss2d_v1), ':', 'Color', currentColor, 'LineWidth', 1.5, ...
        'DisplayName', ['Maskv2 f_0=', num2str(f0)]);
hold on;
legend('Location', 'Best');
xlabel('Frequency [MHz]');
ylabel('Power [dB]');
title('Power Spectra for Different f_0 Values');
grid on;

figure,
% --- Plot Mask Spectrum v1---
plot(freq_MHz, (ps_gauss2d_v1), '--', 'Color', currentColor, 'LineWidth', 1.5, ...
        'DisplayName', ['Maskv1, f_0=', num2str(f0)]);
hold on;

% --- Plot Mask Spectrum v2 ---
plot(freq_MHz, (ps_gauss2d_v1), '--', 'Color', currentColor, 'LineWidth', 1.5, ...
        'DisplayName', ['Maskv2 f_0=', num2str(f0)]);
hold on;
legend('Location', 'Best');
xlabel('Frequency [MHz]');
ylabel('Power');
title('Power Spectra for Different f_0 Values');
grid on;


