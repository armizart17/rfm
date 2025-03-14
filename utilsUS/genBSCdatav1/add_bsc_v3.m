function [rf_output] = add_bsc_v3(rf_input, fs, geometry, bsc_model, coeff_1, coeff_2)
% function [rf_output] = add_bsc_v3(rf_input, fs, geometry, bsc_model, coeff_1, coeff_2)
% Add BSC or Gaussian model modifications to RF data line per line
%
% rf_input: Input 2D RF data
% fs: Sampling frequency
% geometry: Type of geometry ('two_layer' or 'homo')
% bsc_model: Model type 'power-law' (b.f^n) or 'gauss'(g.exp(-s.f^2))
% coeff_1, coeff_2: Model coefficients
%   For 'power-law': coeff_1 = [b_left, b_right], coeff_2 = [n_left, n_right]
%   For 'gauss'    : coeff_1 = [g_left, g_right], coeff_2 = [s_left, s_right]
% Author: EMZ, based on AC's code

    % Initialize output
    rf_output   = zeros(size(rf_input));
    [M, N]      = size(rf_input);

    % Create mask based on geometry
    mask = zeros(1, N);
    if strcmp(geometry, 'two_layer')
        % Divide the area into two vertical layers
        mask(1:round(N/2))     = 1;  % First half
        mask(round(N/2)+1:end) = 2;  % Second half
    elseif strcmp(geometry, 'homo')
        % Apply same modification to all lines
        mask(:) = 1;
    else
        error('Unsupported geometry type');
    end

    % FFT parameters
    NFFT = 2^15;
    freq = (-NFFT/2:NFFT/2-1)/NFFT * fs * 1e-6;  % Frequency in MHz

    % Loop through each column (RF line)
    for ii = 1:N
        rf = rf_input(:, ii);

        % Determine coefficients based on mask and geometry
        if mask(ii) == 1
            c1 = coeff_1(1);
            c2 = coeff_2(1);
        else
            c1 = coeff_1(2);
            c2 = coeff_2(2);
        end

        % Original power spectrum
        ps = abs(fftshift(fft(rf, NFFT))).^2;

        if strcmp(bsc_model, 'power-law')
            % Modify power spectrum using power-law model
            ps_new = ps .* (c1 * abs(freq').^c2);
            ps_new = ps_new / max(ps_new);
            % Modify signal in frequency domain
            newfft = fftshift(fft(rf, NFFT)) .* sqrt(c1 * abs(freq').^c2);

        elseif strcmp(bsc_model, 'gauss')
            % Modify power spectrum using Gaussian model
            ps_new = ps .* (c1 * exp(c2 * freq'.^2));
            ps_new = ps_new / max(ps_new);
            % Modify signal in frequency domain
            newfft = fftshift(fft(rf, NFFT)) .* sqrt(c1 * exp(c2 * freq'.^2));

        else
            error('Unsupported BSC model type');
        end

        % Inverse FFT to get modified RF signal
        newrf               = ifft(fftshift(newfft), NFFT);
        rf_output(:, ii)    = newrf(1:M);  % Ensure size matches input

        % Optional debugging plot for one column
        % if ii == fix(N/2) % at half
        %     figure, plot(freq, 20*log10(ps / max(ps))); hold on;
        %     plot(freq, 10*log10(ps_new));
        %     title('Power Spectrum Analysis');
        %     legend('Original', 'Modified');
        % 
        %     figure, plot(rf); hold on;
        %     plot(real(newrf(1:M)));
        %     title('RF Signal Comparison');
        %     legend('Original', 'Modified');
        % end
    end
end
