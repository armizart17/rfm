function [rf_output] = add_bsc_v2(rf_input, b_vals, n_vals, fs, geometry, template)
    % Add BSC with different parameters for regions defined by geometry
    % rf_input: Input 2D RF data
    % b_vals, n_vals: Arrays of BSC parameters [b1, b2] and [n1, n2]
    % fs: Sampling frequency
    % geometry: Type of geometry ('two_layer' or 'circle')

    % Initialize output
    rf_output = zeros(size(rf_input));
    [M, N] = size(rf_input);

    % Create mask based on geometry
    mask = zeros(M, N);
    if strcmp(geometry, 'two_layer')
        % Divide the area into two vertical layers
        mask(:, 1:round(N/2)) = 1;  % First half
        mask(:, round(N/2)+1:end) = 0;  % Second half
    elseif strcmp(geometry, 'circle')
        % Create a circular mask in the center
        % [X, Y] = meshgrid(1:N, 1:M);
        % center_x = N / 2;
        % center_y = M / 2;
        % radius = min(M, N) / 4;
        % mask(((X - center_x).^2 + (Y - center_y).^2) <= radius^2) = 1;
        % mask(mask == 0) = 2;  % Outside circle
    
        % temporal for simulation
        mask = template;
    else
        error('Unsupported geometry type');
    end

    % Loop through each column (RF line)
    for ii = 1:N
        rf = rf_input(:,ii);
        
        % FFT parameters
        NFFT = 2^15;
        freq = (-NFFT/2:NFFT/2-1)/NFFT*fs*1e-6;

        % Get current mask region (1 or 2)
        region = mask(:, ii);

        % Apply BSC based on region
        if all(region == 1)
            b = b_vals(1);
            n = n_vals(1);
        else
            b = b_vals(2);
            n = n_vals(2);
        end

        % Original power spectrum
        ps = power(abs(fftshift(fft(rf, NFFT))), 2);

        % Modify power spectrum
        ps_new = ps .* (b * abs(freq').^n);
        ps_new = ps_new / max(ps_new);

        % Modify signal in frequency domain
        newfft = fftshift(fft(rf, NFFT)) .* sqrt(b * abs(freq').^n);
        newps = power(abs(newfft), 2);

        % Inverse FFT to get modified RF signal
        newrf = ifft(fftshift(newfft), NFFT);
        rf_output(:, ii) = newrf(1:M);  % Ensure size matches input

        % Optional plots for debugging
        % if ii == fix(N/2) % at half
        %     figure, plot(freq, 20*log10(ps / max(ps))); hold on;
        %     plot(freq, 10*log10(ps_new));
        %     plot(freq, 10*log10(b * abs(freq').^n));
        %     plot(freq, 10*log10(newps / max(newps)));
        %     title('Power Spectrum Analysis');
        % 
        %     figure, plot(rf); hold on;
        %     plot(real(newrf(1:M)));
        %     title('RF Signal Comparison');
        % end
    end
end
