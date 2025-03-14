function analyze_filter(mask, title_str, dx)
    % Function to compute and visualize autocorrelation and FFT of a 2D filter

    figure; imagesc(mask);
    colorbar; title(['Mask: ', title_str]);

    % Compute autocorrelation
    Rmask = conv2(mask, mask, 'same');

    figure; imagesc(Rmask);
    colorbar; title(['Autocorrelation: ', title_str]);

    % Compute 2D FFT
    NFFT = 2^10;
    SRmask = fftshift(abs(fft2(Rmask, NFFT, NFFT)));

    figure; imagesc(20*log10(SRmask));
    colorbar; title(['FFT(Autocorrelation): ', title_str]);

    % Compute spatial frequency
    ft_space = (1540/2) / dx;
    freq = ((-NFFT/2):(NFFT/2-1)) / NFFT * ft_space * 1e-6;

    figure; plot(freq, 10*log10(SRmask(:, NFFT/2)), 'k');
    xlabel('Frequency (MHz)');
    ylabel('Power (dB)');
    title(['Frequency Response: ', title_str]);

    % Log-slope estimation
    log_bsc = log(SRmask(:,NFFT/2));
    fL = 3;
    fH = 7.5;
    f0 = mean([fL fH]);

    % [~, idx_fL] = min(abs(freq - fL));
    % [~, idx_fH] = min(abs(freq - fH));

    idx_fL = find(freq >= fL, 1, 'first');
    idx_fH = find(freq >= fH, 1, 'first');

    new_freq = freq(idx_fL:idx_fH);
    new_log_bsc = log_bsc(idx_fL:idx_fH);
    slope = polyfit(new_freq, new_log_bsc, 1);

    n = f0 * slope(1);
    disp(['Slope (n/f0) for ', title_str, ': ', num2str(n)]);
end

