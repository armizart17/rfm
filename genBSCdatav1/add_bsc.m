function [rf_output] = add_bsc(rf_input, b, n, fs)
% function [rf_output] = add_bsc(rf_input, b, n, fs)

% rf_output = zeros(size(rf_input));
[M, N] = size(rf_input);



for ii = 1:N

    disp(ii)
    rf = rf_input(:,ii);

    %figure; plot(rf);

    NFFT = 2^15;
    freq = (-NFFT/2:NFFT/2-1)/NFFT*fs*1e-6;

    ps = power(abs(fftshift(fft(rf,NFFT))),2);

    ps_new = ps.*(b*abs(freq').^n);
    ps_new = ps_new/max(ps_new);

    newfft = fftshift(fft(rf,NFFT)) .* sqrt(b*abs(freq').^n);
    newps = power( abs(newfft),2 );

    newrf = ifft(fftshift(newfft),NFFT);

    % rf_output(:,ii) = newrf; % old COILA
    rf_output(:,ii) = newrf(1:M); % €MZ


    continue
    figure,
    plot(freq, 20*log10(ps/max(ps))); hold on;
    plot(freq, 10*log10(ps_new) );
    plot(freq, 10*log10(b*abs(freq').^n))
    plot(freq, 10*log10(newps./max(newps)))

   
    figure,
    plot(rf); hold on;
    plot(newrf);
    %keyboard



    keyboard

end



end