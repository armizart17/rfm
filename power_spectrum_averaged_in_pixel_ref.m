function AVERAGED_DSP = power_spectrum_averaged_in_pixel_ref(RF_PIXEL,P,SMOOTH, SPAN)

THOMPSON = 0;
RATIO = 0.50;
NB_SUBSPECTRUM = 20;

if THOMPSON == 1
     nb_sample = fix(RATIO*size(RF_PIXEL,1));
    starting_sample = ...
        1 : round((size(RF_PIXEL,1) - nb_sample)/NB_SUBSPECTRUM) : round((size(RF_PIXEL,1) - nb_sample));
    window_type = 2;
    wind = window_choice(nb_sample,window_type);
    mat_window = ( ones(size(RF_PIXEL,2), 1) * wind' )';
    sub_average_dsp = zeros(P/2,length(starting_sample));
    for j = 1 :  size(RF_PIXEL,3)
        for i = 1 : length(starting_sample)
            sig_wind = squeeze(RF_PIXEL((starting_sample(i):(starting_sample(i)+nb_sample-1)),:,j)).*mat_window;
            tf_sub_pixel = power(abs(fft(sig_wind,P,1)),2);
            sub_average_dsp(:,i) = mean(tf_sub_pixel(1:P/2,:),2);
        end
        AVERAGED_DSP_1(:,j) = mean(sub_average_dsp,2);
    end
    AVERAGED_DSP = mean(AVERAGED_DSP_1,2); 
%     figure(22)
%     plot(10*log10(AVERAGED_DSP_1(:,1)), 'k'); hold on;
%     plot(10*log10(AVERAGED_DSP_1(:,2)), 'b'); hold on;
%     plot(10*log10(AVERAGED_DSP_1(:,3)), 'm'); hold on;
%     plot(10*log10(AVERAGED_DSP_1(:,4)), 'k'); hold on;
%     plot(10*log10(AVERAGED_DSP), 'r', 'linewidth', 3); hold on;
%     plot(10*log10(AVERAGED_DSP_ini), '--k', 'linewidth', 3); hold off;
%     pause
    
elseif THOMPSON == 0
    tf_pixel = zeros(P, size(RF_PIXEL,2), size(RF_PIXEL,3));
    averaged = zeros(P/2, size(RF_PIXEL,3));
    AVERAGED_DSP = zeros(1, P/2);
    
    for i = 1:size(RF_PIXEL,3) % Iter for all frames
        tf_pixel(:,:,i) = power(abs(fft(squeeze(RF_PIXEL(:,:,i)),P,1)),2);
        averaged(:,i) = mean(squeeze(tf_pixel(1:P/2,:,i)),2);     
    end
    AVERAGED_DSP = mean(averaged,2); % Average for all frames
    
    % figure(25),
    % plot(10*log10(tf_pixel(:,44,i)), 'r', 'Linewidth', 2);hold on;
    % plot(10*log10(averaged(:,i)), 'k', 'Linewidth', 1);hold on;
    % plot(10*log10(AVERAGED_DSP), 'g', 'Linewidth', 0.5);hold off;
end

if SMOOTH == 1
    order = 1;
    AVERAGED_DSP = smooth(AVERAGED_DSP, SPAN, 'sgolay', order);
end