function spectralData = calc_powerSpectra(data, pars)
% function spectralData = calc_powerSpectra(data, pars)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION
% Provides the PowerSpectra (S) and compensates it accordingly
% -------------------------------------------------------------------------
% INPUTS
%       data:   Structure with the following fields:
%             - rf or RF        : RF data from sample phantom
%             - x               : Distance between lines in lateral axis [m]
%             - z               : Distance between lines in lateral axis [m]
%             - fs              : Sampling frequency [Hz]
%
%       pars:   Structure with the following fields:
%             - bw              : Bandwidth [MHz]
%             - overlap         : Fraction of overlap samples 
%             - blocksize       : datablock in wavelengths
%             - z_roi           : Roi oordinates of depth   [z_ini, z_end]
%             - x_roi           : Roi oordinates of lateral [x_ini, x_end]
%             - saran_layer     : True or false if needed it
%             - window_type     : (optional) Type of window (1) Hanning, (2) Tukey, (3) Hamming, (4) Tchebychev
% -------------------------------------------------------------------------
% OUTPUTS
%       spectralData:  Structure with the following fields:
%             - powerSpectra    : Power spectrum (rows, cols, freqChannels)
%             - depth           : depth coordinates for PowerSpectra   [m]
%             - lateral         : lateral coordinates for PowerSpectra [m]
%             - band            : band of frequencies [MHz]
%             - delta_snr       : Delta SNR
%             - rf_roi          : RF ROI
% -------------------------------------------------------------------------
% AUTHOR: EMZ based on LIM repository Chahuara's code compute_spectral_data.m
% Note: Remember PowerLaw 
% SR        = log(S_sam ./ S_ref); (rows, cols, freqChannels)
% SR_rpltv  = permute(SR, [3,1,2]); (freqChannels, rows, cols)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Reading experiment settings parameters
    bw              = pars.bw;
    overlap         = pars.overlap;
    blocksize_wv    = pars.blocksize;
    z_ini           = pars.z_roi(1);
    z_end           = pars.z_roi(2);
    x_ini           = pars.x_roi(1);
    x_end           = pars.x_roi(2);
    saran_layer     = pars.saran_layer;

    % Default window
    window_type     = 2; %  (1) Hanning, (2) Tukey, (3) Hamming, (4) Tchebychev
    if isfield(pars, 'window_type')
        window_type = pars.window_type;      
    end
    % position_roi    = pars.position_roi;
    
    % Reading phantoms parameters
    if isfield(data, 'rf')
        rfdata_sam   = data.rf;
    end
    if isfield(data, 'RF')
        rfdata_sam   = data.RF;
    end
    z            = data.z;
    x            = data.x;
    fs           = data.fs;

    dx = x(2)-x(1);
    dz = z(2)-z(1);
    
    c0 = 1540; % [m/s]
    lambda = c0/mean(bw)*1E-6;
    if isfield(pars, 'lambda')
        lambda   = pars.lambda;
    end
    
    fprintf('Power Spectra vHC');
    fprintf('\nCalculation bandwidth:\nBW = %g - %g MHz, Fc = %g MHz\nwavelength (wl) = %g m\n',bw(1),bw(2),mean([bw(1),bw(2)]),lambda);


    %%%%%%
    env_rfdata_sam = abs(hilbert(rfdata_sam));

    %
    ratio_zx = 1.0;
    nz = 2*round(blocksize_wv*lambda/dz /2 * ratio_zx); % Window size
    NFFT = 2^(nextpow2(nz/2)+2);

    % NFFT = 2^10;
    axis_f = (0:(NFFT/2-1))'*fs/NFFT;  % [Hz] (so 0-fs/2), as it should be

    freq_L = bw(1)*1E6; % [Hz] 
    freq_H = bw(2)*1E6; % [Hz]

    ind_f = axis_f >= freq_L & axis_f <= freq_H ;   
    band = axis_f(ind_f)*1E-6; % [MHz]
    bandFull = axis_f*1E-6; % [MHz]

    % f_ini = find(Xfreq >= bw(1), 1, 'first');
    % f_end = find(Xfreq <= bw(2), 1, 'last');
    % band_ind = Fi:Ff;
    % band  = axis_f(f_ini:f_end);  %band = axis_f(band_ind)

    % Saran Layer compensation 
    t_saran = 1;
    if (saran_layer)
        t_saran = transmission_saran(band);
    end

    %% ROI CROP v1

    % axis_x = x;   % [m]
    % axis_z = z;   % [m]

    % Transforming from distance to pixel values
    % position_roi(:,1) = 1 + (length(axis_x)-1)/(axis_x(end)-axis_x(1))*(position_roi(:,1)-axis_x(1));
    % position_roi(:,2) = 1 + (length(axis_z)-1)/(axis_z(end)-axis_z(1))*(position_roi(:,2)-axis_z(1));

    % roi_mask                        = poly2mask(position_roi(:,1),position_roi(:,2),size(rfdata_sam,1),size(rfdata_sam,2));
    % [dim_depth_roi,dim_lateral_roi] = find(roi_mask == 1);
    % 
    % axis_depth_roi   = axis_depth(min(dim_depth_roi):max(dim_depth_roi));
    % axis_lateral_roi = axis_lateral(min(dim_lateral_roi):max(dim_lateral_roi));
    % 
    % rfdata_sam_roi = rfdata_sam(min(dim_depth_roi):max(dim_depth_roi),min(dim_lateral_roi):max(dim_lateral_roi),:);
    % 
    % env_rfdata_sam_roi = env_rfdata_sam(min(dim_depth_roi):max(dim_depth_roi),min(dim_lateral_roi):max(dim_lateral_roi),:);

    %% ROI CROP V2
    ind_x = x_ini <= x & x <= x_end;
    ind_z = z_ini <= z & z <= z_end;

    axis_depth_roi = z(ind_z);
    axis_lateral_roi = x(ind_x);

    rfdata_sam_roi = rfdata_sam(ind_z, ind_x);
    env_rfdata_sam_roi = env_rfdata_sam(ind_z, ind_x);

    %%
    [I2,J2,~] = size(rfdata_sam_roi);

    nb_lambda_axial     = ratio_zx*blocksize_wv; % wavelengths (slightly bigger)
    axial_gate_length   = nb_lambda_axial*lambda;              % [m]
    nb_sample_roi_axial = floor(axial_gate_length/dz);    % [sample]

    % Axial overlap [samples]
    nb_sample_overlap_axial = fix(nb_sample_roi_axial*overlap); % [sample]
    nb_roi_axial            = ceil((I2-nb_sample_overlap_axial)/(nb_sample_roi_axial-nb_sample_overlap_axial));

    % Lateral windows length
    nb_lambda_lateral   = blocksize_wv; % wavelengths
    lateral_gate_length = nb_lambda_lateral*lambda; % [m]

    nb_line_lateral     = ceil(nb_lambda_lateral*lambda/dx);

    % idx_lat = find((axis_lateral_roi-axis_lateral_roi(1)) >= lateral_gate_length);
    % nb_line_lateral = idx_lat(1); clear idx_lat;
                                     
    % Recouvrement interfenetre lateral
    nb_sample_overlap_lateral = fix(nb_line_lateral*overlap); % [sample]
    nb_roi_lateral = fix((J2-nb_sample_overlap_lateral)/(nb_line_lateral-nb_sample_overlap_lateral));

    fprintf('\nAxial gate length = %d wl = %g m = %g samp.\nNumber of axial region = %g\n',nb_lambda_axial,axial_gate_length,nb_sample_roi_axial,nb_roi_axial);
    fprintf('\nLateral gate length = %2.1f wl = %g m = %g lines\nNumber of lateral region = %g\n\n',nb_lambda_lateral,lateral_gate_length,nb_line_lateral,nb_roi_lateral);


    % Initialization
    X_ROI = zeros(nb_roi_lateral, 1);
    Z_ROI = zeros(nb_roi_axial, 1);
    tf_mat_sam = zeros(nb_roi_axial, nb_roi_lateral, length(band));
    S_full = zeros(nb_roi_axial, nb_roi_lateral, length(bandFull));

    X_ROI_v2 = zeros(nb_roi_lateral, 1);
    Z_ROI_v2 = zeros(nb_roi_axial, 1);

    delta_snr = zeros(nb_roi_axial, nb_roi_lateral);
    SNRopt = sqrt(1/(4/pi - 1)); % 1.91 


    for i = 1:nb_roi_axial

        block_i = (i-1)*(nb_sample_roi_axial-nb_sample_overlap_axial) + (1:nb_sample_roi_axial);
        window_vec = window_choice(length(block_i),window_type);

        if i == nb_roi_axial
            u = 0;
            while block_i(end) > I2
                u = u + 1;
                clear block_i;
                clear mat_window;
                clear wind;
                block_i = (i-1)*(nb_sample_roi_axial-nb_sample_overlap_axial) + (1:(nb_sample_roi_axial-u));
                warning('\nFor i = %d, ROI final sample (%d) is larger than the RF-signal size (%d)\n',i,block_i(end),I2);
            end
            window_vec = window_choice(length(block_i),window_type);
        end

        for j = 1:nb_roi_lateral

            block_j = (j-1)*(nb_line_lateral-nb_sample_overlap_lateral) + (1:nb_line_lateral);
            window_mat = window_vec*ones(1,length(block_j));

            if j == nb_roi_lateral
                w = 0;
                while block_j(end) > J2
                    w = w + 1;
                    clear block_j;
                    block_j = (j-1)*(nb_line_lateral-nb_sample_overlap_lateral) + (1:nb_line_lateral-w);
                    warning('\nFor j = %d, ROI final sample (%d) is larger than the ROW size (%d)\n',j,block_j(end),J2);
                end
                window_mat = window_vec*ones(1,length(block_j));
            end

            Z_ROI(i) = axis_depth_roi(round((block_i(end)+block_i(1))/2));   % [m]
            X_ROI(j) = axis_lateral_roi(round((block_j(end)+block_j(1))/2)); % [m]

            Z_ROI_v2(i) = axis_depth_roi(round(block_i(1)));   % [m]
            X_ROI_v2(j) = axis_lateral_roi(round(block_j(1)));   % [m]

            % Selection of data in data-block
            rfdata_block  = rfdata_sam_roi(block_i,block_j,:);

            % Remove the continuous component
            rfdata_block = remove_DC_in_pixel_gen(squeeze(mean(rfdata_block)), rfdata_block, window_mat);

            % Averaged power spectra
            tf_block  = power_spectrum_averaged_in_pixel_gen(rfdata_block, NFFT, 0, 0);
           
            % Selection of the useful bandwidth
            tf_mat_sam(i,j,:) = tf_block(ind_f) ./ t_saran;
            
            S_full(i,j,:) = tf_block;

            %%%%%%%%%
            env_data_block_sam  = env_rfdata_sam_roi(block_i,block_j,:);
            SNR = mean(env_data_block_sam(:))/std(env_data_block_sam(:));
  
            delta_snr(i,j) = 100*abs(SNR - SNRopt)/SNRopt;

        end

    end

    % Packing
    spectralData.powerSpectra      = tf_mat_sam;
    spectralData.depth             = Z_ROI(:); % € old 
    spectralData.lateral           = X_ROI(:); % € old 
    
    % spectralData.depth2            = Z_ROI_v2;
    % spectralData.lateral2          = X_ROI_v2;

    spectralData.band              = band(:);
    spectralData.axial_gate_length = axial_gate_length;
    spectralData.delta_snr         = delta_snr;

    spectralData.rf_roi            = rfdata_sam_roi;
    spectralData.z_roi             = axis_depth_roi;
    spectralData.x_roi             = axis_lateral_roi;

    spectralData.Sfull             = mean(S_full, 4);
    spectralData.bandFull          = bandFull;

return

%% POWER SPECTRUM
function AVERAGED_DSP = power_spectrum_averaged_in_pixel_gen(RF_PIXEL, P, SMOOTH, SPAN)

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
        
    else
        tf_pixel = zeros(P, size(RF_PIXEL,2), size(RF_PIXEL,3));
        averaged = zeros(P/2, size(RF_PIXEL,3));
        % AVERAGED_DSP = zeros(1, P/2);
        
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
return

%% REMOVE DC LEVEL
function RF_PIXEL_PROCESSED = remove_DC_in_pixel_gen(MEAN_PIXEL, RF_PIXEL, WINDOW)

    RF_PIXEL_PROCESSED = zeros(size(RF_PIXEL,1), size(RF_PIXEL,2), size(RF_PIXEL,3));
    
    for i = 1 : size(RF_PIXEL,3)
        O   = ones(1,size(RF_PIXEL,1));
        MAT_MEAN = O'*squeeze(MEAN_PIXEL(:,i)');
        RF_PIXEL_PROCESSED(:,:,i) = (RF_PIXEL(:,:,i) - MAT_MEAN).*WINDOW;
    end

return
%% WINDOW CHOICE (OPTIONAL)
% function window = window_choice(len_window, choice)
% 
%     switch choice
%         case 1
%             window = hann(len_window);
%         case 2
%             window = tukeywin(len_window,0.25);
%         case 3
%             window = hamming(len_window);
%         case 4
%             window = chebwin(len_window);
%     end
% 
% return

%% TRANSMISSION SARAN (OPTIONAL)
% function t = transmission_saran(band)
% 
%     c_s = 2400; %[m/s]
%     c_p = 1540; %[m/s]
%     rho_s = 1690; %[kg/m^3]
%     rho_p = 1000; %[kg/m^3]
%     Z_p = c_p*rho_p;
%     Z_s = c_s*rho_s;
% 
%     f0 = band*1e6;  %[Hz]
%     alpha = 5*(band).^1.5; %[Np/MHz^1.5/m]
%     L = 25.1e-6; %[m]
% 
% 
%     t = abs(     2*Z_p ./...
%         (    (2*Z_p)*cos((2*pi*f0/c_s - 1i*alpha)*L) ...
%         + 1i *(Z_s+Z_p^2/Z_s)*sin((2*pi*f0/c_s - 1i*alpha)*L) )) .^4;
% 
% return
