function [BSC] = calculateBSC_RPM_fast(DATA, REF, pars)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fast version of calculateBSC_RPM
% This version only computes and stores the required BSC parameters,
% eliminating unnecessary calculations for improved efficiency.
% 
% INPUTS: 
% - DATA : Sample data structure
% - REF  : Reference phantom structure
% - pars : Parameters structure
% 
% OUTPUTS: 
% - BSC: Object containing selected BSC data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Read parameters
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

    
    P               = pars.P; %for number of points NFFT
    nb_lambda_axial = blocksize_wv; % [wavelength]   
    overlap_axial   = overlap;   
    overlap_lateral = overlap;  

  
    % Reading data
    z               = DATA.z;
    x               = DATA.x;
    fs              = DATA.fs;

    if isfield(DATA, 'rf')
        rfdata_sam   = DATA.rf;
    end
    if isfield(DATA, 'RF')
        rfdata_sam   = DATA.RF;
    end
    if isfield(REF, 'rf')
        rfdata_ref   = REF.rf;
    end
    if isfield(REF, 'RF')
        rfdata_ref   = DATA.RF;
    end
    m_sam = 1;
    m_ref = 1;
    if isfield(DATA, 'alpha_power')
        m_sam   = DATA.alpha_power;
    end
    if isfield(REF, 'alpha_power')
        m_ref   = REF.alpha_power;
    end

    bmode_sam = db(abs(hilbert(rfdata_sam)));
    bmode_sam = bmode_sam - max(bmode_sam(:));

   
%     latlines = 256; % The number of lines in each bmode_sam
%     transducer_footprint = 38; % The transducer footprint in mm
%     delta_lat = transducer_footprint/latlines; % The distance in between adjacent lines
    
    % Begin the estimation process frame per frame
    % nb_frame = size(rfdata_sam, 3);

tic;
    % Define frequency axis
    c0 = 1540; % [m/s]
    fct = 1e6; %for MHz                         % [MHz]
    lambda = c0/(mean([bw(1) bw(2)])*1e6);% [wavelength]
    fprintf('\nCalculation bandwidth:\nBW = %g - %g MHz, Fc = %g MHz\nwavelength (wl) = %g mm\n', bw(1), bw(2), mean([bw(1) bw(2)]), lambda*1e3);        
    

    dz = z(2)-z(1);
    ratio_zx = 1;
    nz = 2*round(blocksize_wv*lambda/dz /2 * ratio_zx); % Window size
    NFFT = 2^(nextpow2(nz/1)+1);
    P = NFFT;

    Xfreq = (0 : (P/2-1))*fs/P/fct;    % [MHz] 
    Fi       = find(Xfreq >= bw(1), 1, 'first');
    Ff       = find(Xfreq <= bw(2), 1, 'last');
    band_ind = Fi:Ff;
    % band_ind = Xfreq >= bw(1) & Xfreq <= bw(2); find(Xfreq >= bw(1) & Xfreq <= bw(2)) 
    band     = Xfreq(band_ind);       % [MHz]
    band_k   = 2*pi*band*1e6/c0;   % [rad]

    T = transmission_saran(band);
    
    % Crop ROI
    ind_x = x_ini <= x & x <= x_end; 
    ind_z = z_ini <= z & z <= z_end;

    Zaxis = z(ind_z);
    Xaxis = x(ind_x);

    DATA_ROI_seg = rfdata_sam(ind_z, ind_x);
    PHAN_ROI_seg = rfdata_ref(ind_z, ind_x);
    [I2, J2] = size(DATA_ROI_seg);
    
    % Define axial and lateral windowing
    axial_gate_length = blocksize_wv * lambda;
    delta_z = z(2) - z(1);
    nb_sample_ROI_axial = floor(axial_gate_length/delta_z);
    nb_sample_overlap_axial = fix(nb_sample_ROI_axial * overlap);
    nb_ROI_axial = fix((I2-nb_sample_overlap_axial)/(nb_sample_ROI_axial - nb_sample_overlap_axial));
    
    nb_lambda_lateral = blocksize_wv;
    lateral_gate_length = nb_lambda_lateral * lambda;
    nb_line_ROI_lateral = find((Xaxis-Xaxis(1)) >= lateral_gate_length, 1);
    nb_sample_overlap_lateral = fix(nb_line_ROI_lateral * overlap);
    nb_ROI_lateral = fix((J2-nb_sample_overlap_lateral)/(nb_line_ROI_lateral - nb_sample_overlap_lateral));
    
    % Initialize output matrices
    BSC_w_uniform_map = nan(nb_ROI_axial, nb_ROI_lateral, length(band));
    
    for i = 1:nb_ROI_axial
        pixi = (i-1)*(nb_sample_ROI_axial - nb_sample_overlap_axial) + (1:nb_sample_ROI_axial-1);
        wind = window_choice(length(pixi),window_type);
        if i == nb_ROI_axial
            u = 0;
            while pixi(end) > I2
                u = u + 1; clear pixi mat_window wind;
                pixi = (i-1)*(nb_sample_ROI_axial - nb_sample_overlap_axial) ...
                    + (1:(nb_sample_ROI_axial-u));
                warning('\nFor i = %d, ROI final sample (%d) is larger than the RF-signal size (%d)\n', ...
                    i, pixi(end), I2);
            end
            wind = window_choice(length(pixi),window_type);
        end
        for j = 1:nb_ROI_lateral
            pixj = (j-1)*(nb_line_ROI_lateral - nb_sample_overlap_lateral) + (1:nb_line_ROI_lateral);
            mat_window = (ones(length(pixj),1)*wind')';
            if j == nb_ROI_lateral
                w = 0;
                while pixj(end) > J2
                    w = w + 1; clear pixj;
                    pixj = (j-1)*(nb_line_ROI_lateral - nb_sample_overlap_lateral) + (1:nb_line_ROI_lateral-w);
                    warning('\nFor j = %d, ROI final sample (%d) is larger than the ROW size (%d)\n', ...
                        j, pixj(end), J2);
                end
                mat_window = (ones(length(pixj),1)*wind')';
            end


            Z_ROI(i) = Zaxis( round((pixi(end) + pixi(1))/2 ) )*1e2; % [cm]
            X_ROI(j) = Xaxis( round((pixj(end) + pixj(1))/2 ) )*1e2; % [cm]

            % Selection of data in data-block
            phan_block = PHAN_ROI_seg(pixi, pixj,:);
            data_block = DATA_ROI_seg(pixi, pixj);
            
            there_is_nan = double(isnan(data_block));
            [z_find_nan,x_find_nan] = find(there_is_nan == 1);
            
            if  isempty(z_find_nan) == 0 && isempty(x_find_nan) == 0
                
            elseif (isempty(z_find_nan) == 1 && isempty(x_find_nan) == 1)

                % Remove the continuous component
                phan_block = remove_DC_in_pixel_gen(squeeze(mean(phan_block)), phan_block, mat_window);
                data_block = remove_DC_in_pixel_gen(mean(data_block), data_block, mat_window);
                
                % Averaged power spectra
                tf_phan = power_spectrum_averaged_in_pixel_gen(phan_block,P,0,0);
                tf_data = power_spectrum_averaged_in_pixel_gen(data_block,P,0,0);
                
                aux_phan = 10*log10(tf_phan/(max(tf_phan)));
                aux_data = 10*log10(tf_data/(max(tf_data)));

                % if (i==1 && j==1)
                % 
                %     figure (111),
                %     subplot (211),
                %     plot(Xfreq, aux_data' ), hold on, grid minor;
                %     plot(Xfreq, aux_phan' ); hold off;
                %     xlabel('Frequency [MHz]'), ylabel('Power Spectrum INI [dB]')
                %     legend('Sample', 'Ref');
                %     axis([0 max(Xfreq) -60 0]); %xlim 0-20MHZ y ylim -60dB 0dB
                %     title ('Initial Power Spectrum (min depth)')
                % 
                % end
                % 
                % if (i==nb_ROI_axial && j == nb_ROI_lateral)
                % 
                %     figure (111),
                %     subplot (212),
                %     plot(Xfreq, aux_data' ), hold on, grid minor;
                %     plot(Xfreq, aux_phan' ); hold off;
                %     xlabel('Frequency [MHz]'), ylabel('Power Spectrum END [dB]')
                %     legend('Sample', 'Ref');
                %     axis([0 max(Xfreq) -60 0]); %xlim 0-20MHZ y ylim -60dB 0dB    
                %     title ('Final Power Spectrum (max depth)')
                % 
                % end

                 % Selection of the useful bandwidth
                
                %tf_phan = tf_phan(band_ind)./T(:); % If both phantoms have Saran wrap, T = 1
                tf_phan = tf_phan(band_ind);
                tf_data = tf_data(band_ind);
                
                % Spectral ratio for BSC estimation
                SR(i,j,:) = tf_data ./ tf_phan;
                
                % Compensation of the attenuation
                X_ACS_ind = find(1e2*x >= X_ROI(j));
                Z_ACS_ind = find(1e2*z >= Z_ROI(i));

                 %% Cumulative attenution for BSC compensation for uniform media
                % See Eq. (5) in 10.1016/j.ultras.2021.106376
                %Z_ROI(i) = Zaxis( round((pixi(end) + pixi(1))/2 ) )*1e2; % [cm]
                Z_ROI_ini(i) = Zaxis( pixi(1) ) *1e2; % [cm]
                L = Zaxis( pixi(end) ) *1e2 - Zaxis( pixi(1) ) *1e2; % [cm]

                % % Attenuation of the sample phantom
                % ATT_sam = attenuation_phantoms_Np(band, pars.SAM_num, pars.acs_sam); % [Np/cm]
                % cumul_att_data_uniform =  exp(-4 * Z_ROI_ini(i) * ATT_sam) .* ( (1 - exp(-ATT_sam * L))./(ATT_sam * L) ).^2; % [Np]  Uniform medium            
                % 
                % 
                % % Attenuation of the reference phantom
                % ATT_ref = attenuation_phantoms_Np(band, pars.REF_num, pars.acs_ref); % [Np/cm]
                % cumul_att_phan_uniform =  exp(-4 * Z_ROI_ini(i) * ATT_ref) .* ( (1 - exp(-ATT_ref * L))./(ATT_ref * L) ).^2; % [Np]  Uniform medium    

                %% CUMULATIVE OK
                X_ACS_ind = find(1e2*DATA.x >= X_ROI(j));
                Z_ACS_ind = find(1e2*DATA.z >= Z_ROI(i));

                ACS_sam = ones(size(rfdata_sam)).*DATA.acs/8.686;
           
                cumul_att_data =  exp(-4 * (delta_z*1e2)*sum(ACS_sam(1:Z_ACS_ind(1), X_ACS_ind(1)) ).* band.^m_sam); % [Np/cm]

                ACS_ref = REF.acs/8.686;
                

                cumul_att_phan = exp(-4 * Z_ROI(i) * ACS_ref * band.^m_ref);
                
                %% Backscatter coefficient estimation
                if  ~isfield(REF, 'BSC_ref')
                    BSC_data_w_uniform_map = squeeze(squeeze(SR(i,j,:)))'.* cumul_att_phan ./ cumul_att_data;
                else
                    BSC_data_w_uniform_map = squeeze(squeeze(SR(i,j,:)))'.* REF.BSC_ref(band) .* cumul_att_phan ./ cumul_att_data;
                end

                BSC_w_uniform_map(i,j,:)   =  BSC_data_w_uniform_map;
            
                % Calcul of BSC parameters: Slope, Intercept,
                % Slope evaluation and intergrated BSC
                [Sparam, Cparam, Midband, ~, iBSC] = BSC_parameters(band', BSC_data_w_uniform_map );
                Slope_BSC_Uni(i,j) = Sparam;
                Inter_BSC_Uni(i,j) = Cparam;
                Midba_BSC_Uni(i,j) = Midband;
                Integ_BSC_Uni(i,j) = iBSC;
                clear Sparam Cparam Midband iBSC;
            

            end
        end
    end
    
time_t = toc;
fprintf('Exec Time: %.4f \n', time_t); 

    Integ_BSC_UniPrev = Integ_BSC_Uni;

    medianBSC   = median(Integ_BSC_Uni(:), 'omitnan');
    stdBSC      = std(Integ_BSC_Uni(:), 'omitnan');

    Integ_BSC_Uni(Integ_BSC_UniPrev < medianBSC-2*stdBSC)= NaN;
    Integ_BSC_Uni(Integ_BSC_UniPrev > medianBSC+2*stdBSC)= NaN;
    clear iBSC_in_line CV;

    countBSC_Uni = 0;
        for i = 1 : size(BSC_w_uniform_map,1)
            for j = 1 : size(BSC_w_uniform_map,2)
                if  isnan(Integ_BSC_Uni(i,j)) == 0
                    countBSC_Uni = countBSC_Uni + 1;
                    BSC_in_line_Uni(:, countBSC_Uni)    = squeeze(squeeze(BSC_w_uniform_map(i,j,:)));
                    pBSC_Uni(countBSC_Uni)              = Slope_BSC_Uni(i,j);
                    cBSC_Uni(countBSC_Uni)              = Inter_BSC_Uni(i,j);
                    mBSC_Uni(countBSC_Uni)              = Midba_BSC_Uni(i,j);
                end
                
            end
        end

        Integ_BSC_Uni(Integ_BSC_Uni == 0) = NaN;
        CData_Uni = log10(Integ_BSC_Uni);
        
        Integ_BSC_in_line = Integ_BSC_Uni(:);
        Integ_BSC_in_line(isnan(Integ_BSC_in_line) == 1 ) = [];
        av_iBSC_Uni     = mean(Integ_BSC_in_line);
        md_iBSC_Uni     = median(Integ_BSC_in_line);
        sd_iBSC_Uni     = std(Integ_BSC_in_line);
        snr_iBSC_Uni    = av_iBSC_Uni/sd_iBSC_Uni;
        clear Integ_BSC_in_line;

        
        av_BSC_lesion_Uni = mean(BSC_in_line_Uni,2);
        md_BSC_lesion_Uni = median(BSC_in_line_Uni,2);
        if size(BSC_in_line_Uni,2) > 1
            sd_BSC_lesion_Uni  = std(BSC_in_line_Uni')';
        else
            sd_BSC_lesion_Uni = nan(length(md_BSC_lesion_Uni),1);
        end
        
        % figure(1002)
        % figure, 
        % semilogy(band, md_BSC_lesion_Uni, 'linewidth', 3); hold off;
        % set(gca, 'fontsize', 18);
        % xlabel('\bfFrequency [MHz]');
        % ylabel('\bfBSC [cm^{-1}\cdot sr^{-1}]');
        % title('\bfBSC')
        % grid on;
        % %axis([bw(1) bw(2) 1e-3 1e0]);
        % xlim(bw);

        BSC.X_ROI                   = X_ROI;
        BSC.Z_ROI                   = Z_ROI;
        
        %BSC.ACS_map_uniform   = ACS_map_uniform;
        
        BSC.BSC_w_uniform_map       = BSC_w_uniform_map ;
        BSC.countBSC_Uni            = countBSC_Uni;
        BSC.iBSC_Uni                = [av_iBSC_Uni md_iBSC_Uni sd_iBSC_Uni];
        BSC.pBSC_Uni                = [mean(pBSC_Uni) median(pBSC_Uni) std(pBSC_Uni)];
        BSC.mBSC_Uni                = [mean(mBSC_Uni) median(mBSC_Uni) std(mBSC_Uni)];
        BSC.cBSC_Uni                = [mean(cBSC_Uni) median(cBSC_Uni) std(cBSC_Uni)];
        BSC.BSCcurve_Uni            = [av_BSC_lesion_Uni md_BSC_lesion_Uni sd_BSC_lesion_Uni];
        
        BSC.band                    = band(:);
        BSC.nb_lambda               = nb_lambda_axial;
        BSC.overlap                 = overlap_axial;
    
return

%%
function Neff = effective_lines(data_block)

%data_block = randn(size(blockA));
%x = data_block(:,end);
%data_block = x*ones(1,N);


[~, N] =size(data_block);
RHO = corrcoef(data_block);
%v = 0:N-1;
%factor =toeplitz([v(1) fliplr(v(2:end))], v);

%val = factor.*(RHO.^2);
%val = sum(val);

rho = diagSumCalc(RHO,1);
rho = rho(1:N-1)./(1:N-1);

val = (rho.^2)*(1:N-1)';
Neff = N./( 1 + 2/N*( val ) );
%[mean(Neff) std(Neff) median(Neff) mode(Neff)]

% for ii = 1:d2-1,
%    y = data_block(:,ii);
%    bb = corrcoef(x,y);
%    rho(ii) = bb(1,2);
%    factor(ii) = ii;
%    %rho(ii) =  sum( xcorr(x,y') )/( norm(x)*norm(y) );
%
% end
%
% N = d2
% Neff = N/( 1 + 2/N*( factor*(rho').^2 ) )

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