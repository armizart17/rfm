function [BSC] = calculateBSC_RPM(DATA, REF, pars, Bmode)
% function [BSC] = caculateBSC_RPM(DATA, REF, pars, Bmode)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS: 
%         - DATA : Sample data
%                 - DATA.Bmode
%                 - DATA.RF
%                 - DATA.fc
%                 - DATA.fs
%                 - DATA.x
%                 - pars.z
%         - REF  : Reference phantom
%                 - REF.Bmode
%                 - REF.RF
%                 - REF.fc
%                 - REF.fs
%                 - REF.x
%                 - REF.z
%         - pars : 
%                 - pars.z_ini; ROI ini axial (depth)
%                 - pars.z_end; ROI end axial (depth)
%                 - pars.x_ini; ROI ini lateral
%                 - pars.x_end; ROI ini lateral
%                 - pars.BW
%                 - pars.nb_lambda_axial
%                 - pars.overlap_axial;
%                 - pars.overlap_lateral;
%                 - pars.P (number of points FFT)
%                 - pars.REF_num
%                 - pars.BSC_ref
%                 - pars.SAM_num
%         - Bmode : 
% OUTPUTS: 
%         - BSC: Object containing data BSC
% AUTHORs: Edmundo Arom Miranda & Andres Coila, based on LIM repository
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning('off');
fact_transparency = 0.5;

    % The bandwidth:
    BW = pars.BW;
    lambdac = (REF.c/1e3)/mean(BW); % The wavelength at the transducer center frequency
   
%     latlines = 256; % The number of lines in each Bmode
%     transducer_footprint = 38; % The transducer footprint in mm
%     delta_lat = transducer_footprint/latlines; % The distance in between adjacent lines
    
    nb_lambda_axial = pars.nb_lambda_axial; % [wavelength]
   
    overlap_axial = pars.overlap_axial;   
    overlap_lateral = pars.overlap_lateral;
    
                   
    % Begin the estimation process frame per frame
    nb_frame = size(DATA.RF, 3);
    for frame  = 1 : nb_frame
        tic;
        %Bande frequentiel de calcul
        P = pars.P; %for number of points 
        fct = 1e6; %for MHz
        Xfreq = (0 : (P/2-1))*DATA.fs/P/fct;    % [MHz]                           % [MHz]
        lambda = REF.c/(mean([BW(1) BW(2)])*1e6);% [wavelength]
        fprintf('\nCalculation bandwidth:\nBW = %g - %g MHz, Fc = %g MHz\nwavelength (wl) = %g mm\n', BW(1), BW(2), mean([BW(1) BW(2)]), lambda*1e3);
        Fi = find(Xfreq >= BW(1));
        Ff = find(Xfreq >= BW(2));
        band_ind = Fi(1) : Ff(1);               % [sample]
        band = band_ind*REF.fs/fct/P;          % [MHz]
        band_k = 2*pi*band*1e6/REF.c;   % [rad]
        
        % Choice of window type
        window_type = 2; % 1. Hanning, 2. Tuckey-0.25, 3. Hamming, 4. Chebyshev
        
        % Transmisson coefficient for Saran layer
        T = transmission_saran(band);

        % Attenuation of th referece phantom

        REF_num = pars.REF_num;
        alpha_phan_freq = attenuation_phantoms_Np(band, REF_num, pars.acs_ref); % [Np/cm]
        
        % Backscatter
        BSC_ref = pars.BSC_ref;       

    
    Zinitial = pars.z_ini;
    Zfinal   = pars.z_end;

    Xinitial = pars.x_ini;
    Xfinal = pars.x_end;
     
    Jmin = find(DATA.x > Xinitial);  
    Jmax = find(DATA.x > Xfinal);
    
    Imin = find(DATA.z > Zinitial);  
    Imax = find(DATA.z > Zfinal);  

    Zaxis = DATA.z(Imin(1):Imax(1));
    Xaxis = DATA.x(Jmin(1):Jmax(1));
    
    
    DATA_ROI_seg = DATA.RF (Imin(1):Imax(1), Jmin(1):Jmax(1), 1);
    Bmode_ROI_seg = DATA.Bmode(Imin(1):Imax(1), Jmin(1):Jmax(1), 1);       
    PHAN_ROI_seg = REF.RF (Imin(1):Imax(1), Jmin(1):Jmax(1), 1);
    
    [I2,J2] = size(DATA_ROI_seg);

    % Axial window lenght
    axial_gate_length   = nb_lambda_axial*lambda;          % [m]
    Lo2 = axial_gate_length*1e2/2;                         % [cm]
    delta_z = (Zaxis(end)-Zaxis(1))/length(Zaxis);         % [m]
    nb_sample_ROI_axial = floor(axial_gate_length/delta_z);% [sample]
    
    % Axial overlap
    nb_sample_overlap_axial = fix(nb_sample_ROI_axial*overlap_axial);% [sample]
    nb_ROI_axial = fix((I2-nb_sample_overlap_axial)/(nb_sample_ROI_axial - nb_sample_overlap_axial));
    delta_interpixel_axial = (Zaxis(end)-Zaxis(1))/ nb_ROI_axial;
    fprintf('\nAxial gate length = %d wl = %g mm = %g samp.\nNumber of axial region = %g\n', ...
        nb_lambda_axial, axial_gate_length*1e3, nb_sample_ROI_axial, nb_ROI_axial);
    
    % Lateral windows length
%     nb_lambda_lateral = nb_line_ROI_lateral*delta_lat/lambdac;
    nb_lambda_lateral = nb_lambda_axial;
    lateral_gate_length  = nb_lambda_lateral*lambda;               % [m]
    indice_lateral = find((Xaxis-Xaxis(1)) >= lateral_gate_length);
    nb_line_ROI_lateral = indice_lateral(1); clear indice_lateral;

    % Recouvrement inter-fenetre lateral
    nb_sample_overlap_lateral = fix(nb_line_ROI_lateral*overlap_lateral);    % [sample]
    nb_ROI_lateral = fix((J2-nb_sample_overlap_lateral)/(nb_line_ROI_lateral - nb_sample_overlap_lateral));
    delta_interpixel_lateral = (Xaxis(end)-Xaxis(1))/ nb_ROI_lateral;
    fprintf('\nLateral gate length = %2.1f wl = %g mm = %g lines\nNumber of lateral region = %g\n',...
    nb_lambda_lateral, lateral_gate_length*1e3, nb_line_ROI_lateral, nb_ROI_lateral );
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X_ROI   = zeros(1, nb_ROI_lateral);
    Z_ROI   = zeros(1, nb_ROI_axial);
    Z_ROI_start   = zeros(1, nb_ROI_axial);
    
    BSC_w_uniform_map   = nan(nb_ROI_axial, nb_ROI_lateral, length(band));
    
    Slope_BSC_Uni = nan(nb_ROI_axial, nb_ROI_lateral);
    Inter_BSC_Uni = nan(nb_ROI_axial, nb_ROI_lateral);
    Midba_BSC_Uni = nan(nb_ROI_axial, nb_ROI_lateral);
    Integ_BSC_Uni = nan(nb_ROI_axial, nb_ROI_lateral);
    ESD_Expo_BSC_Uni = nan(nb_ROI_axial, nb_ROI_lateral);
    EAC_Expo_BSC_Uni = nan(nb_ROI_axial, nb_ROI_lateral);
    ESD_Gauss_BSC_Uni = nan(nb_ROI_axial, nb_ROI_lateral);
    EAC_Gauss_BSC_Uni = nan(nb_ROI_axial, nb_ROI_lateral);
    gof_Expo_BSC_Uni = nan(nb_ROI_axial, nb_ROI_lateral);
    gof_Gauss_BSC_Uni = nan(nb_ROI_axial, nb_ROI_lateral);
    
    countP = 0;
    

           
    % Main loop -> data block per data block
    % where averaged power spectrum are calculated and employed for
    % the parameter estimation
    for i = 1 : nb_ROI_axial
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
        
        for j = 1 : nb_ROI_lateral
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
                phan_block = remove_DC_in_pixel_ref(squeeze(mean(phan_block)), phan_block, mat_window);
                data_block = remove_DC_in_pixel    (mean(data_block), data_block, mat_window);
                
                % Averaged power spectra
                tf_phan = power_spectrum_averaged_in_pixel_ref(phan_block,P,0,0);
                tf_data = power_spectrum_averaged_in_pixel(data_block,P,0,0);
                
                aux_phan = 10*log10(tf_phan/(max(tf_phan)));
                aux_data = 10*log10(tf_data/(max(tf_data)));

                if (i==1 && j==1)
                
                    figure (111),
                    subplot 211,
                    plot(Xfreq, aux_data' ), hold on, grid minor;
                    plot(Xfreq, aux_phan' ); hold off;
                    xlabel('Frequency [MHz]'), ylabel('Power Spectrum INI [dB]')
                    legend('Sample', 'Ref');
                    axis([0 max(Xfreq) -60 0]); %xlim 0-20MHZ y ylim -60dB 0dB
                    title ('Initial Power Spectrum (min depth)')
                
                end
                
                if (i==nb_ROI_axial && j == nb_ROI_lateral)
                
                    figure (111),
                    subplot 212,
                    plot(Xfreq, aux_data' ), hold on, grid minor;
                    plot(Xfreq, aux_phan' ); hold off;
                    xlabel('Frequency [MHz]'), ylabel('Power Spectrum END [dB]')
                    legend('Sample', 'Ref');
                    axis([0 max(Xfreq) -60 0]); %xlim 0-20MHZ y ylim -60dB 0dB    
                    title ('Final Power Spectrum (max depth)')
                    
                end
%                 figure (111),
%                 plot(Xfreq, aux_data' ), hold on, grid minor;
%                 plot(Xfreq, aux_phan' ); hold off;
%                 xlabel('Frequency [MHz]'), ylabel('Power Spectrum [dB]')
%                 legend('Data', 'Sample');
%                 axis([0 max(Xfreq) -60 0]); %xlim 0-20MHZ y ylim -60dB 0dB
% %                 pause(0.007);

%                 figure,
%                 plot(band, 10*log10( squeeze( squeeze( SR(i,j,:))) )' ), title('SR')



                % Selection of the useful bandwidth
                
                %tf_phan = tf_phan(band_ind)./T(:); % If both phantoms have Saran wrap, T = 1
                tf_phan = tf_phan(band_ind);
                tf_data = tf_data(band_ind);
                
                % Spectral ratio for BSC estimation
                SR(i,j,:) = tf_data ./ tf_phan;
                
                % Compensation of the attenuation
                X_ACS_ind = find(1e2*DATA.x >= X_ROI(j));
                Z_ACS_ind = find(1e2*DATA.z >= Z_ROI(i));
                
                l1 = length( 1:Z_ACS_ind(1) );
                l2 = length( 1: X_ACS_ind(1) );

               
                %ACS_map_uniform = 0.7359*ones(l1, l2);
                
                
                %% Cumulative attenution for BSC compensation for uniform media
                % See Eq. (5) in 10.1016/j.ultras.2021.106376
                %Z_ROI(i) = Zaxis( round((pixi(end) + pixi(1))/2 ) )*1e2; % [cm]
                Z_ROI_ini(i) = Zaxis( pixi(1) ) *1e2; % [cm]
                L = Zaxis( pixi(end) ) *1e2 - Zaxis( pixi(1) ) *1e2; % [cm]
%                 if j == 1 
%                     Z_ROI(i)
%                     Z_ROI
%                     keyboard
%                 end
                % Attenuation of the sample phantom
                ATT_sam = attenuation_phantoms_Np(band, pars.SAM_num, pars.acs_sam); % [Np/cm]
                cumul_att_data_uniform =  exp(-4 * Z_ROI_ini(i) * ATT_sam) .* ( (1 - exp(-ATT_sam * L))./(ATT_sam * L) ).^2; % [Np]  Uniform medium            

                % Attenuation of the reference phantom
                ATT_ref = attenuation_phantoms_Np(band, pars.REF_num, pars.acs_ref); % [Np/cm]
                cumul_att_phan_uniform =  exp(-4 * Z_ROI_ini(i) * ATT_ref) .* ( (1 - exp(-ATT_ref * L))./(ATT_ref * L) ).^2; % [Np]  Uniform medium    
              
                %%
                % Backscatter coefficient estimation
%                 total_att_uniform  = cumul_att_data_uniform - alpha_phan_freq*Z_ROI(i);
                %total_att_uniform  =  + alpha_phan_freq*Z_ROI(i);
                BSC_data_w_uniform_map = squeeze(squeeze(SR(i,j,:)))'.* BSC_ref .* cumul_att_phan_uniform ./ cumul_att_data_uniform;
                BSC_w_uniform_map(i,j,:)        =  BSC_data_w_uniform_map;
                
                % Calcul of BSC parameters: Slope, Intercept,
                % Slope evaluation and intergrated BSC
                [Sparam, Cparam, Midband, ~, iBSC] = BSC_parameters(band', BSC_data_w_uniform_map );
                Slope_BSC_Uni(i,j) = Sparam;
                Inter_BSC_Uni(i,j) = Cparam;
                Midba_BSC_Uni(i,j) = Midband;
                Integ_BSC_Uni(i,j) = iBSC;
                clear Sparam Cparam Midband iBSC;
                
                params.gauss = 1; params.fluid = 0; params.expo = 1;
                params.plot = 0;
                [a_gauss, a_expo, ~, EAC_gauss, EAC_expo, ~, ~, ~, ~, gof] = QUS_estim(BSC_data_w_uniform_map, band_k, params);
                ESD_Expo_BSC_Uni(i,j) = 2*a_expo*1e6;
                EAC_Expo_BSC_Uni(i,j) = EAC_expo;
                gof_Expo_BSC_Uni(i,j) = gof.expo;
                
                ESD_Gauss_BSC_Uni(i,j) = 2*a_gauss*1e6;
                EAC_Gauss_BSC_Uni(i,j) = EAC_gauss;
                gof_Gauss_BSC_Uni(i,j) = gof.gauss;
            end
        end
        
        time_t = toc;
    end
        fprintf('Time: %.4f \n', time_t); 
        
        %% For the uniform map
        % Removing outliers:
        
        Integ_BSC_UniPrev = Integ_BSC_Uni;
%         medianBSC = nanmedian(Integ_BSC_Uni(:));
        medianBSC = median(Integ_BSC_Uni(:), 'omitnan');
%         stdBSC = nanstd(Integ_BSC_Uni(:));
        stdBSC = std(Integ_BSC_Uni(:), 'omitnan');
        Integ_BSC_Uni(Integ_BSC_UniPrev<medianBSC-2*stdBSC)= NaN;
        Integ_BSC_Uni(Integ_BSC_UniPrev>medianBSC+2*stdBSC)= NaN;
        clear iBSC_in_line CV;
        
        ESD_Expo_BSC_UniPrev = ESD_Expo_BSC_Uni;
        EAC_Expo_BSC_UniPrev = EAC_Expo_BSC_Uni;
        gof_Expo_BSC_UniPrev = gof_Expo_BSC_Uni;
        ESD_Gauss_BSC_UniPrev = ESD_Gauss_BSC_Uni;
        EAC_Gauss_BSC_UniPrev = EAC_Gauss_BSC_Uni;
        gof_Gauss_BSC_UniPrev = gof_Gauss_BSC_Uni;
        
%         medianESD = nanmedian(ESD_Expo_BSC_Uni(:));
        medianESD = median(ESD_Expo_BSC_Uni(:),'omitnan');
%         stdESD = nanstd(ESD_Expo_BSC_Uni(:));
        stdESD = std(ESD_Expo_BSC_Uni(:),'omitnan');

        ESD_Expo_BSC_Uni(ESD_Expo_BSC_UniPrev<medianESD-2*stdESD)=NaN;
        ESD_Expo_BSC_Uni(ESD_Expo_BSC_UniPrev>medianESD+2*stdESD)=NaN;
        EAC_Expo_BSC_Uni(ESD_Expo_BSC_UniPrev<medianESD-2*stdESD)=NaN;
        EAC_Expo_BSC_Uni(ESD_Expo_BSC_UniPrev>medianESD+2*stdESD)=NaN;
        gof_Expo_BSC_Uni(ESD_Expo_BSC_UniPrev<medianESD-2*stdESD)=NaN;
        gof_Expo_BSC_Uni(ESD_Expo_BSC_UniPrev>medianESD+2*stdESD)=NaN;
        
%         medianESD = nanmedian(ESD_Gauss_BSC_Uni(:));
        medianESD = median(ESD_Gauss_BSC_Uni(:),'omitnan');
%         stdESD = nanstd(ESD_Gauss_BSC_Uni(:));
        stdESD = std(ESD_Gauss_BSC_Uni(:),'omitnan');

        ESD_Gauss_BSC_Uni(ESD_Gauss_BSC_UniPrev<medianESD-2*stdESD)=NaN;
        ESD_Gauss_BSC_Uni(ESD_Gauss_BSC_UniPrev>medianESD+2*stdESD)=NaN;
        EAC_Gauss_BSC_Uni(ESD_Gauss_BSC_UniPrev<medianESD-2*stdESD)=NaN;
        EAC_Gauss_BSC_Uni(ESD_Gauss_BSC_UniPrev>medianESD+2*stdESD)=NaN;
        gof_Gauss_BSC_Uni(ESD_Gauss_BSC_UniPrev<medianESD-2*stdESD)=NaN;
        gof_Gauss_BSC_Uni(ESD_Gauss_BSC_UniPrev>medianESD+2*stdESD)=NaN;
        
        countBSC_Uni = 0; clear i j;
        for i = 1 : size(BSC_w_uniform_map,1)
            for j = 1 : size(BSC_w_uniform_map,2)
                if      isnan(Integ_BSC_Uni(i,j)) == 0
                    countBSC_Uni = countBSC_Uni + 1;
                    BSC_in_line_Uni(:, countBSC_Uni)    = squeeze(squeeze(BSC_w_uniform_map(i,j,:)));
                    pBSC_Uni(countBSC_Uni)         = Slope_BSC_Uni(i,j);
                    cBSC_Uni(countBSC_Uni)         = Inter_BSC_Uni(i,j);
                    mBSC_Uni(countBSC_Uni)         = Midba_BSC_Uni(i,j);
                end
                
            end
        end
        
        if exist('pBSC_Uni', 'var') == 0
            pBSC_Est = nan;
            cBSC_Est = nan;
            mBSC_Est = nan;
        end
        
        Integ_BSC_Uni(Integ_BSC_Uni == 0) = NaN;
        CData_Uni = log10(Integ_BSC_Uni);
        
        Integ_BSC_in_line = Integ_BSC_Uni(:);
        Integ_BSC_in_line(isnan(Integ_BSC_in_line) == 1 ) = [];
        av_iBSC_Uni = mean(Integ_BSC_in_line);
        md_iBSC_Uni = median(Integ_BSC_in_line);
        sd_iBSC_Uni = std(Integ_BSC_in_line);
        snr_iBSC_Uni = av_iBSC_Uni/sd_iBSC_Uni;
        clear Integ_BSC_in_line;
        
        figure(1000)
        transparency_ibsc = ones(size(CData_Uni,1), size(CData_Uni,2));
        transparency_ibsc(isnan(CData_Uni) == 0) = 1;
        transparency_ibsc(isnan(CData_Uni) == 1) = 0;
        subimage(1e2*DATA.x,1e2*DATA.z, 256*mat2gray(Bmode, [-60 0]), gray); hold on;
%         h = subimage(X_ROI, Z_ROI, 256*mat2gray(CData_Uni, [-4 -1.8]), jet); hold off;
        h = subimage(X_ROI, Z_ROI, 256*mat2gray(CData_Uni, [-6 -4]), jet);
%          h = subimage(X_ROI, Z_ROI, 256*mat2gray(CData_Uni), jet); hold off;
        ibsc_trans = transparency_ibsc*fact_transparency;
        set( h, 'AlphaData', ibsc_trans) ;
        xlabel('\bfLateral distance [cm]','fontsize',12);
        ylabel('\bfAxial distance [cm]','fontsize',12)
        title(['\bfiBSC = ' num2str(av_iBSC_Uni, '%.2s' ) ' \pm' num2str(sd_iBSC_Uni, '%.2s' ) ' sr^{-1}cm^{-1}']);
        set(gca,'fontsize',12);
        colormap jet
        h2 = colorbar;
        clim([-0.0006 0.004])
        ylabel(h2,'sr^{-1}cm^{-1}','FontSize', 14);
%         saveas(gcf,[data_folder  'frame#' int2str(frame) '_BSC_image.fig']);
%         saveas(gcf,[data_folder  'frame#' int2str(frame) '_BSC_image.png'], 'png');
        
        figure(2002)
        transparency_esd = ones(size(ESD_Expo_BSC_Uni,1), size(ESD_Expo_BSC_Uni,2));
        transparency_esd(isnan(ESD_Expo_BSC_Uni) == 0) = 1;
        transparency_esd(isnan(ESD_Expo_BSC_Uni) == 1) = 0;
        subimage(1e2*DATA.x,1e2*DATA.z, 64*mat2gray(Bmode), gray(64)); hold on;
        h = subimage(X_ROI, Z_ROI, 64*mat2gray(ESD_Expo_BSC_Uni, [1 200]), jet(64)); hold off
        esd_trans = transparency_esd*fact_transparency;
        set( h, 'AlphaData', esd_trans) ;
        xlabel('\bfLateral distance [cm]','fontsize',12);
        ylabel('\bfAxial distance [cm]','fontsize',12)
%         title(['\bfESD = ' num2str(nanmean(ESD_Expo_BSC_Uni(:)), 4 ) ' \pm' num2str(nanstd(ESD_Expo_BSC_Uni(:)), 4 ) ' \mu m']);
        title(['\bfESD = ' num2str(mean(ESD_Expo_BSC_Uni(:), 'omitnan'), 4 ) ' \pm' num2str(std(ESD_Expo_BSC_Uni(:), 'omitnan'), 4 ) ' \mu m']);
        set(gca,'fontsize',12);
%         saveas(gcf,[data_folder  'frame#' int2str(frame) '_ESD_image.fig']);
%         saveas(gcf,[data_folder  'frame#' int2str(frame) '_ESD_image.png'], 'png');
        
        figure
        subplot(121)
        imagesc(X_ROI,Z_ROI,CData_Uni)
        xlabel('\bfLateral distance [cm]','fontsize',12);
        ylabel('\bfAxial distance [cm]','fontsize',12)
        title(['\bfiBSC = ' num2str(av_iBSC_Uni, '%.2s' ) ' \pm' num2str(sd_iBSC_Uni, '%.2s' ) ' sr^{-1}cm^{-1}']);
        colorbar
        subplot(122)
        imagesc(X_ROI,Z_ROI,(Integ_BSC_Uni))
        xlabel('\bfLateral distance [cm]','fontsize',12);
        ylabel('\bfAxial distance [cm]','fontsize',12)
        title(['\bfPrev-iBSC = ' num2str(av_iBSC_Uni, '%.2s' ) ' \pm' num2str(sd_iBSC_Uni, '%.2s' ) ' sr^{-1}cm^{-1}']);
        colorbar
        
        figure(1002)
        av_BSC_lesion_Uni = mean(BSC_in_line_Uni,2);
        md_BSC_lesion_Uni = median(BSC_in_line_Uni,2);
        if size(BSC_in_line_Uni,2) > 1
            sd_BSC_lesion_Uni  = std(BSC_in_line_Uni')';
        else
            sd_BSC_lesion_Uni = nan(length(md_BSC_lesion_Uni),1);
        end
        
        
        semilogy(band, md_BSC_lesion_Uni, 'linewidth', 3); hold off;
        set(gca, 'fontsize', 18);
        xlabel('\bfFrequency [MHz]');
        ylabel('\bf<BSC> \sigma [1/(sr.cm)]');
        grid on;
        %axis([BW(1) BW(2) 1e-3 1e0]);
        axis([BW(1) BW(2) 1e-6 1e-1]);
        
        figure, 
        plot(Integ_BSC_Uni(:,10))

%         saveas(gcf,[data_folder  'frame#' int2str(frame) '_BSC_image2.fig']);
        
        
        
        % Best-fit: least square minimization - UNI
        [a_gaus, a_expo, a_fluid, EAC_gaus, EAC_expo, EAC_fluid, fit_gaus, fit_expo, fit_fluid, gof] = QUS_estim(md_BSC_lesion_Uni', band_k);
        ESD_Uni = 2*[a_gaus a_expo a_fluid]; clear a_gaus a_expo a_fluid;
        EAC_Uni = [EAC_gaus EAC_expo EAC_fluid]; clear EAC_gaus EAC_expo EAC_fluid;
        FIT_Uni = [fit_gaus; fit_expo; fit_fluid]; clear fit_gaus fit_expo fit_fluid;
        
        figure(1001)
        semilogy(band, md_BSC_lesion_Uni, 'linewidth', 3); hold on;
        semilogy(band, FIT_Uni, ':','linewidth', 2); hold off;
        xlabel('\bfFrequency [MHz]','fontsize', 12);
        ylabel('\bf<BSC> \sigma [1/(sr.cm)]','fontsize', 12);
        title('Att. Compensation : uniform value')
        legend('Estimated BSC', 'Gaussian form factor', 'Exponential form factor', 'Sp. Fluid form factor')
%         saveas(gcf,[data_folder  'frame#' int2str(frame) '_BSC_form factor_lesion.fig']);
        
        clear BSC;
        

        BSC.X_ROI = X_ROI;
        BSC.Z_ROI = Z_ROI;
        
        %BSC.ACS_map_uniform   = ACS_map_uniform;
        
        BSC.BSC_w_uniform_map = BSC_w_uniform_map ;
        BSC.countBSC_Uni = countBSC_Uni;
        BSC.iBSC_Uni = [av_iBSC_Uni md_iBSC_Uni sd_iBSC_Uni];
        BSC.pBSC_Uni = [mean(pBSC_Uni) median(pBSC_Uni) std(pBSC_Uni)];
        BSC.mBSC_Uni = [mean(mBSC_Uni) median(mBSC_Uni) std(mBSC_Uni)];
        BSC.cBSC_Uni = [mean(cBSC_Uni) median(cBSC_Uni) std(cBSC_Uni)];
        BSC.BSCcurve_Uni = [av_BSC_lesion_Uni md_BSC_lesion_Uni sd_BSC_lesion_Uni];
        BSC.ESD_Uni = ESD_Uni;
        BSC.FIT_Uni = FIT_Uni;
        BSC.EAC_Uni = EAC_Uni;
        BSC.GoF = gof;
        
        BSC.band = band;
        BSC.nb_lambda = nb_lambda_axial;
        BSC.overlap = overlap_axial;
%         BSC.State =  char(patient_list(num_patient,4));
%         BSC.Preset = char(patient_list(num_patient,3));
        
        BSC.ESD_Expo_Map = ESD_Expo_BSC_Uni;
        BSC.EAC_Expo_Map = EAC_Expo_BSC_Uni;
        BSC.gof_Expo_Map = gof_Expo_BSC_Uni;
        BSC.ESD_Expo_Map_Outliers = ESD_Expo_BSC_UniPrev;
        BSC.EAC_Expo_Map_Outliers = EAC_Expo_BSC_UniPrev;
        BSC.gof_Expo_Map_Outliers = gof_Expo_BSC_UniPrev;
        
        BSC.Integ_BSC_Map = Integ_BSC_Uni;
        BSC.ESD_Gauss_Map = ESD_Gauss_BSC_Uni;
        BSC.EAC_Gauss_Map = EAC_Gauss_BSC_Uni;
        BSC.gof_Gauss_Map = gof_Gauss_BSC_Uni;
        BSC.Integ_BSC_Map_Outliers = Integ_BSC_UniPrev;
        BSC.ESD_Gauss_Map_Outliers = ESD_Gauss_BSC_UniPrev;
        BSC.EAC_Gauss_Map_Outliers = EAC_Gauss_BSC_UniPrev;
        BSC.gof_Gauss_Map_Outliers = gof_Gauss_BSC_UniPrev;
             
%         save(out_BSC_estimation, 'BSC');
        
        BSC.CData_Uni = CData_Uni;
        BSC.Integ_BSC_Uni = Integ_BSC_Uni;
    end
    
end

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

end
