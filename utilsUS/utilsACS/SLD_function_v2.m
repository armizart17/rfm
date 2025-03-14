function [SLD, DATA, REF] = SLD_function_v2(DATA, REF, pars)
% function [SLD] = SLD_function_v2(DATA, REF, pars)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION : SLD technique withouth SNR evaluation
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
%                 - pars.cs % speed of sound SoS in [m/s]
%                 - pars.z_ini; % ROI ini axial (depth) % in [m]
%                 - pars.z_end; % ROI end axial (depth) % in [m]
%                 - pars.x_ini; % ROI ini lateral % in [m]
%                 - pars.x_end; % ROI ini lateral % in [m]
%                 - pars.bw BANDWITH i.e [ 4 10 ] % in [MHz]
%                 - pars.nb_lambda_axial % datablock in wavelengths
%                 - pars.overlap_axial; 
%                 - pars.overlap_lateral;
%                 - pars.P % (number of points FFT)
%                 - pars.REF_num % CHOICE OF REFERENCE PHANTOM CHECK attenuation_phantoms_Np.m 
%                 - pars.SLOPE % In case of linear attenuation (simulation)
% OUTPUTS: 
%         - SLD : 
%                 - SLD_term: Spectral Log Difference Term of nºchannels of frequency 3D array (axial, lateral, band)
%                 - band: Frequency band vector [MHz]
%                 - A : A matrix por SLD as 2D-Inverse problem [-4 zp-zd f, speye(m*n)  | 1*f speye(m*n)]
%                 - x_ori : DATA.x 
%                 - z_ori : DATA.z 
%                 - x : X_ROI
%                 - z : Z_ROI
%                 - SLD.Bmode = DATA.Bmode;
% AUTHORs: Edmundo Arom Miranda & Andres Coila, based on LIM repository 
%          (Julien Rouyer Code TVSLD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % DEFINE DATA & SAMPLE      
    
    % ROI
%     Zinitial = 1.85e-2; % in [m]
%     Zfinal   = 4.25e-2; % in [m]
% 
%     Xinitial = 0.25e-2; % in [m];
%     Xfinal = 3.5e-2; % in [m];
    

    %%%%%%%%%%%%% LATER (UNPACKAGE pars) %%%%%%%%%%%%%
    warning('off');
    close(figure(17)), close(figure(171)), close(figure(172))
    
    z_ini           = pars.z_roi(1);
    z_end           = pars.z_roi(2);
    x_ini           = pars.x_roi(1);
    x_end           = pars.x_roi(2);
    bw              = pars.bw; % [MHz]
    overlap         = pars.overlap;
    blocksize_wv    = pars.blocksize;
    
    Zinitial = z_ini; % in [m]
    Zfinal   = z_end; % in [m]
    Xinitial = x_ini;  % in [m];
    Xfinal   = x_end; % in [m];

    % P = pars.P; % number of points
    P = 1024;
    
    nb_lambda_axial   = blocksize_wv;
    nb_lambda_lateral = blocksize_wv;
    overlap_axial     = overlap;
    overlap_lateral   = overlap;

    % Default window
    window_type     = 2; %  (1) Hanning, (2) Tukey, (3) Hamming, (4) Tchebychev
    if isfield(pars, 'window_type')
        window_type = pars.window_type;
    end

    cs = 1540;
    REF_num = pars.REF_num;
    
    if REF_num == 111 % CASE LINEAR DEPENDENCY
        SLOPE = pars.SLOPE; % [dB/cm-Hz]
    else
        SLOPE = NaN; % DISCARD
    end
    
    %%%%%%%%%%%%% LATER (UNPACKAGE pars) %%%%%%%%%%%%%

    %Bande frequentiel de calcul
    fct = 1e6;
    Xfreq = (0 : (P/2-1))*REF.fs/P/fct;    % [MHz]
   
    lambda = cs/(mean([bw(1) bw(2)])*1e6); % [wavelength]
    fprintf('\nCalculation bandwidth:\nBW = %g - %g MHz, Fc = %g MHz\nwavelength (wl) = %g mm\n', bw(1), bw(2), mean([bw(1) bw(2)]), lambda*1e3);
    Fi       = find(Xfreq >= bw(1), 1, 'first');
    Ff       = find(Xfreq <= bw(2), 1, 'last');
    band_ind = Fi:Ff;
    band     = Xfreq(band_ind);       % [MHz]
    
    % Choice of window type
%     window_type = 2; % 1. Hanning, 2. Tuckey-0.25, 3. Hamming, 4. Chebyshev
    
    % Transmisson coefficient for Saran layer
%     T = transmission_saran(band);
    T = 1;

    % Attenuation of th referece phantomedm
    alpha_phan_freq = attenuation_phantoms_Np(band, REF_num, SLOPE); % [Np/cm]
    

    [alpha4, ~, y_lin_ref4, ~] = fit_linear(band,alpha_phan_freq',1);  % [Np/cm/MHz]
    
    
    %%%%%%%%%%%%%%%% FIND %%%%%%%%%%%%%%%%
    Imin  = find(DATA.z > Zinitial);  
    Imax  = find(DATA.z > Zfinal);  

    Jmin = find(DATA.x > Xinitial);  
    Jmax = find(DATA.x > Xfinal);

    Zaxis = DATA.z(Imin(1):Imax(1));
    Xaxis = DATA.x(Jmin(1):Jmax(1));
    ar    = 0;
    
    DATA_ROI_seg = DATA.rf    (Imin(1):Imax(1), Jmin(1):Jmax(1), 1);
    % Bmode_ROI_seg = DATA.Bmode(Imin(1):Imax(1), Jmin(1):Jmax(1), 1);
    REF_ROI_seg = REF.rf (Imin(1):Imax(1), Jmin(1):Jmax(1), :);
    
    %%%%%%%%%%%%%%%% FIND %%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%% ROI CROP V2 %%%%%%%%%%%%
    ind_x = Xinitial <= DATA.x & DATA.x <= Xfinal;
    ind_z = Zinitial <= DATA.z & DATA.z <= Zfinal;

    Zaxis2 = DATA.z(ind_z);
    Xaxis2 = DATA.x(ind_x);

    DATA_ROI_seg2 = DATA.rf (ind_z, ind_x);
    REF_ROI_seg2 = REF.rf(ind_z, ind_x);

    %%%%%%%%%%%%%%%% ROI CROP V2 %%%%%%%%%%%%

    env_rfdata_sam = abs(hilbert( DATA.rf ));
    env_rfdata_sam_roi = env_rfdata_sam (Imin(1):Imax(1), Jmin(1):Jmax(1), 1);

    [I2,J2] = size(DATA_ROI_seg);
    
    % Preparation for the ROI (window) segmentation
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
%     nb_lambda_lateral = nb_lambda_axial; %round(nb_lambda_axial/2);
    lateral_gate_length  = nb_lambda_lateral*lambda;     % [m]
    indice_lateral = find((Xaxis-Xaxis(1)) >= lateral_gate_length);
    nb_line_ROI_lateral = indice_lateral(1); clear indice_lateral;
    
    % Recouvrement inter-fenetre lateral                             
    nb_sample_overlap_lateral = fix(nb_line_ROI_lateral*overlap_lateral);    % [sample]
    nb_ROI_lateral = fix((J2-nb_sample_overlap_lateral)/(nb_line_ROI_lateral - nb_sample_overlap_lateral));
    delta_interpixel_lateral = ( Xaxis(end)-Xaxis(1) )/ nb_ROI_lateral;
    fprintf('\nLateral gate length = %d wl = %g mm = %g lines\nNumber of lateral region = %g\n',...
        nb_lambda_lateral, lateral_gate_length*1e3, nb_line_ROI_lateral, nb_ROI_lateral );
    
    
    %%%%%%%%%%%%%% DISPLAY INFORMATION (COILA CODE) %%%%%%%%%%%%%% 
% 		fileID = fopen('exp.txt','w');
%         
%         fclose(fileID); 
%     
%         fprintf(['Frequency range: ',num2str(bw(1),'%3.1f'),' - ',num2str(bw(2),'%3.1f'),' MHz. cs: ',...
%                     num2str(cs,'%4.1f'),' m/s. Wavelength: ',num2str(lambda*1e3,'%2.2f'),' mm\n']);
% %         fprintf(['Blocksize. x: ',num2str(nx*dx*1e3,'%4.2f'),'mm, z: ',num2str(nz*dz*1e3,'%4.2f'),'mm, overlap: ',num2str(overlap*1e2,'%4.0f'),'%']);
%         fprintf(['Blocksize in wavelengths: ',num2str(new_blocksize,'%3.1f')]);
%         fprintf(['Blocksize in pixels. nf: ',num2str(p,'%i'),' nx: ',num2str(nx,'%i'),', nz: ',num2str(nz,'%i'),', nw: ',num2str(nw,'%i')]);
%         fprintf(['Region of interest. columns: ',num2str(ncol,'%i'),', rows: ',num2str(nrow,'%i')]);
    %%%%%%%%%%%%%% DISPLAY INFORMATION (COILA CODE) %%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%% RIGH X & Z AXIS SIZECOILA CODE %%%%%%%%%%%%%%   
                    
    % Number of repetitions of datablocks within a single
                % datablock
%                 RSLD.rpt = 1/(1-overlap_axial);   % r = rep = 1/(1-overlap)
%                 RSLD.rpt = round(RSLD.rpt);
%                 RSLD.overlap = 1 - (1/RSLD.rpt);
%                 
%                 % Part without lateral overlap
%                 RSLD.wx  = round((blocksize*wl)/(dx*RSLD.rpt));   
%                 % Number of lines laterally = repetitions * block without repetition
%                 nx  = rpt*wx;   
%                 new_blocksize = round(nx*dx/(wl));
%                 
%                 % RF data columns
%                 L2   = size(sam1,2);   
%                 % Blocksize colums
%                 ncol = floor((L2-(rpt-1)*wx)/wx);     
%                 sam1 = sam1(:,1:wx*(ncol+rpt-1));
%                 % Actual rf data columns
%                 L2 = size(sam1,2);   
%                 x  = x(1:L2);
%                 
%                 xi = 1;
%                 xf = L2;
%                 x0 = (xi:wx:xf+1-nx);
%                 x_ACS = x(1,x0+round(nx/2));
%                 nx_ACS = nx*dx;
%                 
%                 n  = length(x0);
%                 
%                 wz = floor(nx*dx/(dz*rpt));
%                 nz = rpt*wz;
%                 
%                 % winsize: Percentage of window (max 0.5)
%                 % nw: Samples of each window axially
%                 nw = 2*floor(winsize*nx*dx/(2*dz)) - 1 ;  
%                 L = (nz - nw)*dz*100;   % (cm)
%                                 
%                 NFFT = 2^(nextpow2(nw)+2);
%                 band = fs*linspace(0,1,NFFT)';   % [Hz] Band of frequencies
%                 rang = (floor(freq_L/fs*NFFT)+1:round(freq_H/fs*NFFT));   % useful frequency range
%                 f  = band(rang)*1e-6; % [MHz]
%                 L3 = length(rang);
%                 p = L3;
%                 
%                 L1   = size(sam1,1);   % RF data: rows
%                 nrow = floor((L1-(rpt-1)*wz)/wz);        % Blocksize rows
%                 sam1 = sam1(1:wz*(nrow+rpt-1),:);
%                 L1   = size(sam1,1);   % RF data: rows
%                 z    = z(1:L1);
%                 
%                 zi = 1;
%                 zf = L1;
%                 z0 = (zi:wz:zf+1-nz);
%                 m  = length(z0);
%                 
%                 z_ACS = z(z0+round(nz/2));
%                 nz_ACS = nz*dz;
%                 
%                 z0p = z0 + (nw-1)/2;
%                 z0d = z0 + (nz-1) - (nw-1)/2;

    %%%%%%%%%%%%%% RIGH X & Z AXIS SIZECOILA CODE %%%%%%%%%%%%%% 

    % Initialisation des matrices
    X_ROI   = zeros(1, nb_ROI_lateral);
    Z_ROI   = zeros(1, nb_ROI_axial);
    
    SLD_term = zeros(nb_ROI_axial, nb_ROI_lateral, length(band));
    SpectSNR = zeros(nb_ROI_axial, nb_ROI_lateral);
    TempoSNR = zeros(nb_ROI_axial, nb_ROI_lateral);
    TempoSNR_out = zeros(nb_ROI_axial, nb_ROI_lateral);

    og_SNR = zeros(nb_ROI_axial, nb_ROI_lateral);
    delta_SNR = zeros(nb_ROI_axial, nb_ROI_lateral);
    delta_SNR_noAbs = zeros(nb_ROI_axial, nb_ROI_lateral);
    
    for i = 1 : nb_ROI_axial
        
        pixi = (i-1)*(nb_sample_ROI_axial - nb_sample_overlap_axial) + (1:nb_sample_ROI_axial-1);
        
        pixi_pro = pixi(1:round(length(pixi)/2));
        pixi_dis = pixi(round((length(pixi)/2+1)):end);
        wind_pro = window_choice(length(pixi_pro),window_type);
        wind_dis = window_choice(length(pixi_dis),window_type);
        
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
            
            wind_pro = window_choice(length(pixi_pro),window_type);
            wind_dis = window_choice(length(pixi_dis),window_type);
        end
        
        
        for j = 1 : nb_ROI_lateral
            
            pixj = (j-1)*(nb_line_ROI_lateral - nb_sample_overlap_lateral) + (1:nb_line_ROI_lateral);
            mat_window = (ones(length(pixj),1)*wind')';
            
            mat_window_pro = (ones(length(pixj),1)*wind_pro')';
            mat_window_dis = (ones(length(pixj),1)*wind_dis')';
            
            
            if j == nb_ROI_lateral
                w = 0;
                while pixj(end) > J2
                    w = w + 1; clear pixj;
                    pixj = (j-1)*(nb_line_ROI_lateral - nb_sample_overlap_lateral) + (1:nb_line_ROI_lateral-w);
                    warning('\nFor j = %d, ROI final sample (%d) is larger than the ROW size (%d)\n', ...
                        j, pixj(end), J2);
                end
                mat_window = (ones(length(pixj),1)*wind')';
                
                mat_window_pro = (ones(length(pixj),1)*wind_pro')';
                mat_window_dis = (ones(length(pixj),1)*wind_dis')';
                
            end
            
            Z_ROI(i) = Zaxis( round((pixi(end) + pixi(1))/2 ) ); % [m]
            X_ROI(j) = Xaxis( round((pixj(end) + pixj(1))/2 ) ); % [m]
            
            Z_subROI_pro = Zaxis( round( (pixi_pro(end) + pixi_pro(1))/2 ) )*1e2;
            Z_subROI_dis = Zaxis( round( (pixi_dis(end) + pixi_dis(1))/2 ) )*1e2;
                        
            % Selection of data in data-block
            ref_block_pro = REF_ROI_seg(pixi_pro, pixj,:);
            ref_block_dis = REF_ROI_seg(pixi_dis, pixj,:);
            data_block_pro = DATA_ROI_seg(pixi_pro, pixj);
            data_block_dis = DATA_ROI_seg(pixi_dis, pixj);

            ref_block_pro = remove_DC_in_pixel_gen(squeeze(mean(ref_block_pro)), ref_block_pro, mat_window_pro);
            ref_block_dis = remove_DC_in_pixel_gen(squeeze(mean(ref_block_dis)), ref_block_dis, mat_window_dis);
            data_block_pro = remove_DC_in_pixel_gen    (mean(data_block_pro), data_block_pro, mat_window_pro);
            data_block_dis = remove_DC_in_pixel_gen    (mean(data_block_dis), data_block_dis, mat_window_dis);
                
            % Averaged power spectra
            tf_ref_pro = power_spectrum_averaged_in_pixel_gen(ref_block_pro,P,0,0);
            tf_ref_dis = power_spectrum_averaged_in_pixel_gen(ref_block_dis,P,0,0);
            tf_data_pro = power_spectrum_averaged_in_pixel_gen    (data_block_pro,P,0,0);
            tf_data_dis = power_spectrum_averaged_in_pixel_gen    (data_block_dis,P,0,0);
                       
            %% MULTIPLE PLOT V2.0 (COILA)

            if ( j == floor(2.5*nb_ROI_lateral/5) )
                if ( i == ceil(nb_ROI_axial/10))

                    
                    %%%% PROXIMAL %%%%
                    figure(171), 
                    plot( Xfreq,  10*log10(tf_data_pro/(max(tf_data_pro))), 'k', 'Linewidth', 2), hold on;
                    plot( Xfreq,  10*log10(tf_ref_pro/(max(tf_ref_pro))), 'k--', 'Linewidth', 1), grid on;
                    
                    xlabel('Frequency [MHz]'), ylabel('Power Spectrum [dB]')
                    axis([0 max(Xfreq) -60 0]); %xlim 0-20MHZ y ylim -60dB 0dB
                    title ('Power Spectrum Prox')

                    %%%% DISTAL %%%%
                    figure(172), 
                    plot( Xfreq,  10*log10(tf_data_dis/(max(tf_data_dis))), 'k', 'Linewidth', 2);hold on;
                    plot( Xfreq,  10*log10(tf_ref_dis/(max(tf_ref_dis))), 'k--', 'Linewidth', 1), grid on;
                    
                    xlabel('Frequency [MHz]'), ylabel('Power Spectrum [dB]')
                    axis([0 max(Xfreq) -60 0]); %xlim 0-20MHZ y ylim -60dB 0dB
                    title ('Power Spectrum Dist')
                    
                elseif ( i == round(nb_ROI_axial/2))

                    %%%% PROXIMAL %%%%
                    figure(171), 
                    plot( Xfreq,  10*log10(tf_data_pro/(max(tf_data_pro))), 'r', 'Linewidth', 2), hold on;
                    plot( Xfreq,  10*log10(tf_ref_pro/(max(tf_ref_pro))), 'r--', 'Linewidth', 1), grid on;
                                        

                    %%%% DISTAL %%%%
                    figure(172), 
                    plot( Xfreq,  10*log10(tf_data_dis/(max(tf_data_dis))), 'r', 'Linewidth', 2), hold on;
                    plot( Xfreq,  10*log10(tf_ref_dis/(max(tf_ref_dis))), 'r--', 'Linewidth', 1), grid on;
                    
                    xlabel('Frequency [MHz]'), ylabel('Power Spectrum [dB]')
                    axis([0 max(Xfreq) -60 0]); %xlim 0-20MHZ y ylim -60dB 0dB
                    

                elseif ( i == floor(9*nb_ROI_axial/10))

                    %%%% PROXIMAL %%%%
                    figure(171), 
                    plot( Xfreq,  10*log10(tf_data_pro/(max(tf_data_pro))), 'b', 'Linewidth', 2), hold on;
                    plot( Xfreq,  10*log10(tf_ref_pro/(max(tf_ref_pro))), 'b--', 'Linewidth', 1), grid on;
                    
                    % HORIZONTAL LINES REFERENCE
                    plot( Xfreq, -20*ones(size(Xfreq)),'k--');
                    plot( Xfreq, -15*ones(size(Xfreq)),'k--');
                    plot( Xfreq, -10*ones(size(Xfreq)),'k--');
                    xlabel('Frequency [MHz]'), ylabel('Power Spectrum [dB]')
                    axis([0 max(Xfreq) -60 0]); %xlim 0-20MHZ y ylim -60dB 0dB
                    

                    %%%% DISTAL %%%%
                    figure(172),  
                    plot( Xfreq,  10*log10(tf_data_dis/(max(tf_data_dis))), 'b', 'Linewidth', 2), hold on;
                    plot( Xfreq,  10*log10(tf_ref_dis/(max(tf_ref_dis))), 'b--', 'Linewidth', 1), grid on;
                    
                    % HORIZONTAL LINES REFERENCE
                    plot( Xfreq, -20*ones(size(Xfreq)),'k--');
                    plot( Xfreq, -15*ones(size(Xfreq)),'k--');
                    plot( Xfreq, -10*ones(size(Xfreq)),'k--');
                    xlabel('Frequency [MHz]'), ylabel('Power Spectrum [dB]')
                    axis([0 max(Xfreq) -60 0]); %xlim 0-20MHZ y ylim -60dB 0dB
                    
                end
            end

            % MULTIPLE PLOT V1.0 (EAMZ)
            if (i==ceil(nb_ROI_axial/10) && j==floor(2.5*nb_ROI_lateral/5))
                
                    figure (17),
                    subplot 211,
                    plot( Xfreq,  10*log10(tf_ref_pro/(max(tf_ref_pro))), 'k', 'Linewidth', 1), hold on, grid on;
%                     plot( Xfreq,  10*log10(tf_ref_dis/(max(tf_ref_dis))), ':k', 'Linewidth', 1);hold on;
        
                    plot( Xfreq,  10*log10(tf_data_pro/(max(tf_data_pro))), 'r', 'Linewidth', 2);
%                     plot( Xfreq,  10*log10(tf_data_dis/(max(tf_data_dis))), ':r', 'Linewidth', 2);hold on;

                    % HORIZONTAL LINES REFERENCE
                    plot( Xfreq, -20*ones(size(Xfreq)),'k--');
                    plot( Xfreq, -15*ones(size(Xfreq)),'k--');
                    plot( Xfreq, -10*ones(size(Xfreq)),'k--'), hold off;

                    xlabel('Frequency [MHz]'), ylabel('Power Spectrum INI [dB]')
%                     legend('REF Prox','REF Dist', 'DATA Prox', 'DATA Dist', 'Location','NorthEast');
                    legend('REF Prox','DATA Prox', 'Location','NorthEast');
                    axis([0 max(Xfreq) -60 0]); %xlim 0-20MHZ y ylim -60dB 0dB
                    title ('Initial Power Spectrum (min depth)')
                
                end
                
             if (i==floor(9*nb_ROI_axial/10) && j == floor(2.5*nb_ROI_lateral/5))
                
                    figure (17),
                    subplot 212,
                    plot( Xfreq,  10*log10(tf_ref_pro/(max(tf_ref_pro))), 'k', 'Linewidth', 1), hold on, grid on;
%                     plot( Xfreq,  10*log10(tf_ref_dis/(max(tf_ref_dis))), ':k', 'Linewidth', 1);hold on;
        
                    plot( Xfreq,  10*log10(tf_data_pro/(max(tf_data_pro))), 'r', 'Linewidth', 2);
%                     plot( Xfreq,  10*log10(tf_data_dis/(max(tf_data_dis))), ':r', 'Linewidth', 2);hold on;

                    % HORIZONTAL LINES REFERENCE
                    plot( Xfreq, -20*ones(size(Xfreq)),'k--');
                    plot( Xfreq, -15*ones(size(Xfreq)),'k--');
                    plot( Xfreq, -10*ones(size(Xfreq)),'k--'), hold off;
                    
                    xlabel('Frequency [MHz]'), ylabel('Power Spectrum END [dB]')
%                     legend('REF Prox','REF Dist', 'DATA Prox', 'DATA Dist', 'Location','NorthEast');
                    legend('REF Prox','DATA Prox', 'Location','NorthEast');
                    axis([0 max(Xfreq) -60 0]); %xlim 0-20MHZ y ylim -60dB 0dB    
                    title ('Final Power Spectrum (max depth)')
                    
                end

            
            % Selection of the useful bandwidth
            tf_ref_pro = tf_ref_pro(band_ind)'./T;
            tf_data_pro = tf_data_pro(band_ind);
            tf_ref_dis = tf_ref_dis(band_ind)'./T;
            tf_data_dis = tf_data_dis(band_ind);
            
            % Inter proxi distal length
            zp_zd = Z_subROI_pro - Z_subROI_dis; % [cm]
            clear Z_subROI_pro Z_subROI_dis
            
            % Spectral log difference
            SR_pro = tf_data_pro' ./ tf_ref_pro;
            SR_dis = tf_data_dis' ./ tf_ref_dis;
            SR_term(i,j,:) = (log(SR_pro) - log(SR_dis));
            SLD_term(i,j,:) = (log(SR_pro) - log(SR_dis)) - 4*alpha_phan_freq * zp_zd ;
            
            %                 figure(12456)
            %                 plot( band, squeeze(squeeze(SLD_term(i,j,:))), 'k', 'Linewidth', 2); hold off;
            %                 xlabel('Frequency [MHz]','fontsize', 14);
            %                 ylabel('SLD term [Np]','fontsize', 14);
            %                 %legend('Low diff log', 'Low fit', 'High diff fit', 'High fit', 'Average diff log', 'Average fit');
            %                 axis([band(1) band(end) min(SLD_term(i,j,:)) max(SLD_term(i,j,:))])
            
            
            env_data_block_sam  = env_rfdata_sam_roi(pixi,pixj,:);
            SNR = mean(env_data_block_sam)/std(env_data_block_sam);
            og_SNR(i,j) = SNR;
            delta_SNR(i,j) = abs(SNR - 1.91).*100./1.91;
            delta_SNR_noAbs(i,j) = (SNR - 1.91).*100./1.91;
        end
        
    end
            figure(171), legend('Top SAM','Top REF','Half SAM','Half REF', 'Bottom SAM','Bottom REF');
            figure(172), legend('Top SAM','Top REF','Half SAM','Half REF', 'Bottom SAM','Bottom REF');
            
    band_T = band';        
    A1 = kron( 4*-zp_zd*band_T , speye(nb_ROI_lateral*nb_ROI_axial) );
    A2 = kron( ones(size(band_T)) , speye(nb_ROI_lateral*nb_ROI_axial) );
    A = [A1 A2];
                
%   figure, imagesc()
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%% PACKAGE EMZ  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SLD.SLD_term = SLD_term;
   SLD.band = band_T;
   SLD.A = A;
   SLD.A1 = A1;
   SLD.A2 = A2;
   SLD.x_ori = DATA.x;
   SLD.z_ori = DATA.z;

   SLD.zp_zd = zp_zd;
   SLD.DATA_ROI = DATA_ROI_seg;
   SLD.DATA_Bmode_ROI = my_RF2Bmode(DATA_ROI_seg);
   
   SLD.DATA_Bmode = my_RF2Bmode(DATA.rf);

   % SNR H. CHAHUARA'S CODE
   SLD.delta_SNR = delta_SNR;
   SLD.og_SNR = og_SNR;
   SLD.delta_SNR_noAbs = delta_SNR_noAbs;

   SLD.z =  Zaxis; % OK
   SLD.x =  Xaxis; % OK

   SLD.z_v2 = Z_ROI; % NOT OK
   SLD.x_v2 = X_ROI; % NOT OK

    %   %%%%%%%%%%%%% FOR WEIGHTS THROUGH FREQUENCY %%%%%%%%%%%%%  
    results_multifreq = zeros(size(SLD.SLD_term,3),3);
    muX = zeros(1,size(SLD.SLD_term,3));
        for i = 1 : size(SLD.SLD_term,3)

            temp = SLD.SLD_term(:,:,i); % SLD IMAGE
            temp = temp(:);
            temp(temp==0)=NaN;

            results_multifreq(i,1) = mean(temp, 'omitnan');
            results_multifreq(i,2) = std(temp, 'omitnan');
            results_multifreq(i,3) = results_multifreq(i,2)/abs(results_multifreq(i,1));

            muX(i) = results_multifreq(i,3);
                   
        end
        SLD.muX = muX;
    figure, 
    plot(SLD.band', SLD.muX, 'b.-'), title('Weights Julien'), grid, 
    xlabel('Frecuency [MHz]'), ylabel('Weights 1/SNR');
%     %%%%%%%%%%%%% FOR WEIGHTS THROUGH FREQUENCY %%%%%%%%%%%%%   

    %     % OPTIMIZED January 2023
    reshapedSLD = reshape(SLD.SLD_term, [], size(SLD.SLD_term, 3));
    meanValues = mean(reshapedSLD, 1, 'omitnan');
    stdValues = std(reshapedSLD, 0, 1, 'omitnan');
    
    ratios = abs(meanValues) ./ stdValues;
    SLD.SNR_ratios = ratios;
    SLD.SR_term = SR_term;
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