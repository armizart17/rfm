function spectralData = calc_powerSpectra_vSimple(data, pars)
% function spectralData = calc_powerSpectra_vSimple(data, pars)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION
% Provides the PowerSpectra (S) and compensates it accordingly (type SLD)
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
%       spectralData:  Structure with the following fields:
%             - powerSpectra    : Power spectrum (rows, cols, freqChannels)
%             - depth           : depth coordinates for PowerSpectra   [m]
%             - lateral         : lateral coordinates for PowerSpectra [m]
%             - band            : band of frequencies [MHz]
%             - delta_snr       : Delta SNR
%             - rf_roi          : RF ROI
% -------------------------------------------------------------------------
% AUTHOR: EMZ based on LIM repository SLD method 
% Note: Remember PowerLaw 
% SR        = log(S_sam ./ S_ref); (rows, cols, freqChannels)
% SR_rpltv  = permute(SR, [3,1,2]); (freqChannels, rows, cols)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Extra specific position spectrum
    if isfield(pars, 'z_target')
        z_target = pars.z_target;
    end
    if isfield(pars, 'x_target')
        x_target = pars.x_target;
    end

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
    
    %% Cropping and finding sample sizes
    % Limits for ACS estimation
    ind_x = x_ini <= x & x <= x_end;
    ind_z = z_ini <= z & z <= z_end;
    x = x(ind_x);
    z = z(ind_z);
    rfdata_sam_roi = rfdata_sam(ind_z,ind_x);
    
    % Lateral samples
    % wx = round(blocksize(1)*(1-overlap)/dx);  % Between windows OLD
    % nx = round(blocksize(1)/dx);  % Window size OLD
    wx = round(blocksize_wv*lambda*(1-overlap)/dx);  % Between windows  
    nx = round(blocksize_wv*lambda/ dx);

    x0 = 1:wx:length(x)-nx;
    x_ACS = x(1,x0+round(nx/2));
    n  = length(x0);
    
    % Axial samples
    % wz = round(blocksize(2)*(1-overlap)/dz); % Between windows OLD
    % nz = 2*round(blocksize(2)/dz /2); % Window size OLD
    % ratio_zx = 1.5; % ^
    % ratio_zx = 1.5; % @
    ratio_zx = 1;
    if isfield(ratio_zx, 'pars')
        ratio_zx = pars.ratio_zx;
    end
    wz = round(blocksize_wv*lambda*(1-overlap)/dz * ratio_zx); % Between windows
    nz = 2*round(blocksize_wv*lambda/dz /2 * ratio_zx); % Window size
   
    z0 = 1:wz:length(z)-nz;
    z_ACS = z(z0+nz/2);
    m  = length(z0);
     
    % Frequency samples
    % NFFT = 2^(nextpow2(nz/2)+1); % old
    NFFT = 2^(nextpow2(nz/1)+1);
    % axis_f = (0:(NFFT/2-1))'*fs/NFFT;  % [Hz] (so 0-fs/2),  it should be
    axis_f = (0:NFFT-1)'/NFFT * fs;   % [Hz] freq axis as default because "spectra" function
    freq_L = bw(1)*1E6; % [Hz] 
    freq_H = bw(2)*1E6; % [Hz]

    ind_f = axis_f >= freq_L & axis_f <= freq_H ;   

    band  = axis_f(ind_f)*1E-6; % [MHz]
    p = length(band);

    bandFull = axis_f*1E-6; % [MHz]
    pFull = length(bandFull);
    
       
    fprintf('Power Spectra vSimple');
    fprintf('\nBandwith: %.2f - %.2f MHz\n',freq_L*1e-6,freq_H*1e-6)
    fprintf('Blocksize in wavelengths: %i\n',blocksize_wv)
    fprintf('Blocksize x: %.2f mm, z: %.2f mm\n',nx*dx*1e3,nz*dz*1e3)
    fprintf('Blocksize in pixels nx: %i, nz: %i\n',nx,nz);
    fprintf('Region of interest columns: %i, rows: %i\n\n',m,n);
    
    %%
    % Saran Layer compensation 
    t_saran = 1;
    if (saran_layer)
        t_saran = transmission_saran(band);
    end

    %% Spectrum
    % Windows for spectrum
    % windowing = tukeywin(nz/2,0.25); % OLD
    % windowing = windowing*ones(1,nx); % OLD

    windowing = window_choice(nz, window_type); %@ nz/2 before
    windowing = windowing*ones(1, nx);
    
    nSamples = size(rfdata_sam_roi,3);
    S_ref = zeros(m,n,p,nSamples);

    S_full = zeros(m,n,pFull,nSamples);

    % € 
    SNR = zeros(m,n,nSamples);

    % SPECIFIC POSITIONS
    % [~, i_idx] = min(abs(z_ACS - z_target)); % $
    % [~, j_idx] = min(abs(x_ACS - x_target)); % $

    % loop in case many Acq 
    % (useful if many frames in reference phantom)
    for iRef = 1:nSamples   
        samRef = rfdata_sam_roi(:,:,iRef);

        envelope = abs(hilbert(samRef)); % € 
        for jj=1:n
            for ii=1:m
                xw = x0(jj);   % x window
                zw = z0(ii);
    
                sub_block = samRef(zw:zw+nz-1,xw:xw+nx-1);

                [tempSp,~] = spectra(sub_block,windowing,0,nz,NFFT);

                S_ref(ii,jj,:,iRef) = tempSp(ind_f) ./ t_saran;

                % % $ 
                % if and(ii == i_idx , jj==j_idx)
                %     S_idx = squeeze(S_ref(i_idx,j_idx,:,iRef));
                % end

                % Extra €
                S_full(ii,jj,:,iRef) = tempSp;

                % € 
                sub_env  = envelope(zw:zw+nz-1, xw:xw+nx-1);
                temp = sub_env(:);
                SNR(ii,jj,iRef) = mean(temp) / std(temp);

            end
        end
    end

%% EXTRA TO TEST %%

% % Plot
% [pxx,fpxx] = pwelch(sam1-mean(sam1),nz,nz-wz,nz,fs);
% meanSpectrum = mean(pxx,2);
% meanSpectrum(1) = 0;
% figure,
% plot(fpxx/1E6,db(meanSpectrum/max(meanSpectrum))),grid on
% xlim([0, fs/2E6])
% hold on
% xline(freq_L/1E6, 'k--')
% xline(freq_H/1E6, 'k--')
% hold off
% xlabel('Frequency [MHz]')
% ylabel('Magnitude [dB]')
%%
    S = mean(S_ref, 4); 

    % €
    SNR = mean(SNR, 3);
    SNRopt = sqrt(1/(4/pi - 1)); % 1.91 
    delta_snr = 100*abs(SNR - SNRopt)/SNRopt;

    % Packing
    spectralData.powerSpectra      = S;
    spectralData.depth             = z_ACS(:); % col vector
    spectralData.lateral           = x_ACS(:); % col vector
    spectralData.band              = band(:);  % col vector
    spectralData.delta_snr         = delta_snr;

    spectralData.rf_roi            = rfdata_sam_roi;
    spectralData.z_roi             = z;
    spectralData.x_roi             = x;

    spectralData.Sfull             = mean(S_full, 4);
    spectralData.bandFull          = bandFull(:);  % col vector
    
    % Plot
    % [pxx,fpxx] = pwelch(rfdata_sam_roi-mean(rfdata_sam_roi),nz,nz-wz,nz,fs);
    % meanSpectrum = mean(pxx,2);
    % meanSpectrum(1) = 0;
    % 
    % figure,
    % plot(fpxx/1e6,db(meanSpectrum/max(meanSpectrum))),grid on
    % xlim([0, fs/2e6])
    % hold on
    % xline(freq_L/1e6, 'k--')
    % xline(freq_H/1e6, 'k--')
    % hold off
    % xlabel('Frequency [MHz]')
    % ylabel('Magnitude [dB]')


return

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

%% FULL SPECTRUM 
%% PLOT FULL Spectral
% Get size of Sp
% [m, n, pFull] = size(spectralData_sam.Sfull);
% 
% nLines = 5;
% lin_cen = round(n / 2); 
% lat_range = max(1, lin_cen-fix(nLines/2)):min(n, lin_cen+fix(nLines/2)); 
% 
% Sd_matrix = squeeze(mean(spectralData_sam.Sfull(:, lat_range, :), 2)); % Mean over 2nd dim (lateral)
% 
% freq_Pos   = (0:(pFull/2-1))*SAM.fs/pFull*1e-6; % freq_rangePos = size ( (0:(1/2-1/pFull))'*SAM.fs) ;
% Sd_Pos     = Sd_matrix(:,1:pFull/2);
% 
% Sd_Pos_dB  = pow2db(Sd_Pos ./ max(Sd_Pos, [], 2));
% % Plot the FULL extracted spectra
% figure;
% set(gcf,'units','normalized','outerposition',[0 0.1 0.75 0.5]); box on;
% subplot(1,3,1)
% imagesc(freq_Pos, spectralData_sam.depth*1e3, Sd_Pos);
% xlabel('Frequency [MHz]');
% ylabel('Depth [mm]');
% colorbar;
% title('Power Spectrum');
% 
% subplot(1,3,2)
% imagesc(freq_Pos, spectralData_sam.depth*1e3, Sd_Pos_dB);
% xlabel('Frequency [MHz]');
% ylabel('Depth [mm]');
% h2 = colorbar; 
% ylabel(h2,'dB');
% title('Power Spectrum');
% 
% subplot(1,3,3)
% plot(freq_Pos, Sd_Pos_dB(1, :), 'DisplayName', 'Top')
% hold on, grid on;
% plot(freq_Pos, Sd_Pos_dB(round(m/2), :), 'DisplayName', 'Half')
% plot(freq_Pos, Sd_Pos_dB(end, :), 'DisplayName', 'Bottom')
% yline(-20, 'k--')
% xlim([0 15]);
% hold off;
% xlabel('Frequency [MHz]');
% ylabel('NormMax [dB]');
% title('Power Spectrum');
% legend('Location', 'Best');
% 
% % FIND PEAKS at all depths
% 
% mDepths = size(Sd_Pos, 1);
% f = freq_Pos;
% ratio = db2mag(-20);
% arrayBands = zeros(mDepths,  2);
% for dd = 1:mDepths
%     y = Sd_Pos(dd, :);
%     [fLeft,fRight] = findFreqBand(f, y, ratio);
%     arrayBands(dd,1) = fLeft;
%     arrayBands(dd,2) = fRight;
% end
% 
% figure,
% imagesc(freq_Pos, spectralData_sam.depth*1e3, Sd_Pos_dB);
% hold on
% plot(arrayBands(:,1), spectralData_sam.depth * 1e3, 'r', 'LineWidth', 2); % Left boundary
% plot(arrayBands(:,2), spectralData_sam.depth * 1e3, 'r', 'LineWidth', 2); % Right boundar
% hold off
% xlabel('Frequency [MHz]');
% ylabel('Depth [mm]');
% h2 = colorbar; 
% ylabel(h2,'dB');
% title('Norm Power Spectrum');
