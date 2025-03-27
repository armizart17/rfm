% Reference Frequency Method v2 March
% AMZ 

% clear all, clc, close all;

Np2dB       = 20*log10(exp(1));
dB2Np       = 1/Np2dB;
range_bmode = [-60 0];

mean2d = @(x) mean(x(:));
std2d = @(x) std(x(:));
cv2d = @(x) 100*std(x(:))/mean(x(:));

calc2dStats = {@(x) mean(x(:)), @(x) std(x(:)), @(x) 100 * std(x(:)) / mean(x(:))};

%% LOAD SAM

% DATA NEW AMZ
pathData = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\dataACS_kwave';

%%%%%%%%%%%%%%% NEW MARCH %%%%%%%%%%%%%%%
alpha_sam = 0.7; % ACS 0.4 0.5 0.6 0.7 1
folderDataSam = 'a_pow1p';
rf_sam_name = sprintf('rf1_a_%.2g', alpha_sam);
% rf_sam_name = strrep(rf_sam_name, '.', 'p');
% rf_sam_name = 'sam0';
SAM = load(fullfile(pathData, folderDataSam, rf_sam_name + ".mat"));

alpha_ref = 0.4; % ACS 0.4 0.5 0.6 0.7 1
folderDataRef = 'a_pow1p';
rf_ref_name = sprintf('rf1_a_%.2g', alpha_ref);
% rf_ref_name = strrep(rf_ref_name, '.', 'p');
% rf_ref_name = 'ref0';
REF = load(fullfile(pathData, folderDataRef, rf_ref_name + ".mat"));
%%%%%%%%%%%%%%% NEW MARCH %%%%%%%%%%%%%%%

gt_acs  = alpha_sam;

% B-MODE CHECK
bmode_sam = db(hilbert(SAM.rf));
bmode_sam = bmode_sam - max(bmode_sam(:));

figure,
% subplot(121), 
imagesc(SAM.x*1E3, SAM.z*1E3, bmode_sam, range_bmode), axis("image");
xlabel('Lateral [mm]'), ylabel('Depth [mm]');
cb = colorbar;
cb.Label.String = 'dB'; % Add the label "dB"
title('SAM')
colormap('gray')

%% SPECTRAL PARAMETERS
pars.bw          = [3 9]; % [MHz]
pars.overlap     = 0.8;
pars.blocksize   = 10; % wavelengths
pars.z_roi       = [5 45]*1E-3; % [m] 
pars.x_roi       = [-17 17]*1E-3; % [m] 
pars.saran_layer = false;
pars.ratio_zx    = 1.25;
pars.window_type = 3; %  (1) Hanning, (2) Tuckey, (3) Hamming, (4) Tchebychev


blocksize_wv_r = 12;
%% RFM V1 
% Reading experiment settings parameters
bw              = pars.bw;
overlap         = pars.overlap;
blocksize_wv    = pars.blocksize;
z_ini           = pars.z_roi(1);
z_end           = pars.z_roi(2);
x_ini           = pars.x_roi(1);
x_end           = pars.x_roi(2);
saran_layer     = pars.saran_layer;
ratio_zx        = pars.ratio_zx;

% Spectral Ratio INIT
z            = SAM.z;
x            = SAM.x;
fs           = SAM.fs;

% Default window
window_type     = 2; %  (1) Hanning, (2) Tukey, (3) Hamming, (4) Tchebychev
if isfield(pars, 'window_type')
    window_type = pars.window_type;
end

% Reading phantoms parameters
if isfield(SAM, 'rf')
    rfdata_sam   = SAM.rf;
end
if isfield(SAM, 'RF')
    rfdata_sam   = SAM.RF;
end

% figure,
% imagesc(SAM.x*1E3, SAM.z*1E3, bmode_sam), axis("image");
% rectangle('Position', 1E3*[pars.x_roi(1) pars.z_roi(1) pars.x_roi(2)-pars.x_roi(1) pars.z_roi(2)-pars.z_roi(1)], ...
%         'EdgeColor','r', 'LineWidth', 2, 'LineStyle','--'), hold off;
% xlabel('Lateral [mm]'), ylabel('Depth [mm]');
% title('SAM')
% colormap('gray')

dx = x(2)-x(1);
dz = z(2)-z(1);
c0 = 1540; % [m/s]
lambda = c0/mean(bw)*1E-6;

% Cropping and finding sample sizes
% Limits for ACS estimation
ind_x = x_ini <= x & x <= x_end;
ind_z = z_ini <= z & z <= z_end;
x = x(ind_x);
z = z(ind_z);
rfdata_sam_roi = rfdata_sam(ind_z,ind_x);

% Lateral samples
wx = round(blocksize_wv*lambda*(1-overlap)/dx);  % Between windows  
nx = round(blocksize_wv*lambda/ dx);

x0 = 1:wx:length(x)-nx;
x_ACS = x(x0+round(nx/2));
x_ACS = x_ACS(:); %Ensure col vector
n  = length(x0);

% Axial samples
wz = round(blocksize_wv*lambda*(1-overlap)/dz * ratio_zx); % Between windows
nz = 2*round(blocksize_wv*lambda/dz /2 * ratio_zx); % Window size

z0 = 1:wz:length(z)-nz;
z_ACS = z(1, z0+round(nz/2));
z_ACS = z_ACS(:); %Ensure col vector
m  = length(z0);

% Frequency samples
% NFFT = 2^(nextpow2(nz/2)+1);
NFFT = 2^(nextpow2(nz/2)+2); %**

axis_f = (0:NFFT-1)'/NFFT * fs;   % [Hz] freq axis as default because "spectra" function
freq_L = bw(1)*1E6; % [Hz] 
freq_H = bw(2)*1E6; % [Hz]

ind_f = axis_f >= freq_L & axis_f <= freq_H ;   

band  = axis_f(ind_f)*1E-6; % [MHz]
bandFull  = axis_f*1E-6; % [MHz]

p = length(band);

fprintf('\nBandwith: %.2f - %.2f MHz\n',freq_L*1e-6,freq_H*1e-6)
fprintf('Blocksize in wavelengths: %i\n',blocksize_wv)
fprintf('Blocksize x: %.2f mm, z: %.2f mm\n',nx*dx*1e3,nz*dz*1e3)
fprintf('Blocksize in pixels nx: %i, nz: %i\n',nx,nz);
fprintf('Region of interest columns: %i, rows: %i\n\n',m,n);

% Saran Layer compensation 
t_saran = 1;
if (saran_layer)
    t_saran = transmission_saran(band);
end

% Windows for spectrum

windowing = window_choice(nz, window_type); %@ nz/2 before
windowing = windowing*ones(1,nx);

nSamples = size(rfdata_sam_roi,3);

% € 
SNR = zeros(m,n,nSamples);

% Spectrum Sample
Sp_k = zeros(m,n,NFFT);

samRef = rfdata_sam_roi;

% Spectrum
for jj=1:n
    for ii=1:m
        xw = x0(jj) ;   % x window
        zw = z0(ii);

        sub_block = samRef(zw:zw+nz-1,xw:xw+nx-1);

        [tempSp,~] = spectra(sub_block,windowing,0,nz,NFFT); %@ nz/2 before

        Sp_k(ii,jj,:) = tempSp();
     
    end
end

% RFM
RSp_k = zeros(m,n,NFFT);
RSp_k(:,:, 2:end) = Sp_k(:,:, 2:end) ./ Sp_k(:,:, 1:end-1); clear Sp_k
% For the first slice, keep the ratio the same as the first slice
RSp_k(:,:, 1) = RSp_k(:,:, 2); % Assuming the first slice ratios are 1 as there is no "i-1"

%% Reference Depth

% blocksize_wv_r = blocksize_wv *2;
blocksize_wv_r = blocksize_wv_r;

wz_r = round(blocksize_wv_r*lambda*(1-overlap)/dz * ratio_zx); % Between windows
nz_r = 2*round(blocksize_wv_r*lambda/dz /2 * ratio_zx); % Window size

z0_r = 1:wz_r:length(z)-nz_r;
z_ACS_r = z(z0_r+ nz/2);
z_ACS_r = z_ACS_r(:); %Ensure col vector
m_r  = length(z0_r);

% Frequency samples
% NFFT = 2^(nextpow2(nz/2)+1);
NFFT = 2^(nextpow2(nz/2)+2); %**

axis_f = (0:NFFT-1)'/NFFT * fs;   % [Hz] freq axis as default because "spectra" function
freq_L = bw(1)*1E6; % [Hz] 
freq_H = bw(2)*1E6; % [Hz]

ind_f = axis_f >= freq_L & axis_f <= freq_H ;   

band  = axis_f(ind_f)*1E-6; % [MHz]
bandFull  = axis_f*1E-6; % [MHz]

fprintf('\nBandwith: %.2f - %.2f MHz\n',freq_L*1e-6,freq_H*1e-6)
fprintf('Blocksize in wavelengths: %i\n',blocksize_wv)
fprintf('Blocksize x: %.2f mm, z: %.2f mm\n',nx*dx*1e3,nz*dz*1e3)
fprintf('Blocksize in pixels nx: %i, nz: %i\n',nx,nz);
fprintf('Region of interest columns: %i, rows: %i\n\n',m,n);

% Windows for spectrum
windowing = window_choice(nz, window_type); %@ nz/2 before
windowing = windowing*ones(1,nx);

% Spectrum Sample
Sp_r = zeros(m_r,n,NFFT);

samRef = rfdata_sam_roi;

% figure, 
% set(gcf, 'Position', [0 0 1 1], 'Units', 'normalized');
font = 12;
for jj=1:n
    for ii=1:m_r
        xw = x0(jj) ;   % x window
        zw = z0_r(ii);

        sub_block = samRef(zw:zw+nz-1,xw:xw+nx-1);

        [tempSp,~] = spectra(sub_block,windowing,0,nz,NFFT); %@ nz/2 before

        Sp_r(ii,jj,:) = tempSp();
     
    end
end

% RFM
RSp_r = zeros(m_r,n,NFFT);
RSp_r(:,:, 2:end) = Sp_r(:,:, 2:end) ./ Sp_r(:,:, 1:end-1); clear Sp_r
% For the first slice, keep the ratio the same as the first slice
RSp_r(:,:, 1) = RSp_r(:,:, 2); % Assuming the first slice ratios are 1 as there is no "i-1"

%% ATTEMPT RFM A_local UFR

% UFR strategy
bw_ufr = [3 9];
freqL = bw_ufr(1); freqH = bw_ufr(2);
range = bandFull >= freqL & bandFull <= freqH;

RSp_k_ufr   = RSp_k(:,:,range);
RSp_r_ufr   = RSp_r(:,:,range);

band_ufr    = bandFull(range);
p_ufr       = length(band_ufr);

% Convert depth values to cm
z_ACS_cm = z_ACS * 1e2;      % Convert from meters to cm
z_ACS_r_cm = z_ACS_r * 1e2;  % Convert reference depths to cm

% Delta MHz 
df_MHz = band_ufr(2) - band_ufr(1);

% Preallocate cell arrays for storing x_temp and y_temp
x_temp_all = cell(m, n);
y_temp_all = cell(m, n);

countX = 0; countY = 0;
tic;
for jj = 1:n  % Loop over lateral positions (x_j)
    for ii = 1:m  % Loop over depth positions (z_k)
        
        % Temporary storage for this location
        y_temp = nan(m_r, p_ufr);  % (Reference depth, Frequency)
        x_temp = nan(m_r, p_ufr);  % (Reference depth, Frequency)
        
        for iFreq = 1:p_ufr  % Loop over frequency bins
            for r = 1:m_r  % Loop over reference depths

                % WAY WITH UNITS
                y = ( log(RSp_k_ufr(ii, jj, iFreq)) - log(RSp_r_ufr(r, jj, iFreq)) ) / (4*df_MHz)*Np2dB;
                % if  y == 0
                %     fprintf('ii=%d, jj=%d, r=%d, iFreq=%d \n', ii, jj, r, iFreq); countY = countY +1;
                % end
                X = z_ACS_cm(ii) - z_ACS_r_cm(r);
                % if  X == 0
                %     fprintf('ii=%d, jj=%d, r=%d, iFreq=%d \n', ii, jj, r, iFreq); countX = countX +1;
                % end

                % Store y and X values in the temporary matrices

                y_temp(r, iFreq) = y;  % (Frequency, Reference depth)
                x_temp(r, iFreq) = X;  % (Frequency, Reference depth)
            end

            %%%%%%%%%%%%%%%% TV Denoising %%%%%%%%%%%%%%%%
            % mu = 5;
            % tol = 1e-4;
            % y_og = y_temp(iFreq, :); 
            % % snr1(i) = mean(Y_row)/std(Y_row);
            % [M, N] = size(y_og(:));
            % [y_opt] = IRLS_TV(y_og(:),speye(M*N),mu, M,N,tol,ones(size(M*N)),ones(M*N,1));
            % y_temp(iFreq, :) = y_opt';
            %%%%%%%%%%%%%%%% TV Denoising %%%%%%%%%%%%%%%%
        end
        
        %%%%%%%%%%%%%%%%%%%%%% PLOT FPR %%%%%%%%%%%%%%%%%%%%%%
        % if (ii==fix(m/2) && jj==fix(n/6)) || (ii==fix(m/2) && jj==fix(n/2)) || (ii==fix(m/2) && jj==fix(5*n/6)) % different depths
        if (jj==fix(n/2) && ii==fix(m/6)) || (jj==fix(n/2) && ii==fix(m/2)) || (jj==fix(n/2) && ii==fix(5*m/6)) % different laterals
            freq1 = 3.5; freq2 = 6; freq3 = 8.5;
            idx_f1 = find(band_ufr >= freq1, 1, 'first'); idx_f2 = find(band_ufr >= freq2, 1, 'first'); idx_f3 = find(band_ufr >= freq3, 1, 'first');
            
            figure;
            set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1000, 600]); % [x, y, width, height]
                       
            title_f1 = sprintf('Freq %.2fMHz (a=%.2f)', band_ufr(idx_f1), band_ufr(idx_f1));
            plot(x_temp_data(:, idx_f1), y_temp(:, idx_f1), 'r.-', 'DisplayName', title_f1 ); 
            hold on; 
            title_f2 = sprintf('Freq %.2fMHz (a=%.2f)', band_ufr(idx_f2), band_ufr(idx_f2));
            plot(x_temp_data(:, idx_f2), y_temp(:, idx_f2), 'b.-', 'DisplayName', title_f2 ); 
            title_f3 = sprintf('Freq %.2fMHz (a=%.2f)', band_ufr(idx_f3), band_ufr(idx_f3));
            plot(x_temp_data(:, idx_f3), y_temp(:, idx_f3), 'k.-', 'DisplayName', title_f3 ); 
            hold off; grid on;
            title(sprintf('Data at ii = %d, jj = %d', ii, jj));
            xlabel('X_t [cm]'); ylabel('RS_{norm} [dB/MHz]');
            legend('Location','best')
        end
        %%%%%%%%%%%%%%%%%%%%%% PLOT FPR %%%%%%%%%%%%%%%%%%%%%%
        x_temp_all{ii, jj} = x_temp;
        y_temp_all{ii, jj} = y_temp;

    end
end
t = toc;
fprintf('Loop way Elapsed time %.2f \n', t);

%% ESTIMATION RFM (OLD FOR PRESAVED TV)

% Prellocate a_local
a_rfm = zeros(m, n); % Preallocate local attenuation matrix (depth x lateral)
tic;
for jj = 1:n  % Loop over lateral positions (x_j)
    for ii = 1:m  % Loop over depth positions (z_k)

        x_temp = x_temp_all{ii, jj};
        y_temp = y_temp_all{ii, jj};

        X_vec = x_temp(:);

        % for i = 1:size(y_temp_data, 1)
        % 
        %     Y_row = y_temp_data(i, :);    
        % 
        %     %%%%%%%%%%%%%%%% TV Denoising %%%%%%%%%%%%%%%%
        %     mu = 5;
        %     tol = 1e-4;
        % 
        %     snr1(i) = mean(Y_row)/std(Y_row);
        %     [M, N] = size(Y_row(:));
        %     [y_opt] = IRLS_TV(Y_row(:),speye(M*N),mu,M,N,tol,ones(size(M*N)),ones(M*N,1));
        %     Y_row = y_opt';
        %     y_temp_data(i, :) = Y_row;
        % 
        %     %%%%%%%%%%%%%%%% TV Denoising %%%%%%%%%%%%%%%%
        % 
        % end
        % %

        %%%%%%%%%%%%%%%%%%%%%% PLOT FPR %%%%%%%%%%%%%%%%%%%%%%
        if (ii==fix(m/2) && jj==fix(n/6)) || (ii==fix(m/2) && jj==fix(n/2)) || (ii==fix(m/2) && jj==fix(5*n/6)) % different depths
        % if (jj==fix(n/2) && ii==fix(m/6)) || (jj==fix(n/2) && ii==fix(m/2)) || (jj==fix(n/2) && ii==fix(5*m/6)) % different laterals
            freq1 = 3.5; freq2 = 6; freq3 = 8.5;
            idx_f1 = find(band_ufr >= freq1, 1, 'first'); idx_f2 = find(band_ufr >= freq2, 1, 'first'); idx_f3 = find(band_ufr >= freq3, 1, 'first');
            
            [slope_f1, ~, ~, ~] = fit_linear(x_temp(idx_f1, :), y_temp(idx_f1, :), 2); 
            [slope_f2, ~, ~, ~] = fit_linear(x_temp(idx_f2, :), y_temp(idx_f2, :), 2); 
            [slope_f3, ~, ~, ~] = fit_linear(x_temp(idx_f3, :), y_temp(idx_f3, :), 2); 

            figure;
            set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1000, 600]); % [x, y, width, height]
          
            title_f1 = sprintf('Freq %.2fMHz (a=%.2f)', band_ufr(idx_f1), slope_f1);
            plot(x_temp(idx_f1, :), y_temp(idx_f1, :), 'r.-', 'DisplayName', title_f1 ); 
            hold on; 
            title_f2 = sprintf('Freq %.2fMHz (a=%.2f)', band_ufr(idx_f2), slope_f2);
            plot(x_temp(idx_f2, :), y_temp(idx_f2, :), 'b.-', 'DisplayName', title_f2 ); 
            title_f3 = sprintf('Freq %.2fMHz (a=%.2f)', band_ufr(idx_f3), slope_f3);
            plot(x_temp(idx_f2, :), y_temp(idx_f3, :), 'k.-', 'DisplayName', title_f3 ); 
            hold off; grid on;
            title(sprintf('Data at ii = %d, jj = %d', ii, jj));
            xlabel('X_t [cm]'); ylabel('RS_{norm} [dB/MHz]');
            legend('Location','best')
            set(gca, 'FontSize', 22)
        end
        %%%%%%%%%%%%%%%%%%%%%% PLOT FPR %%%%%%%%%%%%%%%%%%%%%%

        y_vec = y_temp(:);

        a_rfm(ii, jj) = -(X_vec' * X_vec) \ (X_vec' * y_vec);
        % a_rfm(ii, jj) = - cgs(X_vec' * X_vec, X_vec' * y_vec);

    end
end
t = toc;
fprintf('Loop way Elapsed time %.2f \n', t);
% RFM
[m_a, s_a, cv_a] = deal(calc2dStats{1}(a_rfm), calc2dStats{2}(a_rfm), calc2dStats{3}(a_rfm));
caxis_acs = [0 1.1];
fontSize = 14;

figure, 
imagesc(x_ACS * 1e3, z_ACS * 1e3, a_rfm, caxis_acs); % Convert to mm
axis("image")
colorbar; colormap("turbo")
xlabel('Lateral [mm]');
ylabel('Depth [mm]');
hb2=colorbar; ylabel(hb2,'dB\cdotcm^{-1}\cdotMHz^{-1}')
% title('Local Attenuation Coefficient');
title(sprintf('RFM Local AC (GT= %.2f)\n%.3f $\\pm$ %.3f,  \\%%CV = %.2f', ...
               alpha_sam, m_a, s_a, cv_a), ...
      'Interpreter', 'latex');
set(gca,'fontsize',fontSize)

% Apply constraints to stabilize attenuation estimates
% a_min = 0.1; % Lower bound for attenuation
% a_max = 2.0; % Upper bound for attenuation
% a_local = max(min(a_local, a_max), a_min); % Constrain within bounds

%% USE THIS FROM APRIL AND NOW ON
%%%%%%%%%%%%%%%%%%%% FAST WAY %%%%%%%%%%%%%%%%%%%%

% UFR strategy
bw_ufr = [3 9];
freqL = bw_ufr(1); freqH = bw_ufr(2);
range = bandFull >= freqL & bandFull <= freqH;

RSp_k_ufr   = RSp_k(:,:,range);
RSp_r_ufr   = RSp_r(:,:,range);

band_ufr    = bandFull(range);
p_ufr       = length(band_ufr);

% Convert depth values to cm
z_ACS_cm = z_ACS * 1e2;      % Convert from meters to cm
z_ACS_r_cm = z_ACS_r * 1e2;  % Convert reference depths to cm

% Delta MHz 
df_MHz = band_ufr(2) - band_ufr(1);

% Preallocate cell arrays for storing x_temp and y_temp
x_temp_all = cell(m, n);
y_temp_all = cell(m, n);

tic;
for jj = 1:n  % Loop over lateral positions (x_j)
    for ii = 1:m  % Loop over depth positions (z_k)

        % Temporary storage for this location
        y_temp = nan(m_r, p_ufr);  % (Reference depth, Frequency) ** m_r
        x_temp = nan(m_r, p_ufr);  % (Reference depth, Frequency) ** m_r

        for r = 1:m_r  % Loop over reference depths
            % if (ii==1 && r==1)
            y_col = squeeze( ( log(RSp_k_ufr(ii, jj, :)) - log(RSp_r_ufr(r, jj, :)) ) /(4*df_MHz) *Np2dB ); % p_ufr x 1

            X = z_ACS_cm(ii) - z_ACS_r_cm(r);

            y_temp(r, :) = y_col; %**
            x_temp(r, :) = X; % **

        end
        x_temp_all{ii, jj} = x_temp;
        y_temp_all{ii, jj} = y_temp;

    end
end
t = toc;
fprintf('Loop way Elapsed time %.2f \n', t);
%%%%%%%%%%%%%%%%%%%% FAST WAY %%%%%%%%%%%%%%%%%%%%

%% NOW ESTIMATE
% Prellocate a_local
a_rfm = zeros(m, n); % Preallocate local attenuation matrix (depth x lateral)
tic;
for jj = 1:n  % Loop over lateral positions (x_j)
    for ii = 1:m  % Loop over depth positions (z_k)

        x_temp = x_temp_all{ii, jj};
        y_temp = y_temp_all{ii, jj};

        X_vec = x_temp(:);

        %%%%%%%%%%%%%%%%%%%%% TV DENOISING PER CHANNEL %%%%%%%%%%%%%%%%%%%%%
        % for i = 1:size(y_temp, 2)
        % 
        %     y_col = y_temp(:, i);    
        % 
        %     %%%%%%%%%%%%%%%% TV Denoising %%%%%%%%%%%%%%%%
        %     mu = 5;
        %     tol = 1e-4;
        % 
        %     snr1(i) = mean(y_col)/std(y_col);
        %     [M, N] = size(y_col);
        %     [y_opt] = IRLS_TV(y_col(:),speye(M*N),mu,M,N,tol,ones(size(M*N)),ones(M*N,1));
        %     
        %     y_temp(:, i) = y_opt;
        % 
        %     %%%%%%%%%%%%%%%% TV Denoising %%%%%%%%%%%%%%%%
        % 
        % end
        % 
        %%%%%%%%%%%%%%%%%%%%% TV DENOISING PER CHANNEL %%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%% TNV DENOISING JOINT %%%%%%%%%%%%%%%%%%%%%
        lambda   = 1;  % Regularization parameter
        max_iter = 20;
        tol = 2e-3;
        [y_denoised, cost1, err1, fid1, reg1, iter1] = TNV_regularization(y_temp, ...
            lambda, max_iter, tol);
        y_temp = y_denoised;
        %%%%%%%%%%%%%%%%%%%%% TNV DENOISING JOINT %%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%%%%%%% PLOT FPR %%%%%%%%%%%%%%%%%%%%%%
        if (ii==fix(m/2) && jj==fix(n/6)) || (ii==fix(m/2) && jj==fix(n/2)) || (ii==fix(m/2) && jj==fix(5*n/6)) % different depths
        % if (jj==fix(n/2) && ii==fix(m/6)) || (jj==fix(n/2) && ii==fix(m/2)) || (jj==fix(n/2) && ii==fix(5*m/6)) % different laterals
            freq1 = freqL+0.5; freq2 = 0.5*(freqL+freqH); freq3 = freqH-0.5;
            idx_f1 = find(band_ufr >= freq1, 1, 'first'); idx_f2 = find(band_ufr >= freq2, 1, 'first'); idx_f3 = find(band_ufr >= freq3, 1, 'first');
            

            [slope_f1, ~, ~, ~] = fit_linear(x_temp(:, idx_f1), y_temp(:, idx_f1), 2); 
            [slope_f2, ~, ~, ~] = fit_linear(x_temp(:, idx_f2), y_temp(:, idx_f2), 2); 
            [slope_f3, ~, ~, ~] = fit_linear(x_temp(:, idx_f3), y_temp(:, idx_f3), 2); 

            figure;
            set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1000, 600]); % [x, y, width, height]          

            title_f1 = sprintf('Freq %.2fMHz (a=%.2f)', band_ufr(idx_f1), slope_f1);
            plot(x_temp(:, idx_f1), y_temp(:, idx_f1), 'r.-', 'DisplayName', title_f1 ); 
            hold on; 
            title_f2 = sprintf('Freq %.2fMHz (a=%.2f)', band_ufr(idx_f2), slope_f2);
            plot(x_temp(:, idx_f2), y_temp(:, idx_f2), 'b.-', 'DisplayName', title_f2 ); 
            title_f3 = sprintf('Freq %.2fMHz (a=%.2f)', band_ufr(idx_f3), slope_f3);
            plot(x_temp(:, idx_f3), y_temp(:, idx_f3), 'k.-', 'DisplayName', title_f3 ); 
            hold off; grid on;
            title(sprintf('Data at ii = %d, jj = %d', ii, jj));
            xlabel('X_t [cm]'); ylabel('RS_{norm} [dB/MHz]');
            ylim([-20 20]);
            legend('Location','best')
            set(gca, 'FontSize', 22)
        end
        %%%%%%%%%%%%%%%%%%%%%% PLOT FPR %%%%%%%%%%%%%%%%%%%%%%

        y_vec = y_temp(:);

        a_rfm(ii, jj) = -(X_vec' * X_vec) \ (X_vec' * y_vec);
        % a_rfm(ii, jj) = - cgs(X_vec' * X_vec, X_vec' * y_vec);

    end
end
t = toc;
fprintf('Loop way Elapsed time %.2f \n', t);
% RFM
[m_a, s_a, cv_a] = deal(calc2dStats{1}(a_rfm), calc2dStats{2}(a_rfm), calc2dStats{3}(a_rfm));
caxis_acs = [0 1.1];
fontSize = 14;

figure, 
imagesc(x_ACS * 1e3, z_ACS * 1e3, a_rfm, caxis_acs); % Convert to mm
axis("image")
colorbar; colormap("turbo")
xlabel('Lateral [mm]');
ylabel('Depth [mm]');
hb2=colorbar; ylabel(hb2,'dB\cdotcm^{-1}\cdotMHz^{-1}')
% title('Local Attenuation Coefficient');
title(sprintf('RFM Local AC (GT= %.2f)\n%.3f $\\pm$ %.3f,  \\%%CV = %.2f', ...
               alpha_sam, m_a, s_a, cv_a), ...
      'Interpreter', 'latex');
set(gca,'fontsize',fontSize)




%%
keyboard
%% CLASSICAL SLD  METHOD

% SAMPLE
spectralData_sam = calc_powerSpectra_prox_dis(SAM, pars);

Sp_sam = spectralData_sam.Sp;
Sd_sam = spectralData_sam.Sd;

% REFERENCE
spectralData_ref = calc_powerSpectra_prox_dis(REF, pars);

band   = spectralData_sam.band; % [MHz]
zd_zp  = spectralData_sam.zd_zp; % [m]
Sp_ref = spectralData_ref.Sp;
Sd_ref = spectralData_ref.Sd;

% COMPENSATION ATTENUATION
acs_ref     = alpha_ref; % [dB/cm/MHz]
att_ref     = acs_ref*band/Np2dB; % vector [Np/cm]
att_ref_map = reshape(att_ref, 1, 1, []); % 3D array as SLogRatio

SLogRatio = log(Sp_sam./Sd_sam) - ( log(Sp_ref./Sd_ref) - 4*zd_zp*1E2*att_ref_map ); % (rows, cols, freqChannels)
clear('Sp_sam','Sd_sam', 'Sp_ref', 'Sd_ref', 'att_ref_map', 'att_ref');

[m, n, ~] = size(SLogRatio);
A1 = kron( 4*zd_zp*1E2*band , speye(m*n) );
A2 = kron( ones(size(band)) , speye(m*n) );
[ACS_SLD, ~] = cgs_ACS([A1 A2], SLogRatio);

%% FIGURES COMPARISON RFM VS SLD

caxis_acs = [0 1.1];
fontSize = 14;

% ACS_SLD = ACS_RSLD;
figure, 
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1000, 600]); % [x, y, width, height]

% SLD
[m_a, s_a, cv_a] = deal(calc2dStats{1}(ACS_SLD), calc2dStats{2}(ACS_SLD), calc2dStats{3}(ACS_SLD));

subplot(121)
imagesc(x_ACS * 1e3, z_ACS * 1e3, ACS_SLD, caxis_acs); % Convert to mm
axis("image")
colorbar; colormap("turbo")
xlabel('Lateral [mm]');
ylabel('Depth [mm]');
hb2=colorbar; ylabel(hb2,'dB\cdotcm^{-1}\cdotMHz^{-1}')
% title('Local Attenuation Coefficient');
title(sprintf('SLD Local AC (GT= %.2f)\n%.3f $\\pm$ %.3f,  \\%%CV = %.2f', ...
               alpha_sam, m_a, s_a, cv_a), ...
      'Interpreter', 'latex');
set(gca,'fontsize',fontSize)

% RFM
[m_a, s_a, cv_a] = deal(calc2dStats{1}(a_rfm), calc2dStats{2}(a_rfm), calc2dStats{3}(a_rfm));

subplot(122)
imagesc(x_ACS * 1e3, z_ACS * 1e3, a_rfm, caxis_acs); % Convert to mm
axis("image")
colorbar; colormap("turbo")
xlabel('Lateral [mm]');
ylabel('Depth [mm]');
hb2=colorbar; ylabel(hb2,'dB\cdotcm^{-1}\cdotMHz^{-1}')
% title('Local Attenuation Coefficient');
title(sprintf('RFM Local AC (GT= %.2f)\n%.3f $\\pm$ %.3f,  \\%%CV = %.2f', ...
               alpha_sam, m_a, s_a, cv_a), ...
      'Interpreter', 'latex');
set(gca,'fontsize',fontSize)



%%
keyboard


%% IMOVERLAY SIMPLE

% Simple Script example to overlay RESIZED colorImg to bmodeFull
% Function requires previous resizing (i.e. bigImg)

units           = 1E3;
bmodeFull       = bmode_sam;
colorImg        = bigImg(a_local, rfdata_sam_roi );   
range_bmode     = [-60 0];
range_img       = [0 1.1];
transparency    = 0.65;
x_img           = x_ACS*units;
z_img           = z_ACS*units;
xFull           = SAM.x*units;
zFull           = SAM.z*units;

figure, 

[~,hB,hColor] = imOverlaySimple(bmodeFull, colorImg, range_bmode, ...
                range_img, transparency, x_img, z_img, xFull, zFull);

xlabel('Lateral [mm]'), ylabel('Depth [mm]');
hColor.Label.String = 'dB\cdotcm^{-1}\cdotMHz^{-1}';
title(sprintf('Local AC (GT= %.2f)\n%.3f $\\pm$ %.3f,  \\%%CV = %.2f', ...
               alpha_sam, m_a, s_a, cv_a), ...
      'Interpreter', 'latex');

set(gca,'fontsize',fontSize)

%% SAVE FPR terms WITH FIT UPDATE Feb W3

[m, n, p] = size(RSp_k_ufr);
% Make full ACS all freq
z_ACS_cm = z_ACS*1E2;
df_MHz = min(diff(band));

% Initialize arrays to store FPR, linear fits, and slopes
FPR_all     = zeros(length(z_ACS_cm), n, p); % Store original FPR
FPR_fit_all = zeros(length(z_ACS_cm), n, p); % Store linear fit (lin_y)
slopes_all  = zeros(n, p); % Store slope values
ACS_RFM     = zeros(p, n); 
for ff = 1:p 
    for jj = 1:n
        % Compute FPR
        % FPR = log(squeeze(squeeze(RSp(:,jj,ff)))) / (4*df_MHz)*Np2dB;
        FPR = log(squeeze(squeeze(RSp_k_ufr(:,jj,ff)))) / (4*df_MHz)*Np2dB;
        
        % Perform linear fit
        [lin_slope , lin_intercept, lin_y, R_fit] = fit_linear(z_ACS_cm, FPR, 2); 
        
        % Store FPR, linear fit, and slope
        FPR_all(:, jj, ff) = FPR;
        FPR_fit_all(:, jj, ff) = lin_y;
        slopes_all(jj, ff) = lin_slope;
        
        % Store ACS value
        % ACS_RFM(ff, jj) = -lin_slope; % -8.68*m/(4*delta_f) [dB/MHz/cm]
        ACS_RFM(ff, jj) = +lin_slope; % -8.68*m/(4*delta_f) [dB/MHz/cm]
    end
end

%%
% RESULTS IMAGESC
acs_range = [-1.2 1.2];
figure, 
imagesc(x_ACS*1E2, band, ACS_RFM, acs_range), colorbar
% hb2=colorbar; ylabel(hb2,'dB/cm/MHz')
hb2=colorbar; ylabel(hb2,'dB\cdotcm^{-1}\cdotMHz^{-1}')
xlabel('Lateral [cm]')
ylabel('Freq [MHz]')

acs_range = [0 1.2];
figure, 
imagesc(x_ACS*1E2, band, abs(ACS_RFM), acs_range), colorbar
colormap("turbo")
% hb2=colorbar; ylabel(hb2,'dB/cm/MHz')
hb2=colorbar; ylabel(hb2,'dB\cdotcm^{-1}\cdotMHz^{-1}')
xlabel('Lateral [cm]')
ylabel('Freq [MHz]')

freq1 = 3.5; freq2 = 6; freq3 = 8.5;
idx_f1 = find(band >= freq1, 1, 'first');
idx_f2 = find(band >= freq2, 1, 'first');
idx_f3 = find(band >= freq3, 1, 'first');

figure, 
plot(x_ACS', ACS_RFM(idx_f1,:), 'r'), hold on
plot(x_ACS', ACS_RFM(idx_f2,:), 'g'), hold on
plot(x_ACS', ACS_RFM(idx_f3,:), 'b'), hold on
yline(SAM.alpha_value, 'k--')

f_vs_acs= mean((ACS_RFM), 2);
figure, 
plot(band, f_vs_acs)
yline(mean(f_vs_acs), 'b-')
yline(SAM.alpha_value, 'k--')
xlabel('Freq')
ylabel('ACS')


%% WITH SLOPE PRINT

% Find freqs MHz
freq1 = 3.5; freq2 = 6; freq3 = 8.5;
idx_f1 = find(band >= freq1, 1, 'first');
idx_f2 = find(band >= freq2, 1, 'first');
idx_f3 = find(band >= freq3, 1, 'first');

% Define the range of central horizontal positions to average
pos_range = 1:n; % Selected range of horizontal positions

% Compute the mean FPR and the mean of the linear fit over the selected range
FPR_avg_f1 = mean(FPR_all(:, pos_range, idx_f1), 2); 
FPR_avg_f2 = mean(FPR_all(:, pos_range, idx_f2), 2); 
FPR_avg_f3 = mean(FPR_all(:, pos_range, idx_f3), 2); 

FPR_fit_avg_f1 = mean(FPR_fit_all(:, pos_range, idx_f1), 2);
FPR_fit_avg_f2 = mean(FPR_fit_all(:, pos_range, idx_f2), 2);
FPR_fit_avg_f3 = mean(FPR_fit_all(:, pos_range, idx_f3), 2);

% Compute the average slope over the selected positions
slope_avg_f1 = mean(slopes_all(pos_range, idx_f1));
slope_avg_f2 = mean(slopes_all(pos_range, idx_f2));
slope_avg_f3 = mean(slopes_all(pos_range, idx_f3));

% Plot the results
figure,
sgtitle('Frequency Power Ratio')
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 900, 700]); % [x, y, width, height] in pixels

% Plot raw FPR data
plot(z_ACS_cm, FPR_avg_f1, 'r-', 'DisplayName', sprintf('S_{%.2fMHz}/S_{%.2fMHz} (Slope: %.3f)', band(idx_f1), band(idx_f1-1), slope_avg_f1), 'LineWidth', 1.5), hold on
plot(z_ACS_cm, FPR_avg_f2, 'b-', 'DisplayName', sprintf('S_{%.2fMHz}/S_{%.2fMHz} (Slope: %.3f)', band(idx_f2), band(idx_f2-1), slope_avg_f2), 'LineWidth', 1.5), hold on
plot(z_ACS_cm, FPR_avg_f3, 'k-', 'DisplayName', sprintf('S_{%.2fMHz}/S_{%.2fMHz} (Slope: %.3f)', band(idx_f3), band(idx_f3-1), slope_avg_f3), 'LineWidth', 1.5), hold on

% Plot linear fits WITHOUT DisplayName
plot(z_ACS_cm, FPR_fit_avg_f1, 'r--', 'LineWidth', 2), hold on
plot(z_ACS_cm, FPR_fit_avg_f2, 'b--', 'LineWidth', 2), hold on
plot(z_ACS_cm, FPR_fit_avg_f3, 'k--', 'LineWidth', 2), hold on
legend ('Location', 'best');

% Plot raw FPR data
% plot(z_ACS_cm, FPR_avg_f1, 'r-', 'DisplayName', sprintf('S_{%.2fMHz}/S_{%.2fMHz}', band(idx_f1), band(idx_f1-1)), 'LineWidth', 1.5), hold on
% plot(z_ACS_cm, FPR_avg_f2, 'b-', 'DisplayName', sprintf('S_{%.2fMHz}/S_{%.2fMHz}', band(idx_f2), band(idx_f2-1)), 'LineWidth', 1.5), hold on
% plot(z_ACS_cm, FPR_avg_f3, 'k-', 'DisplayName', sprintf('S_{%.2fMHz}/S_{%.2fMHz}', band(idx_f3), band(idx_f3-1)), 'LineWidth', 1.5), hold on

% 
% % Plot linear fits with slope values in the legend
% plot(z_ACS_cm, FPR_fit_avg_f1, 'r--', 'DisplayName', sprintf('Fit: S_{%.2fMHz}/S_{%.2fMHz} (Slope: %.3f)', band(idx_f1), band(idx_f1-1), slope_avg_f1), 'LineWidth', 2), hold on
% plot(z_ACS_cm, FPR_fit_avg_f2, 'b--', 'DisplayName', sprintf('Fit: S_{%.2fMHz}/S_{%.2fMHz} (Slope: %.3f)', band(idx_f2), band(idx_f2-1), slope_avg_f2), 'LineWidth', 2), hold on
% plot(z_ACS_cm, FPR_fit_avg_f3, 'k--', 'DisplayName', sprintf('Fit: S_{%.2fMHz}/S_{%.2fMHz} (Slope: %.3f)', band(idx_f3), band(idx_f3-1), slope_avg_f3), 'LineWidth', 2), hold on

grid on
xlabel('Depth [cm]')
ylabel('FPR [dB/MHz]')
ylabel('FPR [dB\cdotMHz^{-1}]')

set(gca, 'FontSize', 14)




%% HISTOGRAM

% figure;
% histogram(ACS_RFM, 20, 'Normalization', 'probability'); % 20 bins
% xlabel('Value');
% ylabel('Frequency');
% title('Histogram of ACS Elements');
% grid on;

%% BOX PLOT RESULTS RFM

% Compute mean and standard deviation
mu      = mean(ACS_RFM(:));
sigma   = std(ACS_RFM(:));
cv      = sigma/mu;

% Generate the box plot
figure;
boxplot(ACS_RFM(:)); 
grid on;
yline(gt_acs, 'k--', 'LineWidth', 1.5); % Dashed line for GT ACS
ylim([-5 5]);

% Set title with mean and standard deviation 
title(sprintf('RFM ACS: %.3f \\pm %.3f, %%CV=%.2f', mu, sigma, 100*cv));

% Axis labels
xlabel('ACS');
ylabel('Values');

% Add a text box showing GT ACS
text(0.05, 0.95, sprintf('GT ACS = %.2f', gt_acs), ...
    'Units', 'normalized', 'FontSize', 12, ...
    'Color', 'k', 'BackgroundColor', 'w', ...
    'EdgeColor', 'k', 'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'top');




%% SLD
[m, n, p] = size(SLogRatio);
A1 = kron( 4*zd_zp*1E2*band , speye(m*n) );
A2 = kron( ones(size(band)) , speye(m*n) );
[ACS_SLD, ~] = cgs_ACS([A1 A2], SLogRatio);

% Compute mean and standard deviation
mu      = mean(ACS_SLD(:));
sigma   = std(ACS_SLD(:));
cv      = sigma/mu;

% Generate the box plot
figure;
boxplot(ACS_SLD(:)); 
grid on;
yline(gt_acs, 'k--', 'LineWidth', 1.5); % Dashed line for GT ACS
% ylim([-5 5]);

% Set title with mean and standard deviation 
title(sprintf('SLD ACS: %.3f \\pm %.3f, %%CV=%.2f', mu, sigma, 100*cv));

% Axis labels
xlabel('ACS');
ylabel('Values');

% Add a text box showing GT ACS
text(0.05, 0.95, sprintf('GT ACS = %.2f', gt_acs), ...
    'Units', 'normalized', 'FontSize', 12, ...
    'Color', 'k', 'BackgroundColor', 'w', ...
    'EdgeColor', 'k', 'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'top');

%%
keyboard
%% OPTIMIZATION RSLD

% Optimization constants
tol = 1e-3;
mu1 = 10^3.5; 
mu2 = 10^3.5;

[m, n, p] = size(SLogRatio);
mask = ones(m,n,p);
A1 = kron( 4*zd_zp*1E2*band , speye(m*n) );
A2 = kron( ones(size(band)) , speye(m*n) );

[Bn,~] = AlterOpti_ADMM(A1,A2,SLogRatio(:),mu1,mu2,m,n,tol,mask(:));
ACS_RSLD = reshape(Bn*Np2dB,m,n);

% Compute mean and standard deviation
mu      = mean(ACS_RSLD(:));
sigma   = std(ACS_RSLD(:));
cv      = sigma/mu;

% % Generate the box plot
% figure;
% boxplot(ACS_RSLD(:)); 
% grid on;
% yline(gt_acs, 'k--', 'LineWidth', 1.5); % Dashed line for GT ACS
% % ylim([-5 5]);
% 
% % Set title with mean and standard deviation 
% title(sprintf('RSLD ACS: %.3f \\pm %.3f, %%CV=%.2f', mu, sigma, 100*cv));
% 
% % Axis labels
% xlabel('ACS');
% ylabel('Values');
% 
% % Add a text box showing GT ACS
% text(0.05, 0.95, sprintf('GT ACS = %.2f', gt_acs), ...
%     'Units', 'normalized', 'FontSize', 12, ...
%     'Color', 'k', 'BackgroundColor', 'w', ...
%     'EdgeColor', 'k', 'HorizontalAlignment', 'left', ...
%     'VerticalAlignment', 'top');


%% BOX PLOTS ALL

% Find the maximum number of rows among the arrays
maxLength  = max([length(ACS_RFM(:)), length(ACS_SLD(:)), length(ACS_RSLD(:))]);

% Pad each array with NaNs to match the maximum number of rows
ACS_RFM_padded  = [ACS_RFM(:); NaN(maxLength - length(ACS_RFM(:)), 1)];
ACS_SLD_padded  = [ACS_SLD(:); NaN(maxLength - length(ACS_SLD(:)), 1)];
ACS_RSLD_padded = [ACS_RSLD(:); NaN(maxLength - length(ACS_RSLD(:)), 1)];

% Concatenate horizontally
acs_data = [ACS_RFM_padded, ACS_SLD_padded, ACS_RSLD_padded];

labels = {'RFM', 'SLD', 'RSLD'};

% Create the box plot
figure;
boxplot(acs_data, labels);
yline(gt_acs, 'k--', 'LineWidth', 1.5); % Dashed line for GT ACS
xlabel('ACS Method');
ylabel('\alpha [dB\cdotcm^{-1}\cdotMHz^{-1}]');
title('Box Plot of ACS Maps');
% ylim([-2 2])
grid on;

%%
keyboard
%% DENOISING

% WINDOW FILTER 1D (AVERAGE AND MEDFILT)
windowSize = 9; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;

% TV
mu_TV = 3;
iter_TV = 1000;

% SAVE FPR terms WITH FIT UPDATE Feb W3

[m, n, p] = size(RSp);
% Make full ACS all freq
z_ACS_cm = z_ACS*1E2;
df_MHz = min(diff(band));

% Initialize arrays to store FPR, linear fits, and slopes
FPR_all = zeros(length(z_ACS_cm), n, p); % Store original FPR
FPR_fit_all = zeros(length(z_ACS_cm), n, p); % Store linear fit (lin_y)

slopes_all = zeros(p, n); % Store slope values
FPR_snr_all = zeros(p, n); % Store standard deviation of FPR 
ACS_RFM = zeros(p, n);

for ff = 1:p 
    for jj = 1:n
        % Compute FPR
        FPR = log(squeeze(squeeze(RSp(:,jj,ff)))) / (4*df_MHz)*Np2dB;

        % Denoising     
        % FPR   = filter(b, a, FPR); % average filter
        FPR = medfilt1(FPR, windowSize); % median filter
        % FPR = tv_denoise(FPR, mu_TV, iter_TV); % TV-1D
        
        % Perform linear fit
        [lin_slope , lin_intercept, lin_y, R_fit] = fit_linear(z_ACS_cm, FPR, 2); 
        
        % Store FPR, linear fit, and slope
        FPR_all(:, jj, ff) = FPR;
        FPR_fit_all(:, jj, ff) = lin_y;


        slopes_all(ff, jj) = lin_slope;
        
        % Store ACS value
        ACS_RFM(ff, jj) = -lin_slope; % -8.68*m/(4*delta_f) [dB/MHz/cm]

        % std
        FPR_snr_all(ff, jj) = mean(FPR)/std(FPR);
    end
end

% WITH SLOPE PRINT

% Find freqs MHz
freq1 = 3.5; freq2 = 6; freq3 = 8.5;
idx_f1 = find(band >= freq1, 1, 'first');
idx_f2 = find(band >= freq2, 1, 'first');
idx_f3 = find(band >= freq3, 1, 'first');

% Define the range of central horizontal positions to average
pos_range = 1:n; % Selected range of horizontal positions

% Compute the mean FPR and the mean of the linear fit over the selected range
FPR_avg_f1 = mean(FPR_all(:, pos_range, idx_f1), 2); 
FPR_avg_f2 = mean(FPR_all(:, pos_range, idx_f2), 2); 
FPR_avg_f3 = mean(FPR_all(:, pos_range, idx_f3), 2); 

FPR_fit_avg_f1 = mean(FPR_fit_all(:, pos_range, idx_f1), 2);
FPR_fit_avg_f2 = mean(FPR_fit_all(:, pos_range, idx_f2), 2);
FPR_fit_avg_f3 = mean(FPR_fit_all(:, pos_range, idx_f3), 2);

% Compute the average slope over the selected positions
slope_avg_f1 = mean(slopes_all(idx_f1, pos_range), 2);
slope_avg_f2 = mean(slopes_all(idx_f2, pos_range), 2);
slope_avg_f3 = mean(slopes_all(idx_f3, pos_range), 2);

% Plot the results
figure, 
sgtitle('Frequency Power Ratio')
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 900, 700]); % [x, y, width, height] in pixels

% Plot raw FPR data
plot(z_ACS_cm, FPR_avg_f1, 'r-', 'DisplayName', sprintf('S_{%.2fMHz}/S_{%.2fMHz}', band(idx_f1), band(idx_f1-1)), 'LineWidth', 1.5), hold on
plot(z_ACS_cm, FPR_avg_f2, 'b-', 'DisplayName', sprintf('S_{%.2fMHz}/S_{%.2fMHz}', band(idx_f2), band(idx_f2-1)), 'LineWidth', 1.5), hold on
plot(z_ACS_cm, FPR_avg_f3, 'k-', 'DisplayName', sprintf('S_{%.2fMHz}/S_{%.2fMHz}', band(idx_f3), band(idx_f3-1)), 'LineWidth', 1.5), hold on

% Plot linear fits with slope values in the legend
plot(z_ACS_cm, FPR_fit_avg_f1, 'r--', 'DisplayName', sprintf('Fit: S_{%.2fMHz}/S_{%.2fMHz} (Slope: %.3f)', band(idx_f1), band(idx_f1-1), slope_avg_f1), 'LineWidth', 2), hold on
plot(z_ACS_cm, FPR_fit_avg_f2, 'b--', 'DisplayName', sprintf('Fit: S_{%.2fMHz}/S_{%.2fMHz} (Slope: %.3f)', band(idx_f2), band(idx_f2-1), slope_avg_f2), 'LineWidth', 2), hold on
plot(z_ACS_cm, FPR_fit_avg_f3, 'k--', 'DisplayName', sprintf('Fit: S_{%.2fMHz}/S_{%.2fMHz} (Slope: %.3f)', band(idx_f3), band(idx_f3-1), slope_avg_f3), 'LineWidth', 2), hold on

grid on
xlabel('Depth [cm]')
ylabel('FPR [dB/MHz]')
% ylabel('FPR [dB\cdotMHz^{-1}]')
legend('Location', 'Best')
ylim([-10 12]);
set(gca, 'FontSize', 14)

% BOX PLOT RESULTS RFM

% Compute mean and standard deviation
mu      = mean(ACS_RFM(:));
sigma   = std(ACS_RFM(:));
cv      = sigma/mu;

% Generate the box plot
figure;
boxplot(ACS_RFM(:)); 
grid on;
yline(gt_acs, 'k--', 'LineWidth', 1.5); % Dashed line for GT ACS
ylim([-5 5]);

% Set title with mean and standard deviation 
title(sprintf('RFM ACS: %.3f \\pm %.3f, %%CV=%.2f', mu, sigma, 100*cv));

% Axis labels
xlabel('ACS');
ylabel('Values');

% Add a text box showing GT ACS
text(0.05, 0.95, sprintf('GT ACS = %.2f', gt_acs), ...
    'Units', 'normalized', 'FontSize', 12, ...
    'Color', 'k', 'BackgroundColor', 'w', ...
    'EdgeColor', 'k', 'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'top');


%%
pos_range = 1:n;
pos_range = floor(n/2);
% pos_range = 
% Find freqs MHz
freq1 = 3.5; freq2 = 6; freq3 = 8.5;
idx_f1 = find(band >= freq1, 1, 'first');
idx_f2 = find(band >= freq2, 1, 'first');
idx_f3 = find(band >= freq3, 1, 'first');

FPR_avg_f3 = mean(FPR_all(:, pos_range, idx_f3), 2); 

figure, 
plot(z_ACS_cm, FPR_avg_f3, 'k-', 'DisplayName', sprintf('S_{%.2fMHz}/S_{%.2fMHz}', band(idx_f3), band(idx_f3-1)), 'LineWidth', 1.5), hold on
grid on
xlabel('Depth [cm]')
ylabel('FPR [dB/MHz]')
% ylabel('FPR [dB\cdotMHz^{-1}]')
legend('Location', 'Best')
ylim([-10 5]);
set(gca, 'FontSize', 14)

figure,
list_mu_TV = logspace(log10(0.001), log10(7), 15);
iter_TV = 100;

for uu = 1:length(list_mu_TV)
    mu_TV = list_mu_TV(uu);
    FPR_avg_f3_den = tv_denoise(FPR, mu_TV, iter_TV); % TV-1D

    plot(z_ACS_cm, FPR_avg_f3_den, '-', 'DisplayName', sprintf('u = %.3g', mu_TV), 'LineWidth', 1.5), hold on
    hold on
    ylim([-19 5]);
    pause(0.5)
end
grid on
xlabel('Depth [cm]')
ylabel('FPR [dB/MHz]')
% ylabel('FPR [dB\cdotMHz^{-1}]')
legend('Location', 'Best')

set(gca, 'FontSize', 14)

%%
% Define log-spaced values for mu_TV
iter_TV = 100;
mu_TV_list = logspace(log10(0.05), log10(7), 15); 

% Initialize arrays to store ACS results for each mu_TV
ACS_RFM_all = cell(length(mu_TV_list), 1);

% Loop over different mu_TV values
for k = 1:length(mu_TV_list)
    mu_TV = mu_TV_list(k);
    
    % Initialize arrays to store FPR, linear fits, and slopes
    FPR_all = zeros(length(z_ACS_cm), n, p);
    FPR_fit_all = zeros(length(z_ACS_cm), n, p);
    slopes_all = zeros(n, p);

    for ff = 1:p 
        for jj = 1:n
            % Compute FPR
            FPR = log(squeeze(RSp(:,jj,ff))) / (4*df_MHz) * Np2dB;

            % Denoising using current mu_TV
            % FPR = tv_denoise(FPR, mu_TV, iter_TV);
            FPR = tv_denoise(FPR, mu_TV*fact_freq(ff), iter_TV); % TV-1D

            % Perform linear fit
            [lin_slope, lin_intercept, lin_y, R_fit] = fit_linear(z_ACS_cm, FPR, 2); 

            % Store results
            FPR_all(:, jj, ff) = FPR;
            FPR_fit_all(:, jj, ff) = lin_y;
            slopes_all(jj, ff) = lin_slope;
            
            % Store ACS value
            ACS_RFM(ff, jj) = -lin_slope; 
        end
    end
    
    % Store ACS results
    ACS_RFM_all{k} = ACS_RFM(:);
end

%% BOX PLOTS FOR EACH mu_TV
figure;
for k = 1:length(mu_TV_list)
    % Extract ACS values for current mu_TV
    acs_values = ACS_RFM_all{k};
    
    % Compute mean, standard deviation, and coefficient of variation
    mu = mean(acs_values);
    sigma = std(acs_values);
    cv = (sigma / mu) * 100;

    % Create subplot for each mu_TV
    subplot(2, 5, k); % Arrange in 2 rows, 5 columns
    boxplot(acs_values);
    grid on;
    yline(gt_acs, 'k--', 'LineWidth', 1.5); % Dashed line for GT ACS
    ylim([-5 5]);

    % Set title with mu_TV value
    title(sprintf('\\mu_{TV} = %.3f', mu_TV_list(k)));
    
    % Axis labels
    xlabel('ACS');
    ylabel('Values');

    % Add text annotation with mean, std, and CV
    text(0.7, 0.1, sprintf('%.3f \\pm %.3f, CV=%.2f%%', mu, sigma, cv), ...
        'Units', 'normalized', 'FontSize', 10, ...
        'Color', 'k', 'BackgroundColor', 'w', ...
        'EdgeColor', 'k', 'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'bottom');
end

% Global title for the entire figure
sgtitle('Effect of Different \mu_{TV} on ACS Estimation');


%% PLOT SNR vs SHIFT FREQ

figure, 
subplot(121)
imagesc(x_ACS*1e3, band, FPR_snr_all)
% axis("image")
xlabel('Lateral [mm]'), ylabel('Freq [MHz]')
title('\bfSNR FPR');

pos_range = 1:n;
% Compute the average slope over the selected positions
snr_FPR_pos_range = mean(FPR_snr_all(:, pos_range), 2);

subplot(122)
plot(band, snr_FPR_pos_range)
grid on;
xlabel('Freq [MHz]')
ylabel('SNR')
title('\bfSNR FPR');

%%
pos_range = 1:n;
% Compute the average slope over the selected positions
snr_FPR_pos_range = mean(FPR_snr_all(:, pos_range), 2);
snr_FPR_pos_range = rescale(snr_FPR_pos_range, 0.1, 1);
figure, 
plot(band, snr_FPR_pos_range)
grid on;
xlabel('bfFreq [MHz]')
ylabel('\bfSNR FPR ')
title('\bfSNR FPR')

figure, 
plot(band, pow2db(snr_FPR_pos_range))
grid on;
xlabel('bfFreq [MHz]')
ylabel('\bfSNR [dB]')
title('\bfNorm SNR FPR')

figure, 
semilogy(band, snr_FPR_pos_range)
grid on;
xlabel('bfFreq [MHz]')
ylabel('\bfSNR FPR [dB]')
title('\bfSNR FPR')

%% DENOISING

% WINDOW FILTER 1D (AVERAGE AND MEDFILT)
windowSize = 9; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;

% TV
mu_TV = 10;
fact_freq = snr_FPR_pos_range;
iter_TV = 1000;

% SAVE FPR terms WITH FIT UPDATE Feb W3

[m, n, p] = size(RSp);
% Make full ACS all freq
z_ACS_cm = z_ACS*1E2;
df_MHz = min(diff(band));

% Initialize arrays to store FPR, linear fits, and slopes
FPR_all = zeros(length(z_ACS_cm), n, p); % Store original FPR
FPR_fit_all = zeros(length(z_ACS_cm), n, p); % Store linear fit (lin_y)

slopes_all = zeros(p, n); % Store slope values
FPR_snr_all = zeros(p, n); % Store standard deviation of FPR 
ACS_RFM = zeros(p, n);

for ff = 1:p 
    for jj = 1:n
        % Compute FPR
        FPR = log(squeeze(squeeze(RSp(:,jj,ff)))) / (4*df_MHz)*Np2dB;

        % Denoising     
        % FPR   = filter(b, a, FPR); % average filter
        % FPR = medfilt1(FPR, windowSize); % median filter
        % FPR = tv_denoise(FPR, mu_TV, iter_TV); % TV-1D
        FPR = tv_denoise(FPR, mu_TV*fact_freq(ff), iter_TV); % TV-1D
        
        % Perform linear fit
        [lin_slope , lin_intercept, lin_y, R_fit] = fit_linear(z_ACS_cm, FPR, 2); 
        
        % Store FPR, linear fit, and slope
        FPR_all(:, jj, ff) = FPR;
        FPR_fit_all(:, jj, ff) = lin_y;

        slopes_all(ff, jj) = lin_slope;
        
        % Store ACS value
        ACS_RFM(ff, jj) = -lin_slope; % -8.68*m/(4*delta_f) [dB/MHz/cm]

        % std
        FPR_snr_all(ff, jj) = mean(FPR)/std(FPR);
    end
end

% WITH SLOPE PRINT

% Find freqs MHz
freq1 = 3.5; freq2 = 6; freq3 = 8.5;
idx_f1 = find(band >= freq1, 1, 'first');
idx_f2 = find(band >= freq2, 1, 'first');
idx_f3 = find(band >= freq3, 1, 'first');

% Define the range of central horizontal positions to average
pos_range = 1:n; % Selected range of horizontal positions

% Compute the mean FPR and the mean of the linear fit over the selected range
FPR_avg_f1 = mean(FPR_all(:, pos_range, idx_f1), 2); 
FPR_avg_f2 = mean(FPR_all(:, pos_range, idx_f2), 2); 
FPR_avg_f3 = mean(FPR_all(:, pos_range, idx_f3), 2); 

FPR_fit_avg_f1 = mean(FPR_fit_all(:, pos_range, idx_f1), 2);
FPR_fit_avg_f2 = mean(FPR_fit_all(:, pos_range, idx_f2), 2);
FPR_fit_avg_f3 = mean(FPR_fit_all(:, pos_range, idx_f3), 2);

% Compute the average slope over the selected positions
slope_avg_f1 = mean(slopes_all(idx_f1, pos_range), 2);
slope_avg_f2 = mean(slopes_all(idx_f2, pos_range), 2);
slope_avg_f3 = mean(slopes_all(idx_f3, pos_range), 2);

% Plot the results
figure, 
sgtitle('Frequency Power Ratio')
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 900, 700]); % [x, y, width, height] in pixels

% Plot raw FPR data
plot(z_ACS_cm, FPR_avg_f1, 'r-', 'DisplayName', sprintf('S_{%.2fMHz}/S_{%.2fMHz}', band(idx_f1), band(idx_f1-1)), 'LineWidth', 1.5), hold on
plot(z_ACS_cm, FPR_avg_f2, 'b-', 'DisplayName', sprintf('S_{%.2fMHz}/S_{%.2fMHz}', band(idx_f2), band(idx_f2-1)), 'LineWidth', 1.5), hold on
plot(z_ACS_cm, FPR_avg_f3, 'k-', 'DisplayName', sprintf('S_{%.2fMHz}/S_{%.2fMHz}', band(idx_f3), band(idx_f3-1)), 'LineWidth', 1.5), hold on

% Plot linear fits with slope values in the legend
plot(z_ACS_cm, FPR_fit_avg_f1, 'r--', 'DisplayName', sprintf('Fit: S_{%.2fMHz}/S_{%.2fMHz} (Slope: %.3f)', band(idx_f1), band(idx_f1-1), slope_avg_f1), 'LineWidth', 2), hold on
plot(z_ACS_cm, FPR_fit_avg_f2, 'b--', 'DisplayName', sprintf('Fit: S_{%.2fMHz}/S_{%.2fMHz} (Slope: %.3f)', band(idx_f2), band(idx_f2-1), slope_avg_f2), 'LineWidth', 2), hold on
plot(z_ACS_cm, FPR_fit_avg_f3, 'k--', 'DisplayName', sprintf('Fit: S_{%.2fMHz}/S_{%.2fMHz} (Slope: %.3f)', band(idx_f3), band(idx_f3-1), slope_avg_f3), 'LineWidth', 2), hold on

grid on
xlabel('Depth [cm]')
ylabel('FPR [dB/MHz]')
% ylabel('FPR [dB\cdotMHz^{-1}]')
legend('Location', 'Best')
ylim([-10 12]);
set(gca, 'FontSize', 14)

% BOX PLOT RESULTS RFM

% Compute mean and standard deviation
mu      = mean(ACS_RFM(:));
sigma   = std(ACS_RFM(:));
cv      = sigma/mu;

% Generate the box plot
figure;
boxplot(ACS_RFM(:)); 
grid on;
yline(gt_acs, 'k--', 'LineWidth', 1.5); % Dashed line for GT ACS
ylim([-5 5]);

% Set title with mean and standard deviation 
title(sprintf('RFM ACS: %.3f \\pm %.3f, %%CV=%.2f', mu, sigma, 100*cv));

% Axis labels
xlabel('ACS');
ylabel('Values');

% Add a text box showing GT ACS
text(0.05, 0.95, sprintf('GT ACS = %.2f', gt_acs), ...
    'Units', 'normalized', 'FontSize', 12, ...
    'Color', 'k', 'BackgroundColor', 'w', ...
    'EdgeColor', 'k', 'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'top');

%% WITH ERRORS

%% FULL RANGE (I believe with error)
a_local = zeros(m, n); % Preallocate local attenuation matrix (depth x lateral)

tic;
for jj = 1:n  % Loop over lateral positions (x_j)
    for ii = 1:m  % Loop over depth positions (z_k)
        
        y_vec = []; % Initialize y vector for this location
        X_mat = []; % Initialize X matrix for this location
        
        for r = 1:m_r  % Loop over reference depths
            for i = 2:p  % Loop over frequency bins
                % Compute y = log(RSnorm) at this depth & lateral position
                y = log(RSp_k(ii, jj, i)) - log(RSp_r(r, jj, i));
                
                % Define X = -4 * (fi - fi-1) * (zk - zr)
                X = -4 * (band(i) - band(i-1)) * (z_ACS(ii) - z_ACS_r(r));

                % Store values for least squares regression
                y_vec = [y_vec; y(:)];
                X_mat = [X_mat; X];
            end
        end

        % Solve for local attenuation a(z_k, x_j) using least squares
        if ~isempty(y_vec)
            a_local(ii, jj) = (X_mat' * X_mat) \ (X_mat' * y_vec);
        end
    end
end
t = toc;
fprintf('Elapsed time %.2f \n', t);


[m_a, s_a, cv_a] = deal(calc2dStats{1}(a_local), calc2dStats{2}(a_local), calc2dStats{3}(a_local));

caxis_acs = [0 1.1];
fontSize = 14;

figure;
subplot(121)
imagesc(x_ACS * 1e3, z_ACS * 1e3, a_local, caxis_acs); % Convert to mm
axis("image")
colorbar; colormap("turbo")
xlabel('Lateral [mm]');
ylabel('Depth [mm]');
hb2=colorbar; ylabel(hb2,'dB\cdotcm^{-1}\cdotMHz^{-1}')
% title('Local Attenuation Coefficient');
title(sprintf('Local AC (GT= %.2f)\n%.3f $\\pm$ %.3f,  \\%%CV = %.2f', ...
               alpha_sam, m_a, s_a, cv_a), ...
      'Interpreter', 'latex');

set(gca,'fontsize',fontSize)