% ====================================================================== %
% Script for clinical data.
% Created on June 19th, 2024
% Modification for RFM method EMZ
% ====================================================================== %

% reinit
baseDir = 'D:\emirandaz\qus\data\liver\healthy';

sampleDir = fullfile(baseDir,'samples');
% refsDir = fullfile(baseDir,'refs','joined');
resultsDir = 'D:\emirandaz\qus\rfm\tnv\healthyLiver';
figsDir = 'D:\emirandaz\qus\rfm\tnv\healthyLiver';

if ~exist(resultsDir) mkdir(resultsDir); end
if ~exist(figsDir) mkdir(figsDir); end


% acqDir = dir(fullfile(sampleDir,'016-03.mat')); %65*ma
% acqDir = dir(fullfile(sampleDir,'007-05.mat')); %*ma
% acqDir = dir(fullfile(sampleDir,'014-01.mat')); %*ma
% acqDir = dir(fullfile(sampleDir,'016-06.mat')); %emz
acqDir = dir(fullfile(sampleDir,'020-05.mat')); %dv

%% FIRST UTILS
manualroi = false;
denTNVRFM = false;

Np2dB       = 20*log10(exp(1));
dB2Np       = 1/Np2dB;
range_bmode = [-60 0];

mean2d = @(x) mean(x(:));
std2d = @(x) std(x(:));
cv2d = @(x) 100*std(x(:))/mean(x(:));

calc2dStats = {@(x) mean(x(:)), @(x) std(x(:)), @(x) 100 * std(x(:)) / mean(x(:))};
%%
for iFile = 1:length(acqDir)

% samName = sampleFiles(iFile).name(1:end-4);
% fileName = fullfile(pathData, samName);

% iAcq = list_acqDir (iList);
samName = acqDir(iFile).name(1:end-4);
fileName = fullfile(sampleDir, samName);

%% SPECTRAL PARAMETERS
pars.bw          = [1.2 3.9]; % [MHz]
pars.overlap     = 0.8;
pars.blocksize   = 10; % wavelengths
pars.blocklines  = 8;
pars.saran_layer = false;
pars.ratio_zx    = 1;
pars.window_type = 3; %  (1) Hanning, (2) Tuckey, (3) Hamming, (4) Tchebychev

blocksize_wv_r = 18;

fprintf('ID N° %d : %s\n', iFile, samName(1:end-4));

%% Loading file and variables
load(fullfile(sampleDir,samName+".mat"));
fprintf("Loading sample %s \n", samName)
RcvData = cell2mat(RcvData);
n_frame = size(RcvData,3); % Frame selector
RcvData = RcvData(:, :, n_frame); % Select frame of RF Data

% Additional variables
central_freq = Receive(1).demodFrequency*1e6; % Central frequency of pulse
fs = Receive(1).decimSampleRate*1e6; % According to "NS200BW" Acquisition Mode
n_pulses = P.numRays; % number of pulses
n_elements = Trans.numelements; % number of transducer elements
num_samples = Receive(1).endSample - Receive(1).startSample +1; % samples per channel
sound_speed = Resource.Parameters.speedOfSound; % [m/s]
wvl = sound_speed/central_freq; % [m] wavelength in meters
scalemm2wvl = 1/wvl;

% Initialize variables
rf_channel = zeros(num_samples , n_elements, n_pulses);
rx_apods = zeros(1, n_elements, n_pulses);
rf = zeros(num_samples, n_pulses);

% Organize data
for n = 1:n_pulses % Iterate through pulses
    rf_channel(:, :, n) = RcvData(Receive(n).startSample:Receive(n).endSample, :);
end


% Acommodate to time delays and rf signals
focus = 20/1000;
t = (0:(num_samples-1))/fs; % [sec.] time domain 0:T:(N_sample-1)*T
[rx_delays] = getRXDelays(Trans, t, n_elements, n_pulses, sound_speed, wvl);


% Dynamic Aperture
f_num = 3;
z = sound_speed*t/2;
elem_pitch = Trans.spacingMm*1e-3;
maxAprSz = 32;
dyn_aperture = zeros(length(z), n_elements, n_pulses);
for n = 1:n_pulses
    for z_i = 1:length(z)
        a = z(z_i)/(2*f_num);
        hlfAprSz = floor(a / elem_pitch);
        if (hlfAprSz > maxAprSz/2)
            hlfAprSz = floor(maxAprSz / 2);
        end
        a_i = -hlfAprSz: hlfAprSz;    % aperture indices
        fulAprSz = 2*hlfAprSz + 1;
        aper_center = n;
        aper = aper_center + a_i;
        aper = aper(aper>=1);
        aper = aper(aper<=128);
        dyn_aperture(z_i, aper, n) = 1;
    end
end


%% Beamforming

% Delay-and-sum
for n = 1:n_pulses
    % Delaying
    for e = 1:n_elements
        % rx_apods(e, n);
        rf_channel(:, e, n) = ...
            interp1(t, rf_channel(:, e, n), rx_delays(:,e, n), 'linear', 0);
    end
    rf_channel(:, :, n) = rf_channel(:, :, n) .* dyn_aperture(:, :, n);
    % Summing
    rf(:, n) = sum(rf_channel(:, :, n), 2);
end

%% B-Mode and coordinates
b_mode = 20*log10(abs(hilbert(rf)));
b_mode = b_mode-max(b_mode(:));

param = getparam('C5-2v');
siz = size(rf);
z = sound_speed*t/2;
zmax = z(end);
R = param.radius;
p = param.pitch;
N = param.Nelements;
L = 2*R*sin(asin(p/2/R)*(N-1)); % chord length
d = sqrt(R^2-L^2/4); % apothem
z0 = -d;

th = -(linspace(atan2(L/2,d),atan2(-L/2,d),siz(2)))*180/pi;
r = linspace(R+p,-z0+zmax,siz(1));

% To Polar Coordinates
[xPolar,zPolar, z0Polar] = impolgrid(size(b_mode), z(end),param);

%%
BmodeFull = db(hilbert(rf));
BmodeFull = BmodeFull - max(BmodeFull(:));
caption = samName;
% figure,
% set(gcf, 'Units', 'pixels', 'Position', [100, 100, 600, 600]); % [x, y, width, height]
% pcolor(xPolar*1e3, zPolar*1e3, BmodeFull);
% shading interp
% cbar = colorbar;
% ylabel(cbar, 'dB', 'FontSize', 12);
% clim([-70 0]);
% colormap gray               
% ylabel('Depth [mm]', 'FontSize', 14);
% xlabel('Lateral [mm]', 'FontSize', 14);
% axis equal ij tight
% title(caption, 'FontSize', 14);
% set(gca, 'Color', 'k');  % axes background


%% Selecting ROI

xFull = th; % [deg]
r0 = r(1);
zFull = (r-r0)*1e2; % [cm]

SAM.z = zFull*1e-2; % to [m]
SAM.x = xFull;
SAM.fs = fs;
SAM.rf = rf;
SAM.bMode = BmodeFull;
bmode_sam   = SAM.bMode;

%%
% figure('Units','centimeters', 'Position',[5 5 15 15]),
% imagesc(xFull,SAM.z*1e3,BmodeFull,range_bmode);
% % ylim([0 10])
% colormap gray; clim(range_bmode);
% hb2=colorbar; ylabel(hb2,'dB')
% xlabel('Lateral [º]'), ylabel('Depth [mm]'); 
% % title(caption)
% 
% confirmation = '';
% while ~strcmp(confirmation,'Yes')
%     rect = getrect;
%     confirmation = questdlg('Sure?');
%     if strcmp(confirmation,'Cancel')
%         disp(rect)
%         break
%     end
% end
% close,
% 
% pars.x_roi     = [rect(1), rect(1)+rect(3)];      % [º]
% pars.z_roi     = [rect(2), rect(2)+rect(4)]*1E-3; % [m]


pars.z_roi = [0.0415    0.0877]; % 0.0462
pars.x_roi = [1.0492   30.8074]; % 29.75

% NEW
pars.x_roi = [3.0992   28.7574]; % 25.75
pars.z_roi = [0.0435    0.0857]; % 0.0442
%%

% fprintf('ID N° %d : %s\n', iFile, samName);
% 
% if ( strcmp(samName, '020-04') )
%     pars.x_roi = [-15.3560 12.8762];
%     pars.z_roi = [0.0398 0.0776];
% 
% end
%%

% RFM V1 
% Reading experiment settings parameters
bw              = pars.bw;
overlap         = pars.overlap;
blocksize_wv    = pars.blocksize;
blocklines      = pars.blocklines;
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

% Lateral samples (LINEAR)
% wx = round(blocksize_wv*lambda*(1-overlap)/dx);  % Between windows  
% nx = round(blocksize_wv*lambda/ dx);

% Lateral samples (CURVE)
wx = round(blocklines*(1-overlap));  % Between windows
nx = blocklines;                 % Window size

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

band     = axis_f(ind_f)*1E-6; % [MHz]
bandFull = axis_f*1E-6; % [MHz]

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
RSp_k(:,:, 2:end) = Sp_k(:,:, 2:end) ./ Sp_k(:,:, 1:end-1); 
% clear Sp_k
% For the first slice, keep the ratio the same as the first slice
RSp_k(:,:, 1) = RSp_k(:,:, 2); % Assuming the first slice ratios are 1 as there is no "i-1"

RSp_k = log(RSp_k); % @@
% Reference Depth

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
RSp_r(:,:, 2:end) = Sp_r(:,:, 2:end) ./ Sp_r(:,:, 1:end-1); 
% clear Sp_r
% For the first slice, keep the ratio the same as the first slice
RSp_r(:,:, 1) = RSp_r(:,:, 2); % Assuming the first slice ratios are 1 as there is no "i-1"

RSp_r = log(RSp_r); % @@
% USE THIS FROM APRIL AND NOW ON
%%%%%%%%%%%%%%%%%%%% FAST WAY %%%%%%%%%%%%%%%%%%%%

% UFR strategy

bw_ufr = pars.bw;
freqL = bw_ufr(1); freqH = bw_ufr(2);
range = bandFull >= freqL & bandFull <= freqH;

band_ufr    = bandFull(range);
p_ufr       = length(band_ufr);

RSp_k_ufr   = RSp_k(:,:,range);
RSp_r_ufr   = RSp_r(:,:,range);

%
% DENOISING TNV RSP
mu          = 1.3;
tau         = 0.0100;
maxIter     = 1000;
stableIter  = 50;
% tol         = 0.5e-4; % tolerance error
tol         = 1e-3; % tolerance error
RSp_k_ufr_vec = reshape(RSp_k_ufr, [], size(RSp_k_ufr, 3));   
mux_RSp = 1./(abs(mean(RSp_k_ufr_vec, 1, 'omitnan')) ./ std(RSp_k_ufr_vec, 0, 1, 'omitnan') + 1E-5);
% weigthChannels = rescale(mux_RSp, 1, 10);
% weigthChannels = mux_RSp;
weigthChannels = ones(1, p_ufr);
 
% gamma = 1;
% vec_gamma = 1 + (10 - 1) * (1 - linspace(0, 1, p_ufr).^gamma);
% vec_log = logspace(log10(10), log10(1), p_ufr);
% figure, plot(vec_gamma)

% weigthChannels = vec_gamma;

[RSp_k_ufr_opt, cost, error, fid, reg] = pdo_den_wtnv(RSp_k_ufr, mu, tau, maxIter, tol, stableIter, weigthChannels);

[RSp_r_ufr_opt, cost, error, fid, reg] = pdo_den_wtnv(RSp_r_ufr, mu, tau, maxIter, tol, stableIter, weigthChannels);

% Plot settings
% rows = 3; cols = 5;
% total_plots = rows * cols;
% l_Rps = round(linspace(1, length(band_ufr), total_plots));
% 
% % ORIGINAL
% figure;
% t = tiledlayout(rows, cols, 'Padding', 'compact', 'TileSpacing', 'compact');
% for i = 1:total_plots
%     slice_idx = l_Rps(i);
%     nexttile;
%     imagesc(RSp_k_ufr(:, :, slice_idx));
%     axis off;
%     title(['RSp ', num2str(round(band_ufr(slice_idx),2)), 'MHz']);
% end
% 
% % TNV
% figure;
% t = tiledlayout(rows, cols, 'Padding', 'compact', 'TileSpacing', 'compact');
% for i = 1:total_plots
%     slice_idx = l_Rps(i);
%     nexttile;
%     imagesc(RSp_k_ufr_opt(:, :, slice_idx));
%     axis off;
%     title(['TNV RSp ', num2str(round(band_ufr(slice_idx),2)), 'MHz']);
% end
% RFM METHOD

% Convert depth values to cm
z_ACS_cm = z_ACS * 1e2;      % Convert from meters to cm
z_ACS_r_cm = z_ACS_r * 1e2;  % Convert reference depths to cm

% Delta MHz 
df_MHz = band_ufr(2) - band_ufr(1);

% Preallocate cell arrays for storing x_temp and y_temp
x_temp_all = cell(m, n);
y_temp_all = cell(m, n);
y_temp_all2 = cell(m, n);

tic;
for jj = 1:n  % Loop over lateral positions (x_j)
    for ii = 1:m  % Loop over depth positions (z_k)

        % Temporary storage for this location
        y_temp = nan(m_r, p_ufr);  % (Reference depth, Frequency) ** m_r
        y_temp2 = nan(m_r, p_ufr);  % (Reference depth, Frequency) ** m_r
        x_temp = nan(m_r, p_ufr);  % (Reference depth, Frequency) ** m_r

        for r = 1:m_r  % Loop over reference depths
            % if (ii==1 && r==1)
            y_col = squeeze( ( (RSp_k_ufr(ii, jj, :)) - (RSp_r_ufr(r, jj, :)) ) /(4*df_MHz) *Np2dB ); % p_ufr x 1
            y_col2 = squeeze( ( (RSp_k_ufr_opt(ii, jj, :)) - (RSp_r_ufr_opt(r, jj, :)) ) /(4*df_MHz) *Np2dB ); % p_ufr x 1
            
            X = z_ACS_cm(ii) - z_ACS_r_cm(r);

            y_temp(r, :) = y_col; %**
            y_temp2(r, :) = y_col2;
            x_temp(r, :) = X; % **

        end
        x_temp_all{ii, jj} = x_temp;
        y_temp_all{ii, jj} = y_temp;
        y_temp_all2{ii, jj} = y_temp2;

    end
end
t = toc;
fprintf('Loop way Elapsed time %.2f \n', t);
%%%%%%%%%%%%%%%%%%%% FAST WAY %%%%%%%%%%%%%%%%%%%%

% NOW ESTIMATION CGS RFM

a_rfm = zeros(m, n); 
a_rfm2 = zeros(m, n); 
for jj = 1:n  % Loop over lateral positions (x_j)
    for ii = 1:m  % Loop over depth positions (z_k)

        x_temp = x_temp_all{ii, jj};
        y_temp = y_temp_all{ii, jj};
        y_temp2 = y_temp_all2{ii, jj};

        X_vec = x_temp(:);
        y_vec = y_temp(:);
        y_vec2 = y_temp2(:);

        % % Option 1
        % a_rfm(ii, jj) = -(X_vec' * X_vec) \ (X_vec' * y_vec);
        % a_rfm(ii, jj) = abs((X_vec' * X_vec) \ (X_vec' * y_vec));
        % 
        % a_rfm2(ii, jj) = -(X_vec' * X_vec) \ (X_vec' * y_vec2);
        % a_rfm2(ii, jj) = abs((X_vec' * X_vec) \ (X_vec' * y_vec2));


        % % Option 2
        A2_full    = kron(speye(p_ufr), x_temp(:, 1));
        A2t_full   = A2_full';

        slope_full = -(cgs(A2t_full*A2_full, A2t_full*y_vec));
        slope2_full = -(cgs(A2t_full*A2_full, A2t_full*y_vec2));

        a_rfm(ii, jj) = mean(slope_full);
        a_rfm2(ii, jj) = mean(slope2_full);

        % a_rfm(ii, jj) = median(slope_full);
        % a_rfm2(ii, jj) = median(slope2_full);

        % % Option 3 (INTERCEPT)

        % A31_full   = kron( speye(p_ufr), x_temp(:, 1) );
        % A32_full   = kron( speye(p_ufr), ones(size( x_temp(:, 1))) );
        % A3_full    = horzcat(A31_full, A32_full);
        % A3t_full   = A3_full';
        % 
        % x_opt = -cgs(A3t_full*A3_full, A3t_full*y_vec);
        % slope_full = x_opt(1:end/2);
        % 
        % x2_opt = -cgs(A3t_full*A3_full, A3t_full*y_vec2);
        % slope2_full = x2_opt(1:end/2);
        % 
        % a_rfm(ii, jj) = mean(slope_full);
        % a_rfm2(ii, jj) = mean(slope2_full);
        % 
        % % a_rfm(ii, jj) = median(slope_full);
        % % a_rfm2(ii, jj) = median(slope2_full);


        % % Option 4 (LINEAR FIT)
        % slope_fit_fx = zeros(1, p_ufr);
        % slope_fit_fx2 = zeros(1, p_ufr);
        % for pp = 1:p_ufr
        %     [slope_fx, ~, ~, ~] = fit_linear(x_temp(:, pp), y_temp(:, pp), 2);
        %     [slope_fx2, ~, ~, ~] = fit_linear(x_temp(:, pp), y_temp2(:, pp), 2);
        %     slope_fit_fx(1,pp) =    -(slope_fx) ;
        %     slope_fit_fx2(1,pp) =   -(slope_fx2) ;
        % end
        % a_rfm(ii, jj) = mean(slope_fit_fx);
        % a_rfm2(ii, jj) = mean(slope_fit_fx2);
        % 
        % % a_rfm(ii, jj) = median(slope_fit_fx);
        % % a_rfm2(ii, jj) = median(slope_fit_fx2);


    end
end
%
%% RFM
[m_a, s_a, cv_a]    = deal(calc2dStats{1}(a_rfm), calc2dStats{2}(a_rfm), calc2dStats{3}(a_rfm));
[m_a2, s_a2, cv_a2] = deal(calc2dStats{1}(a_rfm2), calc2dStats{2}(a_rfm2), calc2dStats{3}(a_rfm2));

range_acs = [0 1.5];
fact_transp = 0.7;
fontSize = 15;
% 
% figure, 
% imagesc(a_rfm, range_acs); % Convert to mm
% % % axis("image")
% colorbar; colormap("turbo")

% figure('Units', 'pixels', 'Position', [100, 100, 800, 600]); % [x, y, width, height]
% 
% subplot(121)
% imagesc(x_ACS * 1e3, z_ACS * 1e3, a_rfm, range_acs); % Convert to mm
% % axis("image")
% colorbar; colormap("turbo")
% xlabel('Lateral [mm]');
% ylabel('Depth [mm]');
% hb2=colorbar; ylabel(hb2,'dB\cdotcm^{-1}\cdotMHz^{-1}')
% title(sprintf('RFM %.3f$\\pm$%.3f, \\%%CV=%.2f', m_a, s_a, cv_a), ...
%       'Interpreter', 'latex');
% set(gca,'fontsize',fontSize)
% 
% subplot(122)
% imagesc(x_ACS * 1e3, z_ACS * 1e3, a_rfm2, range_acs); % Convert to mm
% % axis("image")
% colorbar; colormap("turbo")
% xlabel('Lateral [mm]');
% ylabel('Depth [mm]');
% hb2=colorbar; ylabel(hb2,'dB\cdotcm^{-1}\cdotMHz^{-1}')
% title(sprintf('TNV-RFM %.3f$\\pm$%.3f, \\%%CV=%.2f', m_a2, s_a2, cv_a2), ...
%       'Interpreter', 'latex');
% set(gca,'fontsize',fontSize)

%

xFull = SAM.x;
zFull = SAM.z;
[X,Z] = meshgrid(xFull, zFull);

roi = X >= x_ACS(1) & X <= x_ACS(end) & Z >= z_ACS(1) & Z <= z_ACS(end);

figure('Units', 'pixels', 'Position', [100, 100, 1200, 600]); % [x, y, width, height]

tiledlayout(1, 3, 'TileSpacing','compact', 'Padding','compact')

t1 = nexttile();
imagesc(xFull,zFull*1e3, bmode_sam, range_bmode); % axis image;
% title('B-mode')
title(['Bmode k', mat2str(round(pars.blocksize,2)), ...
    'r', mat2str(round(blocksize_wv_r,2)), ...
    'L', mat2str(round(pars.blocklines,2))])
% ylim(yLimits)
hold on
contour(xFull,zFull,roi,1,'w--')
hold off
xlabel('Lateral [°]')
ylabel('Axial [mm]')
hBm = colorbar('Ticks',-60:20:0);
hBm.Label.String = 'dB';
% hBm.Location = 'westoutside';
set(gca,'fontsize',fontSize)

nexttile,
[~,hB,hColor] = imOverlayInterp(bmode_sam,a_rfm,range_bmode,range_acs,fact_transp,...
    x_ACS,z_ACS*1e3,roi,xFull,zFull*1e3);
title(sprintf('RFM %.3f$\\pm$%.3f, \\%%CV=%.2f', m_a, s_a, cv_a), ...
      'Interpreter', 'latex');
% title(sprintf('RFM: %.3f±%.3f, %%CV=%.2f', m_a, s_a, cv_a));
% subtitle(['ACS = ',num2str(m_a,2),' dB/cm/MHz'])
% colorbar off
axis normal
% ylim(yLimits)
hold on
contour(xFull,zFull*1e3,roi,1,'w--')
hold off
% axis off
colorbar off
xlabel('Lateral [°]')
set(gca,'fontsize',fontSize)

nexttile,
[~,hB,hColor] = imOverlayInterp(bmode_sam,a_rfm2,range_bmode,range_acs,fact_transp,...
    x_ACS,z_ACS*1e3,roi,xFull,zFull*1e3);
title(sprintf('TNV-RFM %.3f$\\pm$%.3f, \\%%CV=%.2f', m_a2, s_a2, cv_a2), ...
      'Interpreter', 'latex');
% title(sprintf('TNV-RFM: %.3f±%.3f, %%CV=%.2f', m_a2, s_a2, cv_a2));
% subtitle(['ACS = ',num2str(m_a,2),' dB/cm/MHz'])
% colorbar off
axis normal
% ylim(yLimits)
hold on
contour(xFull,zFull*1e3,roi,1,'w--')
hold off
% axis off
xlabel('Lateral [°]')
hColor.Label.String = 'dB\cdotcm\cdotMHz^{-1}';
% hColor.Label.String = 'dB/cm/MHz';
colormap(t1,'gray')
set(gca,'fontsize',fontSize)

% SAVE FIG (RECT)
% exportgraphics(gcf,fullfile(figsDir,samName(1:end-4)+"_rect.png"), ...
%     'Resolution','300')

% keyboard

%
%% MULTI PREV TEST

dirMulti = ['D:\emirandaz\qus\rfm\IUS2025\HL\multil_',num2str(pars.blocklines),'vIUS'];
if ~exist(dirMulti) mkdir(dirMulti); end

% ['Bmode k', mat2str(round(pars.blocksize,2)), ...
%     'r', mat2str(round(blocksize_wv_r,2)), ...
%     'L', mat2str(round(pars.blocklines,2))]
save_all_figures_to_directory(dirMulti,char("M"+pars.blocksize+"fig"))
close all

%%
 dirOut = 'D:\emirandaz\qus\rfm\IUS2025\abstract\hl';
if ~exist(dirOut) mkdir(dirOut); end

%% Plot RFM

% [xPolar,zPolar, z0Polar] = impolgrid(size(b_mode), z(end),param);

[X,Z] = meshgrid(xFull, zFull);
% r0      = SAM.r0;
r0 = 0.0501;
fact_transp = 0.45;

fontSize = 26;
roi = X >= x_ACS(3) & X <= x_ACS(end) & Z >= z_ACS(1) & Z <= z_ACS(end);

[TH_acs,R_acs] = meshgrid(-x_ACS*pi/180 + pi/2,z_ACS + r0);
[xPolarACS, zPolarACS] = pol2cart(TH_acs,R_acs);
zPolarACS = zPolarACS + z0Polar;

fontSize = 26;
% roi = X >= x_ACS(3) & X <= x_ACS(end) & Z >= z_ACS(1) & Z <= z_ACS(end);


figure('Units','pixels', 'Position', [100, 100, 1200, 600]);

% [ax1,~] = imOverlayPolar(bmode_sam,a_rfm,[-62.5 0],range_acs,0.7, ...
%     xPolar,zPolar,xPolarACS,zPolarACS);
[ax1,~] = imOverlayPolar(bmode_sam,a_rfm,range_bmode,range_acs,fact_transp, ...
    xPolar,zPolar,xPolarACS,zPolarACS);

title(ax1, sprintf('RFM: %.2f ± %.2f, CV=%.2f%%', m_a, s_a, cv_a));
xlabel(ax1, 'Lateral [cm]');
ylabel(ax1, 'Axial [cm]');
set(ax1, 'FontSize', fontSize);
colorbar off
ylim([0 17])
xlim([-11.3 11.3])
hold on;
% contour(xPolar*1e2, zPolar*1e2, roi, 1, 'w--');
% hb2=colorbar; ylabel(hb2,'dB\cdotcm\cdotMHz^{-1}', 'FontSize', fontSize)

% Add label box "FL"
annotation('textbox', [0.28, 0.795, 0.05, 0.10], 'String', 'HL', ...
    'EdgeColor', 'w', 'BackgroundColor', 'k', 'Color', 'w','FontSize', fontSize+6, ...
    'FontWeight', 'bold', 'HorizontalAlignment', 'center');
hold off;
set(ax1,'fontsize',fontSize)

% exportgraphics(gcf,fullfile(dirOut,samName(1:end-4)+"_pol.png"), ...
%     'Resolution','300')

% Plot TNV-RFM
figure('Units','pixels', 'Position', [100, 100, 1200, 600]);

% [ax2, ~] = imOverlayPolar(bmode_sam, a_rfm2, [-62.5 0], range_acs, 0.7, ...
    % xPolar, zPolar, xPolarACS, zPolarACS);
[ax2, ~] = imOverlayPolar(bmode_sam, a_rfm2, range_bmode, range_acs, fact_transp, ...
    xPolar, zPolar, xPolarACS, zPolarACS);

title(ax2, sprintf('TNV-RFM: %.2f ± %.2f, CV=%.2f%%', m_a2, s_a2, cv_a2));
xlabel(ax2, 'Lateral [cm]');
% ylabel(ax2, 'Axial [cm]');
set(ax2, 'FontSize', fontSize);
ylim([0 17])
xlim([-11.3 11.3])
hold on;
% contour(xPolar*1e2, zPolar*1e2, roi, 1, 'w--');
hb2=colorbar; ylabel(hb2,'dB\cdotcm\cdotMHz^{-1}', 'FontSize', fontSize)
% Add label box "FL"
annotation('textbox', [0.28, 0.795, 0.05, 0.10], 'String', 'HL', ...
    'EdgeColor', 'w', 'BackgroundColor', 'k', 'Color', 'w','FontSize', fontSize+6, ...
    'FontWeight', 'bold', 'HorizontalAlignment', 'center');
hold off;
set(ax2,'fontsize',fontSize)

% exportgraphics(gcf,fullfile(dirOut,samName(1:end-4)+"_polTNV.png"), ...
%     'Resolution','300')
end

function t_delay = getRXDelays(Trans, t, n_elements, n_pulses, sound_speed, wvl)

t_delay = zeros(length(t), n_elements, n_pulses);
% (x, z) [m] Obtain positions of center of every element
element_pos_x = Trans.ElementPos(:, 1)*wvl;
element_pos_z = Trans.ElementPos(:, 3)*wvl;
phi = Trans.ElementPos(:, 4);

for n = 1:n_pulses
    for e = 1:n_elements
        focus = sound_speed*t(:)/2;
        xfocus = element_pos_x(n) + sin(phi(n)) * focus;
        zfocus = element_pos_z(n) + cos(phi(n)) * focus;
        t_delay(:,e,n) = (focus + sqrt((zfocus- element_pos_z(e)).^2 + ...
            (xfocus - element_pos_x(e)).^2))/sound_speed;
    end
end

end