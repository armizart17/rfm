% Reference Frequency Method v2 March
% AMZ 

clear all, 
close all;

Np2dB       = 20*log10(exp(1));
dB2Np       = 1/Np2dB;
range_bmode = [-60 0];
range_acs   = [0 1.2];
manualroi   = false;
fontSize = 26;

dataVerasonics = true;
dataSonix      = false;

mean2d = @(x) mean(x(:));
std2d = @(x) std(x(:));
cv2d = @(x) 100*std(x(:))/mean(x(:));
calc2dStats = {@(x) mean(x(:)), @(x) std(x(:)), @(x) 100 * std(x(:)) / mean(x(:))};

list_Phantom = ["", "_2", "_3"];
iPhantom = 3;

[ret, pcname] = system('hostname');

if strcmp(pcname(1:end-1), 'C084285') % PC LIM
    baseDir = 'D:\emirandaz\qus\data';
else % PC EMZ
    baseDir = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\';
end
%% DATA VERASONICS
if dataVerasonics
pathData = fullfile(baseDir, 'SavedDataQUSPhantom\bf'); 

folderDataSam = '544'; numPhantomSam = '544'; alpha_sam = 0.53; sos_sam = 1539;
% folderDataSam = '261'; numPhantomSam = '261'; alpha_sam = 0.54; sos_sam = 1509;
samName = numPhantomSam + "_F" + list_Phantom(iPhantom);

folderDataRef = '261'; numPhantomRef = '261'; alpha_ref = 0.48; sos_ref = 1509;
% folderDataRef = '544'; numPhantomRef = '544'; alpha_ref = 0.53; sos_ref = 1539;

% refName = numPhantomRef + "_F";

filesRef = dir(fullfile(pathData, folderDataRef,'*.mat'));
end
%% DATA SONIX
if dataSonix
% folderDataSam = 'ID261V2'; numPhantomSam = 261; 
folderDataSam = 'ID544V2'; numPhantomSam = 544; alpha_sam = 0.53; sos_sam = 1539;

filesSam = dir(fullfile(pathData, folderDataSam,'*.mat'));
% samList = ["16-19-52","16-20-50","16-21-22", "16-21-22"];
samName = filesSam(1).name;

folderDataRef = 'ID544V2'; numPhantomSam = 544; alpha_ref = 0.53; sos_ref = 1539;
folderDataRef = 'ID261V2'; numPhantomSam = 261; alpha_ref = 0.48; sos_ref = 1509;

filesRef = dir(fullfile(pathData, folderDataRef,'*.mat'));
end
%% LOAD DDATA

SAM = load (fullfile(pathData, folderDataSam, samName));
SAM.rf = SAM.rf(:,:,1);
SAM.c0 = sos_sam;

gt_acs  = alpha_sam;

% B-MODE CHECK
bmode_sam = db(hilbert(SAM.rf));
bmode_sam = bmode_sam - max(bmode_sam(:));

numRefs  = length(filesRef); 
REF      = load( fullfile(pathData, folderDataRef, filesRef(1).name ) );
newrf  = nan([size(REF.rf), numRefs], 'like', REF.rf); % Use 'like' for type consistency
for i = 1:numRefs
    newrf(:,:,i) = load(fullfile(pathData,folderDataRef,filesRef(i).name ), 'rf').rf(:,:,1); % Directly extract rf, avoiding redundant variables
end

REF.rf = newrf;
REF.acs = alpha_ref;
REF.c0  = sos_ref;

%% SPECTRAL PARAMETERS

pars.bw          = [3 8.4]; % [MHz]
pars.overlap     = 0.8;
pars.blocksize   = 10; % wavelengths
pars.z_roi       = [4 36]*1E-3; % [m] 
pars.x_roi       = [-15 15]*1E-3; % [m] 
pars.saran_layer = false;
pars.ratio_zx    = 1.25;
pars.window_type = 3; %  (1) Hanning, (2) Tuckey, (3) Hamming, (4) Tchebychev

blocksize_wv_r = 20;

%% ROI SELECTION
if ~manualroi

    % pars.z_roi       = [5 35]*1E-3; % [m] 
    pars.z_roi       = [25 50]*1E-3; % [m] 
    pars.x_roi       = [-17 17]*1E-3; % [m] 
        
else 

    figure('Units','centimeters', 'Position',[5 5 15 15]),
    imagesc(SAM.x*1E3, SAM.z*1E3,bmode_sam,range_bmode);
    colormap gray; clim(range_bmode);
    hb2=colorbar; ylabel(hb2,'dB')
    xlabel('Lateral [mm]'), ylabel('Depth [mm]'); 
    title('Bmode')
    
    confirmation = '';
    while ~strcmp(confirmation,'Yes')
        rect = getrect;
        confirmation = questdlg('Sure?');
        if strcmp(confirmation,'Cancel')
            disp(rect)
            break
        end
    end
    close,

    pars.x_roi     = [rect(1), rect(1)+rect(3)]*1E-3; % [m]
    pars.z_roi     = [rect(2), rect(2)+rect(4)]*1E-3; % [m]
end

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
RSp_k(:,:, 2:end) = Sp_k(:,:, 2:end) ./ Sp_k(:,:, 1:end-1);
% For the first slice, keep the ratio the same as the first slice
RSp_k(:,:, 1) = RSp_k(:,:, 2); % Assuming the first slice ratios are 1 as there is no "i-1"
RSp_k = log(RSp_k); % @@
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
RSp_r(:,:, 2:end) = Sp_r(:,:, 2:end) ./ Sp_r(:,:, 1:end-1);
% For the first slice, keep the ratio the same as the first slice
RSp_r(:,:, 1) = RSp_r(:,:, 2); % Assuming the first slice ratios are 1 as there is no "i-1"
RSp_r = log(RSp_r); % @@
%% ATTEMPT RFM ACS
% USE THIS FROM APRIL AND NOW ON
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
            % y_col = squeeze( ( log(RSp_k_ufr(ii, jj, :)) - log(RSp_r_ufr(r, jj, :)) ) /(4*df_MHz) *Np2dB ); % p_ufr x 1
            y_col = squeeze( ( (RSp_k_ufr(ii, jj, :)) - (RSp_r_ufr(r, jj, :)) ) /(4*df_MHz) *Np2dB ); %@@ p_ufr x 1

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

% NOW ESTIMATION 

[XX_ACS,ZZ_ACS] = meshgrid(x_ACS, z_ACS);
[Xq,Zq] = meshgrid(x,z);

% PRELOCATE A MAPS FOR MU_REF
a_rfm = zeros(m,n);

tic;
for jj = 1:n  % Loop over lateral positions (x_j)
    for ii = 1:m  % Loop over depth positions (z_k)

        x_temp = x_temp_all{ii, jj};
        y_temp = y_temp_all{ii, jj};

        X_vec = x_temp(:);
        y_vec = y_temp(:);
       
        a_rfm(ii, jj) = -(X_vec' * X_vec) \ (X_vec' * y_vec);

    end
end
t = toc;
fprintf('Loop way for mu = %.2f, Elapsed time %.2f \n', lambda, t);


% IMOVERLAY SIMPLE

% Simple Script example to overlay RESIZED colorImg to bmodeFull
% Function requires previous resizing (i.e. bigImg)

% RFM
[m_a, s_a, cv_a] = deal(calc2dStats{1}(a_rfm), calc2dStats{2}(a_rfm), calc2dStats{3}(a_rfm));


zOffset = 1.85;

units           = 1E2;
bmodeFull       = bmode_sam;
colorImg        = bigImg(a_rfm, rfdata_sam_roi );   
range_bmode     = [-60 0];
range_img       = [0 1.1];
transparency    = 0.65;
x_img           = x_ACS*units;
z_img           = z_ACS*units - zOffset;
xFull           = SAM.x*units;
zFull           = SAM.z*units;

figure, 
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1000, 600]); % [x, y, width, height]
[~,hB,hColor] = imOverlaySimple(bmodeFull, colorImg, range_bmode, ...
                range_img, transparency, x_img, z_img, xFull, zFull);

xlabel('Lateral [cm]'), ylabel('Depth [cm]');
hColor.Label.String = 'dB\cdotcm^{-1}\cdotMHz^{-1}';
title(sprintf('RFM: %.2f ± %.2f, CV = %.2f%%', ...
               m_a, s_a, cv_a));
axis("image")
ylim([0.1 3.5])
yticks([0 1 2 3 4])
set(gca,'fontsize',fontSize)

% METRICS
m_RFM = get_metrics_homo_gt(a_rfm, true(size(a_rfm)), alpha_sam, 'RFM');

%% TNV-RFM

% DENOISING TNV RSP
mu          = 1;
tau         = 0.01;
maxIter     = 1000;
stableIter  = 20;
% tol         = 0.5e-4; % tolerance error
tol         = 1e-3; % tolerance error
RSp_k_ufr_vec = reshape(RSp_k_ufr, [], size(RSp_k_ufr, 3));   
mux_RSp = 1./(abs(mean(RSp_k_ufr_vec, 1, 'omitnan')) ./ std(RSp_k_ufr_vec, 0, 1, 'omitnan') + 1E-5);
% weigthChannels = rescale(mux_RSp, 1, 10);
% weigthChannels = mux_RSp;
weigthChannels = ones(1, p_ufr);

% gamma = 1;
% vec_gamma = 1 + (2 - 1) * (1 - linspace(0, 1, p_ufr).^gamma);
% vec_log = logspace(log10(10), log10(1), p_ufr);
% weigthChannels = vec_gamma;

% figure, plot(weigthChannels)

[RSp_k_ufr_opt, cost, error, fid, reg] = pdo_den_wtnv(RSp_k_ufr, mu, tau, maxIter, tol, stableIter, weigthChannels);

[RSp_r_ufr_opt, cost, error, fid, reg] = pdo_den_wtnv(RSp_r_ufr, mu, tau, maxIter, tol, stableIter, weigthChannels);

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
            % y_col = squeeze( ( log(RSp_k_ufr(ii, jj, :)) - log(RSp_r_ufr(r, jj, :)) ) /(4*df_MHz) *Np2dB ); % p_ufr x 1
            y_col = squeeze( ( (RSp_k_ufr_opt(ii, jj, :)) - (RSp_r_ufr_opt(r, jj, :)) ) /(4*df_MHz) *Np2dB ); %@@ p_ufr x 1

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

% NOW ESTIMATION 

[XX_ACS,ZZ_ACS] = meshgrid(x_ACS, z_ACS);
[Xq,Zq] = meshgrid(x,z);

% PRELOCATE A MAPS FOR MU_REF
a_tnv_rfm = zeros(m,n);

tic;
for jj = 1:n  % Loop over lateral positions (x_j)
    for ii = 1:m  % Loop over depth positions (z_k)

        x_temp = x_temp_all{ii, jj};
        y_temp = y_temp_all{ii, jj};

        X_vec = x_temp(:);
        y_vec = y_temp(:);
       
        a_tnv_rfm(ii, jj) = -(X_vec' * X_vec) \ (X_vec' * y_vec);

    end
end
t = toc;
fprintf('Loop way for mu = %.2f, Elapsed time %.2f \n', lambda, t);


% IMOVERLAY SIMPLE

% Simple Script example to overlay RESIZED colorImg to bmodeFull
% Function requires previous resizing (i.e. bigImg)

% RFM
[m_a, s_a, cv_a] = deal(calc2dStats{1}(a_tnv_rfm), calc2dStats{2}(a_tnv_rfm), calc2dStats{3}(a_tnv_rfm));


zOffset = 1.85;

units           = 1E2;
bmodeFull       = bmode_sam;
colorImg        = bigImg(a_tnv_rfm, rfdata_sam_roi );   
range_bmode     = [-60 0];
range_img       = [0 1.1];
transparency    = 0.65;
x_img           = x_ACS*units;
z_img           = z_ACS*units - zOffset;
xFull           = SAM.x*units;
zFull           = SAM.z*units;

figure, 
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1000, 600]); % [x, y, width, height]
[~,hB,hColor] = imOverlaySimple(bmodeFull, colorImg, range_bmode, ...
                range_img, transparency, x_img, z_img, xFull, zFull);

xlabel('Lateral [cm]'), ylabel('Depth [cm]');
hColor.Label.String = 'dB\cdotcm^{-1}\cdotMHz^{-1}';
title(sprintf('TNV-RFM: %.2f ± %.2f, CV = %.2f%%', ...
               m_a, s_a, cv_a));
axis("image")
ylim([0.1 3.5])
yticks([0 1 2 3 4])
set(gca,'fontsize',fontSize)

% METRICS
m_TNVFRM = get_metrics_homo_gt(a_tnv_rfm, true(size(a_tnv_rfm)), alpha_sam, 'TNV-RFM');

%%

% Example: assuming a_rfm and a_tnv_rfm are vectors or matrices
% Flatten them to column vectors if necessary
a_rfm = a_rfm(:);
a_tnv_rfm = a_tnv_rfm(:);

% Combine data and group labels
all_data = [a_rfm; a_tnv_rfm];
group = [repmat({'RFM'}, length(a_rfm), 1); repmat({'TNV-RFM'}, length(a_tnv_rfm), 1)];

% Create boxplot
figure;
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 800, 600]); % [x, y, width, height]
boxplot(all_data, group);
% axis("equal")
yline(alpha_sam, 'k--')

set(gca, 'XTickLabel', {'RFM', 'TNV-RFM'}, 'FontWeight', 'bold');
ylabel('ACS [dB\cdotcm^{-1}\cdotMHz^{-1}]');
% xlabel('\bfMethods');
title('\bfLiver Phantom')
set(gca,'fontsize',fontSize+2)
grid on;


%% SAVE FIG
dirOut = 'D:\emirandaz\qus\rfm\IUS2025\abstract\phantomLIM';
if ~exist(dirOut) mkdir(dirOut); end

nameFig = ['LIM_P', num2str(iPhantom),'_'];

save_all_figures_to_directory(dirOut, nameFig, 'svg')

%% SAVE TABLE METRICS

Metrics = [m_RFM; m_TNVFRM]; 

T        = struct2table(Metrics);
T.method = categorical(T.method);

% Write to Excel
nameExcel = ['LIM_P', num2str(iPhantom)]; 

excelFile = fullfile(dirOut, nameExcel+".xlsx");

writetable(T, excelFile, 'Sheet', 'Metrics', 'WriteRowNames', true);

