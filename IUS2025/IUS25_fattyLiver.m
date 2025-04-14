% Reference Frequency Method vApril
% Here I try a simulation with different mu_values for TNV method
% Method by gpt and EMZ 
clear all;
% clear all, clc, close all;
alpha_gt = NaN;
manualroi = false;
denTNVRFM = false;

Np2dB       = 20*log10(exp(1));
dB2Np       = 1/Np2dB;
range_bmode = [-60 0];
range_acs   = [0 1.2];
fact_transp = 0.7;
fontSize    = 15;

mean2d = @(x) mean(x(:));
std2d = @(x) std(x(:));
cv2d = @(x) 100*std(x(:))/mean(x(:));

calc2dStats = {@(x) mean(x(:)), @(x) std(x(:)), @(x) 100 * std(x(:)) / mean(x(:))};

%% LOAD SAM

% DATA NEW AMZ
% 

% DATA LIM PC
pathData = 'D:\emirandaz\qus\data\liver\bf_M04_D02';
figsDir = 'D:\emirandaz\qus\data\liver\rfm_v2';

if ~exist(figsDir) mkdir(figsDir); end

acqDir = dir(fullfile(pathData,'*.mat'));

list_acqDir = [1 2 5 6 7 8 9 11 12 15 16 17 18 20];
list_acqDir = [ 1 ]; %1 6 11 15
for iList = 1:length(list_acqDir)

% Each 5
% iAcq = 1;
iAcq = list_acqDir (iList);
samName = acqDir(iAcq).name;
fileName = fullfile(pathData, samName);

%%%%%%%%%%%%%%% NEW MARCH %%%%%%%%%%%%%%%
SAM = load(fileName);
%%%%%%%%%%%%%%% NEW MARCH %%%%%%%%%%%%%%%


% B-MODE CHECK (PROBABLY NEED PADDING)
% bmode_sam = db(hilbert(SAM.rf));
% bmode_sam = bmode_sam - max(bmode_sam(:));

bmode_sam   = SAM.bMode;
SAM.x       = SAM.xr;
SAM.z       = SAM.zr;

caption = strrep(samName(1:end-4), '_', ' ');

% figure,
% imagesc(SAM.x, SAM.z*1E3, bmode_sam, range_bmode), 
% xlabel('Degree [º]'), ylabel('Depth [mm]');
% cb = colorbar;
% cb.Label.String = 'dB'; % Add the label "dB"
% title(caption)
% % axis("image");
% colormap('gray')

%% SPECTRAL PARAMETERS
pars.bw          = [1.2 3.9]; % [MHz]
pars.overlap     = 0.8;
pars.blocksize   = 11; % wavelengths
pars.blocklines  = 6;
pars.saran_layer = false;
pars.ratio_zx    = 1;
pars.window_type = 3; %  (1) Hanning, (2) Tuckey, (3) Hamming, (4) Tchebychev

blocksize_wv_r = 24;

%% ROI SELECTION (I DID 14 cases)

fprintf('ID N° %d : %s\n', iAcq, samName(1:end-4));

if ( strcmp(samName(1:end-4), '42460445_IHR_F') )
% pars.z_roi = [0.0534482412002264 0.0965834129422875];
% pars.x_roi = [-23.5124207044383 -2.87523485999935];

% pars.z_roi = [0.0464482412002264 0.0965834129422875];
% pars.x_roi = [-24.0124207044383 -2.87523485999935];

% ORIGINAL
pars.z_roi = [0.0564482412002264 0.0965834129422875]; % 0.0401
pars.x_roi = [-24.4124207044383 -2.67523485999935]; % 21.73

% pars.z_roi = [0.0534 0.0965];
% pars.x_roi = [-24.412 -2.375];
% 
% pars.z_roi = [55.05 96.5]*1e-3;
% pars.x_roi = [-24.412 -2.375];
% BIG IUS SIMILAR SIZE HEALTHY (DEPRECATED)
% Diff z_roi = 0.0462
% Diff x_roi = 29.75
% pars.x_roi = [-27.41   -2.2198];
% pars.z_roi = [0.050    0.0985];

elseif ( strcmp(samName(1:end-4), '42460445_IOLAI_F') )
pars.z_roi = [0.06265 0.13];
pars.x_roi = [-27 -5.8];

elseif ( strcmp(samName(1:end-4), '42460445_PL_F') )
pars.z_roi = [0.04765518886002839 0.0887218844567291];
pars.x_roi = [-20.3157416733283 1.23992789555525];

elseif ( strcmp(samName(1:end-4), '44730576_IHR_F') )
pars.z_roi = [0.0432077721719176 0.0787915326854975];
pars.x_roi = [-25.0846951177715 3.52902554888802];

% pars.z_roi = [0.023802428576483 0.0705162395428045];
% pars.x_roi = [-22.2233230511056 10.587076646664];

elseif ( strcmp(samName(1:end-4), '44730576_IOLAI_F') )
pars.z_roi = [0.0361737730006286 0.0638960050286501];
pars.x_roi = [-8.10722085555353 18.9804347088842];

elseif ( strcmp(samName(1:end-4), '44730576_LCM_F') )
pars.z_roi = [0.0456903601147255 0.074653886114151];
pars.x_roi = [-15.3560300911073 9.0610115444422];

% pars.z_roi = [0.0448628308004562 0.0684474162571313];
% pars.x_roi = [-13.44844871333 10.2055603711086];

elseif ( strcmp(samName(1:end-4), '44730576_LHI_F') )
pars.z_roi = [0.0485867127146681 0.0763089447426896];
pars.x_roi = [-21.4602904999946 4.6735743755544];

elseif ( strcmp(samName(1:end-4), '70141854_IHR_F') )
pars.z_roi = [0.0319664782005137 0.0771364740569589];
pars.x_roi = [-17.4543696066623 12.8761742999968];

% pars.z_roi = [0.0382425962863018 0.0643097696857848];
% pars.x_roi = [-13.2576905755522 18.0266440199955];

elseif ( strcmp(samName(1:end-4), '70141854_IOLAI_F') )
pars.z_roi = [0.0419664782005137 0.0920320017138063];
pars.x_roi = [-13.6392068511077 11.350109197775];

elseif ( strcmp(samName(1:end-4), '70141854_PL_F') ) % (BESSELS)
pars.z_roi = [0.0440353014861869 0.091204472399537];
pars.x_roi = [-8.10722085555353 22.7955974644388];

elseif ( strcmp(samName(1:end-4), '70897309_IHR_F') ) 
pars.z_roi = [0.0382425962863018 0.0725850628284778];
pars.x_roi = [-24.7031788422161 0.476895344444337];

elseif ( strcmp(samName(1:end-4), '70897309_IOLAI_F') ) % LIKE IT
pars.z_roi = [0.042794007514783 0.0870668258281905];
pars.x_roi = [-18.5989184333287 11.5408673355527];

elseif ( strcmp(samName(1:end-4), '70897309_LCM_F') ) 
pars.z_roi = [0.0560344765430918 0.120581763056097];
pars.x_roi = [-25.2754532555492 -1.81220230888842];

elseif ( strcmp(samName(1:end-4), '70897309_PL_F') ) 
pars.z_roi = [0.0432077721719176 0.112720234570539];
pars.x_roi = [-28.7090997355484 -5.62736506444304];

else 

    figure('Units','centimeters', 'Position',[5 5 15 15]),
    imagesc(SAM.x, SAM.z*1E3,bmode_sam,range_bmode);
    colormap gray; clim(range_bmode);
    hb2=colorbar; ylabel(hb2,'dB')
    xlabel('Lateral [mm]'), ylabel('Depth [mm]'); 
    title(caption)

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

    pars.x_roi     = [rect(1), rect(1)+rect(3)];      % [º]
    pars.z_roi     = [rect(2), rect(2)+rect(4)]*1E-3; % [m]
    continue;

end

%%
if manualroi 
figure('Units','centimeters', 'Position',[5 5 15 15]),
imagesc(SAM.x, SAM.z*1E3,bmode_sam,range_bmode);
colormap gray; clim(range_bmode);
hb2=colorbar; ylabel(hb2,'dB')
xlabel('Lateral [mm]'), ylabel('Depth [mm]'); 
title(caption)

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

pars.x_roi     = [rect(1), rect(1)+rect(3)];      % [º]
pars.z_roi     = [rect(2), rect(2)+rect(4)]*1E-3; % [m]
end
% CHECK ROI

% figure,
% imagesc(SAM.x, SAM.z*1E3, bmode_sam,range_bmode), 
% hb2=colorbar; ylabel(hb2,'dB')
% rectangle('Position', [pars.x_roi(1) 1E3*pars.z_roi(1) pars.x_roi(2)-pars.x_roi(1) 1E3*(pars.z_roi(2)-pars.z_roi(1))], ...
%         'EdgeColor','r', 'LineWidth', 2, 'LineStyle','--'), hold off;
% xlabel('Lateral [º]'), 
% ylabel('Depth [mm]');
% title(caption)
% % axis("image");
% colormap('gray')

% (Remove later is to check roi)

% % fprintf('==========(ROI)===========\n')
% fprintf('Axial [mm] = %s;\n', mat2str(round(1e3*pars.z_roi,2)));
% fprintf('Lateral [°] = %s;\n', mat2str(round(pars.x_roi,2)));

%
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
mu          = 1.15;
tau         = 0.01000;
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
% vec_gamma = 1 + (1.1- 1) * (1 - linspace(0, 1, p_ufr).^gamma);
% vec_log = logspace(log10(10), log10(1), p_ufr);
% figure, plot(vec_gamma)
% 
% weigthChannels = vec_gamma;

[RSp_k_ufr_opt, cost, error, fid, reg] = pdo_den_wtnv(RSp_k_ufr, mu, tau, maxIter, tol, stableIter, weigthChannels);

[RSp_r_ufr_opt, cost, error, fid, reg] = pdo_den_wtnv(RSp_r_ufr, mu, tau, maxIter, tol, stableIter, weigthChannels);

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
        a_rfm(ii, jj) = -(X_vec' * X_vec) \ (X_vec' * y_vec);
        % a_rfm(ii, jj) = abs((X_vec' * X_vec) \ (X_vec' * y_vec));
        % 
        a_rfm2(ii, jj) = -(X_vec' * X_vec) \ (X_vec' * y_vec2);
        % a_rfm2(ii, jj) = abs((X_vec' * X_vec) \ (X_vec' * y_vec2));


        % % Option 2
        % A2_full    = kron(speye(p_ufr), x_temp(:, 1));
        % A2t_full   = A2_full';
        % 
        % slope_full = -(cgs(A2t_full*A2_full, A2t_full*y_vec));
        % slope2_full = -(cgs(A2t_full*A2_full, A2t_full*y_vec2));
        % 
        % a_rfm(ii, jj) = mean(slope_full);
        % a_rfm2(ii, jj) = mean(slope2_full);

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
% RFM
[m_a, s_a, cv_a]    = deal(calc2dStats{1}(a_rfm), calc2dStats{2}(a_rfm), calc2dStats{3}(a_rfm));
[m_a2, s_a2, cv_a2] = deal(calc2dStats{1}(a_rfm2), calc2dStats{2}(a_rfm2), calc2dStats{3}(a_rfm2));


% VISUALIZATION RECT
fontSize = 14;
xFull = SAM.x;
zFull = SAM.z;
[X,Z] = meshgrid(xFull, zFull);

roi = X >= x_ACS(1) & X <= x_ACS(end) & Z >= z_ACS(1) & Z <= z_ACS(end);

figure('Units', 'pixels', 'Position', [100, 100, 1200, 600]); % [x, y, width, height]

tiledlayout(1, 3, 'TileSpacing','compact', 'Padding','compact')

t1 = nexttile();
imagesc(xFull,zFull*1e3, bmode_sam, range_bmode); % axis image;
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

keyboard

%% MULTI PREV TEST
dirMulti = ['D:\emirandaz\qus\rfm\IUS2025\FL\multil_',num2str(pars.blocklines)];
if ~exist(dirMulti) mkdir(dirMulti); end

% ['Bmode k', mat2str(round(pars.blocksize,2)), ...
%     'r', mat2str(round(blocksize_wv_r,2)), ...
%     'L', mat2str(round(pars.blocklines,2))]
save_all_figures_to_directory(dirMulti,char("M"+pars.blocksize+"fig"))
close all

%%
dirOut = 'D:\emirandaz\qus\rfm\IUS2025\abstract\fl';
if ~exist(dirOut) mkdir(dirOut); end


%% Plot in cartesian cords
% Plot RFM
xPolar  = SAM.xp;
zPolar  = SAM.zp;
z0Polar = SAM.z0p; 
r0      = SAM.r0;
r0 = 0.0401;
fact_transp = 0.45;

fontSize = 26;


roi = X >= x_ACS(1) & X <= x_ACS(end) & Z >= z_ACS(1) & Z <= z_ACS(end);

% x_ACS2 = x_ACS;
% x_ACS2(1) = x_ACS(1)-0.1;

% z_ACS2 = z_ACS;
% z_ACS2(1) = z_ACS(1)-0.1;

x_ACS2 = linspace(x_ACS(1)-0.1, x_ACS(end)+0.1, length(x_ACS))';
z_ACS2 = linspace(z_ACS(1)-0.0007, z_ACS(end)+0.0007, length(z_ACS))';

% [TH_acs,R_acs] = meshgrid(-x_ACS*pi/180 + pi/2,z_ACS + r0);
[TH_acs,R_acs] = meshgrid(-x_ACS2*pi/180 + pi/2,z_ACS2 + r0);

[xPolarACS, zPolarACS] = pol2cart(TH_acs,R_acs);
zPolarACS = zPolarACS + z0Polar;


figure('Units','pixels', 'Position', [100, 100, 1200, 600]);

[ax1,~] = imOverlayPolar(bmode_sam,a_rfm,range_bmode,range_acs,fact_transp, ...
    xPolar,zPolar,xPolarACS,zPolarACS);

title(ax1, sprintf('RFM: %.2f ± %.2f, CV=%.2f%%', m_a, s_a, cv_a));
xlabel(ax1, 'Lateral [cm]');
ylabel(ax1, 'Axial [cm]');
set(ax1, 'FontSize', fontSize);
ylim([-1 17])
xlim([-12 12])
hold on;
% contour(xPolar*1e2, zPolar*1e2, roi, 1, 'w--');
% hb2=colorbar; ylabel(hb2,'dB\cdotcm\cdotMHz^{-1}', 'FontSize', fontSize)

% Add label box "FL"
annotation('textbox', [0.28, 0.785, 0.05, 0.10], 'String', 'FL', ...
    'EdgeColor', 'w', 'BackgroundColor', 'k', 'Color', 'w','FontSize', fontSize+6, ...
    'FontWeight', 'bold', 'HorizontalAlignment', 'center');
hold off;
set(ax1,'fontsize',fontSize)


% exportgraphics(gcf,fullfile(dirOut,samName(1:end-4)+"_pol.png"), ...
%     'Resolution','300')
% Plot TNV-RFM
figure('Units','pixels', 'Position', [100, 100, 1200, 600]);

[ax2, ~] = imOverlayPolar(bmode_sam, a_rfm2, range_bmode, range_acs, fact_transp, ...
    xPolar, zPolar, xPolarACS, zPolarACS);

title(ax2, sprintf('TNV-RFM: %.2f ± %.2f, CV=%.2f%%', m_a2, s_a2, cv_a2));
xlabel(ax2, 'Lateral [cm]');
% ylabel(ax2, 'Axial [cm]');
set(ax2, 'FontSize', fontSize);
ylim([-1 17])
xlim([-12 12])
hold on;
% contour(xPolar*1e2, zPolar*1e2, roi, 1, 'w--');
% hb2=colorbar; ylabel(hb2,'dB\cdotcm\cdotMHz^{-1}', 'FontSize', fontSize)

% Add label box "FL"
annotation('textbox', [0.28, 0.785, 0.05, 0.10], 'String', 'FL', ...
    'EdgeColor', 'w', 'BackgroundColor', 'k', 'Color', 'w','FontSize', fontSize+6, ...
    'FontWeight', 'bold', 'HorizontalAlignment', 'center');
hold off;
set(ax2,'fontsize',fontSize)

% exportgraphics(gcf,fullfile(dirOut,samName(1:end-4)+"_polTNV.png"), ...
%     'Resolution','300')
end


%%
% keyboard



