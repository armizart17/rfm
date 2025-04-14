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
pars.blocksize   = 12; % wavelengths
pars.blocklines  = 8;
pars.saran_layer = false;
pars.ratio_zx    = 1;
pars.window_type = 3; %  (1) Hanning, (2) Tuckey, (3) Hamming, (4) Tchebychev

blocksize_wv_r = 18;

%% ROI SELECTION (I DID 14 cases)

fprintf('ID N° %d : %s\n', iAcq, samName(1:end-4));

if ( strcmp(samName(1:end-4), '42460445_IHR_F') )
% pars.z_roi = [0.0534482412002264 0.0965834129422875];
% pars.x_roi = [-23.5124207044383 -2.87523485999935];

% pars.z_roi = [0.0464482412002264 0.0965834129422875];
% pars.x_roi = [-24.0124207044383 -2.87523485999935];

pars.z_roi = [0.0564482412002264 0.0965834129422875];
pars.x_roi = [-24.5124207044383 -2.57523485999935];

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

    % figure('Units','centimeters', 'Position',[5 5 15 15]),
    % imagesc(SAM.x, SAM.z*1E3,bmode_sam,range_bmode);
    % colormap gray; clim(range_bmode);
    % hb2=colorbar; ylabel(hb2,'dB')
    % xlabel('Lateral [mm]'), ylabel('Depth [mm]'); 
    % title(caption)
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
%% CHECK ROI

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

%% (Remove later is to check roi)

% fprintf('==========(ROI)===========\n')
fprintf('Axial [mm] = %s;\n', mat2str(round(1e3*pars.z_roi,2)));
fprintf('Lateral [°] = %s;\n', mat2str(round(pars.x_roi,2)));

%%
%% RFM V1 
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
mu          = 1.25;
tau         = 0.0100;
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

        % a_rfm(ii, jj) = mean(slope_full);
        % a_rfm2(ii, jj) = mean(slope2_full);

        a_rfm(ii, jj) = median(slope_full);
        a_rfm2(ii, jj) = median(slope2_full);

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

range_acs = [0 1.5];
fact_transp = 0.7;
fontSize = 15;

figure, 
imagesc(a_rfm, range_acs); % Convert to mm
% % axis("image")
colorbar; colormap("turbo")

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



xFull = SAM.x;
zFull = SAM.z;
[X,Z] = meshgrid(xFull, zFull);

roi = X >= x_ACS(1) & X <= x_ACS(end) & Z >= z_ACS(1) & Z <= z_ACS(end);

figure('Units', 'pixels', 'Position', [100, 100, 1200, 600]); % [x, y, width, height]

tiledlayout(1, 3, 'TileSpacing','compact', 'Padding','compact')

t1 = nexttile();
imagesc(xFull,zFull*1e3, bmode_sam, range_bmode); % axis image;
title('B-mode')
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
%% Plot in cartesian cords
xPolar  = SAM.xp;
zPolar  = SAM.zp;
z0Polar = SAM.z0p; 
r0      = SAM.r0;
r0 = 0.0401;

fontSize = 24;
roi = X >= x_ACS(3) & X <= x_ACS(end) & Z >= z_ACS(1) & Z <= z_ACS(end);

[TH_acs,R_acs] = meshgrid(-x_ACS*pi/180 + pi/2,z_ACS + r0);
[xPolarACS, zPolarACS] = pol2cart(TH_acs,R_acs);
zPolarACS = zPolarACS + z0Polar;

figure('Units','pixels', 'Position', [100, 100, 1200, 600]);


% tiledlayout(1, 2, 'TileSpacing','compact', 'Padding','compact')

% nexttile,
[ax1,~] = imOverlayPolar(bmode_sam,a_rfm,range_bmode,range_acs,0.7, ...
    xPolar,zPolar,xPolarACS,zPolarACS);
% yticks(ax1,[4 8 12 16])
% title(ax1,sprintf('RFM %.3f$\\pm$%.3f, \CV=%.2f%%', m_a, s_a, cv_a), ...
%       'Interpreter', 'latex');
title(ax1,sprintf('RFM %.3f ± %.3f, CV=%.2f%%', m_a, s_a, cv_a));
xlabel(ax1,'Lateral [cm]'), ylabel(ax1,'Axial [cm]')

% xlim([-9 9]),
% ylim([0 15])
hold on
contour(xPolar*1e2, (zPolar)*1e2, roi,1,'w--')
% hb2=colorbar; ylabel(hb2,'dB\cdotcm\cdotMHz^{-1}', 'FontSize', fontSize)
hold off
set(ax1,'fontsize',fontSize)
% exportgraphics(gcf,fullfile(figsDir,samName(1:end-4)+"rfm_polar.png"), ...
%     'Resolution','300')


figure('Units','pixels', 'Position', [100, 100, 1200, 600]);
% nexttile,
[ax1,~] = imOverlayPolar(bmode_sam,a_rfm2,range_bmode,range_acs,0.7, ...
    xPolar,zPolar,xPolarACS,zPolarACS);
colorbar
% yticks(ax1,[4 8 12 16])
% title(ax1,sprintf('TNV-RFM %.3f$\\pm$%.3f, \\%%CV=%.2f', m_a2, s_a2, cv_a2), ...
%       'Interpreter', 'latex');
title(ax1,sprintf('TNV-RFM %.3f ± %.3f, CV=%.2f%%', m_a2, s_a2, cv_a2));
xlabel(ax1,'Lateral [cm]'),
% ylabel(ax1,'Axial [cm]')
set(gca,'fontsize',fontSize)
% xlim([-9 9])
% ylim([0 12])
hold on
contour(xPolar*1e2, (zPolar)*1e2, roi,1,'w--')
hb2=colorbar; ylabel(hb2,'dB\cdotcm\cdotMHz^{-1}', 'FontSize', fontSize)
hold off
set(ax1,'fontsize',fontSize)

% SAVE FIG (RECT)
% exportgraphics(gcf,fullfile(figsDir,samName(1:end-4)+"tnv_polar.png"), ...
%     'Resolution','300')

% close all,
% pause(0.05);
end


%%
% keyboard



%% NOW ESTIMATION TNV

[XX_ACS,ZZ_ACS] = meshgrid(x_ACS, z_ACS);
[Xq,Zq] = meshgrid(x,z);

% Preallocate local attenuation matrix (depth x lateral)
a_rfm1 = zeros(m, n); 
a_rfm2 = zeros(m, n);

% REG
% mu_range = 10.^(0.1:0.05:1);
mu_range = 10.^linspace(0.1, 1.2, 15);

% TNV v1
maxIter1 = 20;
tol1 = 2e-3;

% TNV2 
weights2 = ones(length(band_ufr), 1);
tau2 = 0.01;
maxIter2 = 300;
tol2 = 0.5e-3;
stableIter2 = 10;

% PRELOCATE A MAPS FOR MU_REF
a_rfm_tnv1 = zeros(m,n, length(mu_range));
a_rfm_tnv2 = zeros(m,n, length(mu_range));


tic;
for uu = 1:length(mu_range)
lambda = mu_range(uu);


for jj = 1:n  % Loop over lateral positions (x_j)
    for ii = 1:m  % Loop over depth positions (z_k)

        x_temp = x_temp_all{ii, jj};
        y_temp = y_temp_all{ii, jj};

        X_vec = x_temp(:);

        %%%%%%%%%%%%%%%%%%%%% TNV DENOISING JOINT 1 %%%%%%%%%%%%%%%%%%%%%

        [y_denoised1, cost1, err1, fid1, reg1, iter1] = TNV_regularization(y_temp, ...
            lambda, maxIter1, tol1);
        y_temp1 = y_denoised1;
        %%%%%%%%%%%%%%%%%%%%% TNV DENOISING JOINT 1 %%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%% TNV DENOISING JOINT 2 %%%%%%%%%%%%%%%%%%%%%

        [y_denoised2, cost2, err2, fid2, reg2] = pdo_den_wtnv_1d(y_temp, ...
            lambda, tau2, maxIter2, tol2, stableIter2, weights2);

        y_temp2 = y_denoised2;
        %%%%%%%%%%%%%%%%%%%%% TNV DENOISING JOINT 2 %%%%%%%%%%%%%%%%%%%%%

        y_vec = y_temp1(:);
        a_rfm1(ii, jj) = -(X_vec' * X_vec) \ (X_vec' * y_vec);
        a_rfm1(ii, jj) = abs((X_vec' * X_vec) \ (X_vec' * y_vec));

        y_vec = y_temp2(:);
        a_rfm2(ii, jj) = -(X_vec' * X_vec) \ (X_vec' * y_vec);
        a_rfm2(ii, jj) = abs((X_vec' * X_vec) \ (X_vec' * y_vec));
        % a_rfm(ii, jj) = - cgs(X_vec' * X_vec, X_vec' * y_vec);

    end
end



% SAVE MAP
a_rfm_tnv1(:,:,uu) = a_rfm1;
a_rfm_tnv2(:,:,uu) = a_rfm2;

% BIGGER
AttInterp1 = interp2(XX_ACS,ZZ_ACS, a_rfm1, Xq,Zq);
% AttInterp1 = bigImg(a_rfm1, rfdata_sam_roi); %less std and rmse

AttInterp2 = interp2(XX_ACS,ZZ_ACS, a_rfm2, Xq,Zq);
% AttInterp2 = bigImg(a_rfm2, rfdata_sam_roi); %less std and rmse

% METRICS
r1  = get_metrics_homo_gt(AttInterp1, true(size(AttInterp1)), alpha_sam, 'RFM-TNV1');
r2  = get_metrics_homo_gt(AttInterp2, true(size(AttInterp2)), alpha_sam, 'RFM-TNV2');
r1.mu = lambda;
r2.mu = lambda;

MetricsParam(uu*2-1)   = r1; 
MetricsParam(uu*2)     = r2; 

end

t = toc;
% fprintf('Loop way for mu = %.2f, Elapsed time %.2f \n', lambda, t);
fprintf('Elapsed time %.2f \n', t);

% SAVE DATA
outDir0 = 'D:\emirandaz\qus\rfm\tnv';

outDir = fullfile(outDir0,"a"+alpha_sam);

if ~exist(outDir) mkdir (outDir); end

fileNameOut = "tnv_a"+alpha_sam +".mat";
save(fullfile(outDir, fileNameOut),"a_rfm_tnv1", "a_rfm_tnv2", ...
    "MetricsParam", "mu_range");





%%

acs_range = [0 1.1];
figure, 
sgtitle('RFM TNV 1')
for uu = 1:length(mu_range)
subplot(4, 5, uu)
imagesc(x_ACS, z_ACS, a_rfm_tnv1(:,:,uu), acs_range)
title(['\mu=', num2str(mu_range(uu))]);
colormap("jet")
end

acs_range = [0 1.1];
figure, 
sgtitle('RFM TNV 2')
for uu = 1:length(mu_range)
subplot(4, 5, uu)
imagesc(x_ACS, z_ACS, a_rfm_tnv2(:,:,uu), acs_range)
title(['\mu=', num2str(mu_range(uu))]);
colormap("jet")
end

%%
Tacs        = struct2table(MetricsParam);
Tacs.method = categorical(Tacs.method);

muVec = Tacs(Tacs.method=='RFM-TNV1',:).mu;
tabTNV1 =Tacs(Tacs.method=='RFM-TNV1',:);
tabTNV2 =Tacs(Tacs.method=='RFM-TNV2',:);

colors = lines(2);
lw = 1.5;

figure,
hold on
semilogx((muVec),100*tabTNV1.rmse_homo, '*-', 'LineWidth',lw)
semilogx((muVec),100*tabTNV2.rmse_homo, '*-', 'LineWidth',lw)
xlabel('\mu')
% plot(log10(muVec),100*tabTNV1.rmse_homo, '*-', 'LineWidth',lw)
% plot(log10(muVec),100*tabTNV2.rmse_homo, '*-', 'LineWidth',lw)
% xlabel('log_{10}\mu')
hold off
ylabel('RMSE [%]')
grid on
legend('TNV1','TNV2')
title('RMSE')
% ylim([0 20])
xlim([0 10])


figure,
hold on
errorbar((muVec),tabTNV1.mean_homo,tabTNV1.std_homo, 'LineWidth',lw)
errorbar((muVec),tabTNV2.mean_homo,tabTNV2.std_homo, 'LineWidth',lw)
yline(alpha_sam, 'k--')
% yline(groundTruthTargets(end), '--', 'Color',colors(2,:))
hold off
xlabel('\mu')
ylabel('ACS [dB/cm/MHz]')
grid on
legend('TNV1','TNV2')
title('TNV-RFM')
ylim([0 1.1])




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

SLogRatio_vec = reshape(SLogRatio, [], size(SLogRatio, 3));   
mux_SLogRatio = 1./(abs(mean(SLogRatio_vec, 1, 'omitnan')) ./ std(SLogRatio_vec, 0, 1, 'omitnan') + 1E-5);
weightEstimators = rescale(mux_SLogRatio, 1, max(mux_SLogRatio));

[m, n, ~] = size(SLogRatio);
A1 = kron( 4*zd_zp*1E2*band , speye(m*n) );
A2 = kron( ones(size(band)) , speye(m*n) );
[ACS_SLD, ~] = cgs_ACS([A1 A2], SLogRatio);

%% FIGURES COMPARISON RFM VS SLD

range_acs = [0 1.1];
fontSize = 14;

% ACS_SLD = ACS_RSLD;
figure, 
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1000, 600]); % [x, y, width, height]

% SLD
[m_a, s_a, cv_a] = deal(calc2dStats{1}(ACS_SLD), calc2dStats{2}(ACS_SLD), calc2dStats{3}(ACS_SLD));

subplot(121)
imagesc(x_ACS * 1e3, z_ACS * 1e3, ACS_SLD, range_acs); % Convert to mm
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
imagesc(x_ACS * 1e3, z_ACS * 1e3, a_rfm, range_acs); % Convert to mm
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
ATOT_RFM     = zeros(p, n); 
for ff = 1:p 
    for jj = 1:n
        % Compute FPR
        % FPR = log(squeeze(squeeze(RSp(:,jj,ff)))) / (4*df_MHz)*Np2dB;
        FPR = squeeze(RSp_k_ufr(:,jj,ff)) / (4*df_MHz)*Np2dB;
        
        % Perform linear fit
        [lin_slope , lin_intercept, lin_y, R_fit] = fit_linear(z_ACS_cm, FPR, 2); 
        
        % Store FPR, linear fit, and slope
        FPR_all(:, jj, ff) = FPR;
        FPR_fit_all(:, jj, ff) = lin_y;
        slopes_all(jj, ff) = lin_slope;
        
        % Store ACS value
        % ATOT_RFM(ff, jj) = -lin_slope; % -8.68*m/(4*delta_f) [dB/MHz/cm]
        ATOT_RFM(ff, jj) = -lin_slope; % -8.68*m/(4*delta_f) [dB/MHz/cm]
    end
end

%%
% RESULTS IMAGESC
acs_range = [-1.2 1.2];
figure, 
imagesc(x_ACS*1E2, band, ATOT_RFM, acs_range), colorbar
% hb2=colorbar; ylabel(hb2,'dB/cm/MHz')
hb2=colorbar; ylabel(hb2,'dB\cdotcm^{-1}\cdotMHz^{-1}')
xlabel('Lateral [cm]')
ylabel('Freq [MHz]')

acs_range = [0 1.2];
figure, 
imagesc(x_ACS*1E2, band, abs(ATOT_RFM), acs_range), colorbar
colormap("turbo")
% hb2=colorbar; ylabel(hb2,'dB/cm/MHz')
hb2=colorbar; ylabel(hb2,'dB\cdotcm^{-1}\cdotMHz^{-1}')
xlabel('Lateral [cm]')
ylabel('Freq [MHz]')

freq1 = 1.35; freq2 = 2.55; freq3 = 3.75;
idx_f1 = find(band >= freq1, 1, 'first');
idx_f2 = find(band >= freq2, 1, 'first');
idx_f3 = find(band >= freq3, 1, 'first');

figure, 
plot(x_ACS', ATOT_RFM(idx_f1,:), 'r'), hold on
plot(x_ACS', ATOT_RFM(idx_f2,:), 'g'), hold on
plot(x_ACS', ATOT_RFM(idx_f3,:), 'b'), hold on
% yline(SAM.alpha_value, 'k--')

f_vs_acs= mean((ATOT_RFM), 2);
figure, 
plot(band, f_vs_acs)
yline(mean(f_vs_acs), 'b-')
% yline(SAM.alpha_value, 'k--')
xlabel('Freq')
ylabel('ACS')


%% WITH SLOPE PRINT

% Find freqs MHz
freq1 = 1.35; freq2 = 2.55; freq3 = 3.75;
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
% histogram(ATOT_RFM, 20, 'Normalization', 'probability'); % 20 bins
% xlabel('Value');
% ylabel('Frequency');
% title('Histogram of ACS Elements');
% grid on;

%% BOX PLOT RESULTS RFM

% Compute mean and standard deviation
mu      = mean(ATOT_RFM(:));
sigma   = std(ATOT_RFM(:));
cv      = sigma/mu;

% Generate the box plot
figure;
boxplot(ATOT_RFM(:)); 
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
maxLength  = max([length(ATOT_RFM(:)), length(ACS_SLD(:)), length(ACS_RSLD(:))]);

% Pad each array with NaNs to match the maximum number of rows
ACS_RFM_padded  = [ATT_RFM(:); NaN(maxLength - length(ATT_RFM(:)), 1)];
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
ATT_RFM = zeros(p, n);

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
        ATT_RFM(ff, jj) = -lin_slope; % -8.68*m/(4*delta_f) [dB/MHz/cm]

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
mu      = mean(ATT_RFM(:));
sigma   = std(ATT_RFM(:));
cv      = sigma/mu;

% Generate the box plot
figure;
boxplot(ATT_RFM(:)); 
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

range_acs = [0 1.1];
fontSize = 14;

figure;
subplot(121)
imagesc(x_ACS * 1e3, z_ACS * 1e3, a_local, range_acs); % Convert to mm
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