% Note: Remember SLD 
% acs_ref     = 0.4;                % [dB/cm/MHz]
% att_ref     = acs_ref*band/Np2dB; % vector [Np/cm]
% att_ref_map = ones(size(Sp_sam)) .* reshape(att_ref, 1, 1, []); % 3D array as SLogRatio
% sld         = log(Sp_sam./Sd_sam) - ( log(Sp_ref./Sd_ref) - 4*zd_zp*1E2*att_ref_map ); % (rows, cols, freqChannels)

%%
%%
clear all, clc, close all;

Np2dB       = 20*log10(exp(1));
dB2Np       = 1/Np2dB;
%% LOAD SAM
pathData = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\dataACS_kwave';

% SAM = load(fullfile(pathData, 'sam0.mat'));
% REF = load(fullfile(pathData, 'ref0.mat'));

SAM = load(fullfile(pathData, 'ref0.mat'));

% B-MODE CHECK
bmode_sam = db(hilbert(SAM.rf));
bmode_sam = bmode_sam - max(bmode_sam(:));

figure,

% subplot(121), 
imagesc(SAM.x*1E3, SAM.z*1E3, bmode_sam), axis("image");
xlabel('Lateral [mm]'), ylabel('Depth [mm]');
title('SAM')
colormap('gray')


%%
pars.bw          = [3 9]; % [MHz]
pars.overlap     = 0.8;
pars.blocksize   = 20; % wavelengths
pars.z_roi       = [5 45]*1E-3; % [m] 
pars.x_roi       = [-17 17]*1E-3; % [m] 
pars.saran_layer = false;

pars.window_type = 3; %  (1) Hanning, (2) Tuckey, (3) Hamming, (4) Tchebychev

figure,
imagesc(SAM.x*1E3, SAM.z*1E3, bmode_sam), axis("image");
rectangle('Position', 1E3*[pars.x_roi(1) pars.z_roi(1) pars.x_roi(2)-pars.x_roi(1) pars.z_roi(2)-pars.z_roi(1)], ...
        'EdgeColor','r', 'LineWidth', 2, 'LineStyle','--'), hold off;
xlabel('Lateral [mm]'), ylabel('Depth [mm]');
title('SAM')
colormap('gray')

%%
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
if isfield(SAM, 'rf')
    rfdata_sam   = SAM.rf;
end
if isfield(SAM, 'RF')
    rfdata_sam   = SAM.RF;
end

%% Spectral Ratio INIT
z            = SAM.z;
x            = SAM.x;
fs           = SAM.fs;

dx = x(2)-x(1);
dz = z(2)-z(1);
c0 = 1540; % [m/s]
lambda = c0/mean(bw)*1E-6;

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
ratio_zx = 1; % @
wz = round(blocksize_wv*lambda*(1-overlap)/dz * ratio_zx); % Between windows
nz = 2*round(blocksize_wv*lambda/dz /2 * ratio_zx); % Window size

z0p = 1:wz:length(z)-nz;
z0d = z0p + nz/2;
z_ACS = z(z0p+ nz/2);
m  = length(z0p);

 
% Frequency samples
NFFT = 2^(nextpow2(nz/2)+1);
NFFT = 4096;

% axis_f = (0:(NFFT/2-1))'*fs/NFFT;  % [Hz] (so 0-fs/2),  it should be
axis_f = (0:NFFT-1)'/NFFT * fs;   % [Hz] freq axis as default because "spectra" function
freq_L = bw(1)*1E6; % [Hz] 
freq_H = bw(2)*1E6; % [Hz]

ind_f = axis_f >= freq_L & axis_f <= freq_H ;   

band  = axis_f(ind_f)*1E-6; % [MHz]
p = length(band);

zd_zp = (nz/2)*dz;   % zd_zp = 2*\delta_z = nz/2 =  %  [m]

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

windowing = window_choice(nz/2, window_type);
windowing = windowing*ones(1,nx);

nSamples = size(rfdata_sam_roi,3);
Sp_ref = zeros(m,n,p,nSamples);
Sd_ref = zeros(m,n,p,nSamples);

% € 
SNR = zeros(m,n,nSamples);

%%
% loop in case many Acq 
    % (useful if many frames in reference phantom)
    for iRef = 1:nSamples   
        samRef = rfdata_sam_roi(:,:,iRef);

        envelope = abs(hilbert(samRef)); % € 
        for jj=1:n
            for ii=1:m
                xw = x0(jj) ;   % x window
                zp = z0p(ii);
                zd = z0d(ii);
    
                sub_block_p = samRef(zp:zp+nz/2-1,xw:xw+nx-1);
                sub_block_d = samRef(zd:zd+nz/2-1,xw:xw+nx-1);

                [tempSp,~] = spectra(sub_block_p,windowing,0,nz/2,NFFT);
                [tempSd,~] = spectra(sub_block_d,windowing,0,nz/2,NFFT);
    
                Sp_ref(ii,jj,:,iRef) = tempSp(ind_f) ./ t_saran;
                Sd_ref(ii,jj,:,iRef) = tempSd(ind_f) ./ t_saran;

                % % € 
                % sub_env_p  = envelope(zp:zp+nz/2-1,xw:xw+nx-1);
                % sub_env_d  = envelope(zd:zd+nz/2-1,xw:xw+nx-1);
                % temp = [sub_env_p(:); sub_env_d(:)];
                % SNR(ii,jj,iRef) = mean(temp) / std(temp);

            end
        end
    end

%% Attempt 1 v2024

% Spectrum Sample
Sp = zeros(m,n,p);
Sd = zeros(m,n,p);


% out = load([refDir,'/', refFiles(1).name]);
samRef = REF.rf;
samRef = samRef(ind_z,ind_x); % Cropping

    
f_MHz = axis_f((1:NFFT/2+1))*1e-6;
figure, 
% set(gcf, 'Position', [0 0 1 1], 'Units', 'normalized');
for jj=1:n
    for ii=1:m
        xw = x0(jj) ;   % x window
        zp = z0p(ii);
        zd = z0d(ii);

        sub_block_p = samRef(zp:zp+nz/2-1,xw:xw+nx-1);
        sub_block_d = samRef(zd:zd+nz/2-1,xw:xw+nx-1);

        [tempSp,~] = spectra(sub_block_p,windowing,0,nz/2,NFFT);
        [tempSd,~] = spectra(sub_block_d,windowing,0,nz/2,NFFT);

        % plot(f_MHz, 10*log10(tempSp((1:NFFT/2+1))/max(tempSp((1:NFFT/2+1)))),'k');
        % grid on; pause(0.1);
        xlim([0 14]), ylim([-70 0]);

        Sp(ii,jj,:) = tempSp(ind_f);
        Sd(ii,jj,:) = tempSd(ind_f);
    end
end

% % Step 2: Calculate Spectral Ratio
% RSp = Sp(:,:, 2:end) ./ Sp(:,:, 1:end-1);
% RSd = Sd(:,:, 2:end) ./ Sd(:,:, 1:end-1);
% RSp = cat(3, RSp(:, :, 1), RSp); % Prepend the first element
% RSd = cat(3, RSd(:, :, 1), RSd); % Prepend the first element 

RSp = zeros(m,n,p);
RSd = zeros(m,n,p);

% Compute the ratio along the third dimension
RSp(:, :, 2:end) = Sp(:, :, 2:end) ./ Sp(:, :, 1:end-1);
% For the last slice, keep the ratio the same as the previous slice
RSp(:, :, 1) = 1; % Assuming the first slice ratios are 1 as there is no "i-1"

% Compute the ratio along the third dimension
RSd(:, :, 2:end) = Sd(:, :, 2:end) ./ Sd(:, :, 1:end-1);
% For the last slice, keep the ratio the same as the previous slice
RSd(:, :, 1) = 1; % Assuming the first slice ratios are 1 as there is no "i-1"


b = log(RSp) - log(RSd);
b_vec = b(:);

df = fs/NFFT; % 0.26MHz

zd_zp = (nz/2)*dz;   % zd_zp = 2*\delta_z = nz/2 =  %  [m]
zp_zd = -zd_zp;

% System of eq
A1 = kron( -4*zp_zd*df*ones(p) , speye(m*n) );
A2 = kron( ones(p) , speye(m*n) );

A = [A1 A2];


%%
% axis_f
% 
idx_f_5M = find(band >= 5, 1, 'first');
idx_f_6M = find(band >= 6, 1, 'first');
idx_f_7M = find(band >= 7, 1, 'first');


% Choose half region 
z_half = ceil(m/2);
x_half = ceil(n/2);

FPR_5 = b(1:end, x_half, idx_f_5M);
FPR_6 = b(1:end, x_half, idx_f_6M);
FPR_7 = b(1:end, x_half, idx_f_7M);

figure, 
plot(z_ACS*1E2, FPR_5, 'DisplayName', '5MHz'), hold on, grid on;
plot(z_ACS*1E2, FPR_6, 'DisplayName', '6MHz'), hold on
plot(z_ACS*1E2, FPR_7, 'DisplayName', '7MHz'), hold on

xlabel('Depth [cm]')
ylabel('Freq Ratio [dB/MHz]')
legend('Location','Best')
set(gca, 'FontSize', 14)


%%
x_opt = cgs(A'*A, A'*b(:));

Bn = x_opt(1:m*n);
Cn = x_opt(m*n+1:end);

BR_RFM = reshape(Bn*Np2dB,m,n);
CR_RFM = reshape(Cn,m,n);

acs_map = reshape(BR_RFM*Np2dB,m,n);

figure, imagesc(acs_map), colorbar, colormap("jet");



%% Attempt 1 v2025

windowing = window_choice(nz/2, window_type);
windowing = windowing*ones(1,nx);

% z0 = z0p
z0 = 1:wz:length(z)-nz;
z_ACS = z(z0+nz/2);
m  = length(z0);

% Spectrum Sample
Sp = zeros(m,n,p);
Sd = zeros(m,n,p);


% out = load([refDir,'/', refFiles(1).name]);
samRef = SAM.rf;
samRef = samRef(ind_z,ind_x); % Cropping

    
f_MHz = axis_f((1:NFFT/2+1))*1e-6;
% figure, 
% set(gcf, 'Position', [0 0 1 1], 'Units', 'normalized');
font = 12;
for jj=1:n
    for ii=1:m
        xw = x0(jj) ;   % x window
        zw = z0(ii);

        sub_block = samRef(zw:zw+nz-1,xw:xw+nx-1);

        [tempSp,~] = spectra(sub_block,windowing,0,nz/2,NFFT);

        % plot(f_MHz, 10*log10(tempSp((1:NFFT/2+1))/max(tempSp((1:NFFT/2+1)))),'k');
        % grid on; pause(0.1);
        % xlim([0 14]), ylim([-70 0]);

        Sp(ii,jj,:) = tempSp(ind_f);

        if jj==floor(3*n/6)
                if ii==ceil(m/6)
                    figure(103)
                    set(gca,'FontSize',font);
                    % plot(f_MHz,10*log10(tempSp((1:NFFT/2+1))/max(tempSp((1:NFFT/2+1)))),'k');
                    plot(f_MHz,pow2db(tempSp((1:NFFT/2+1))),'k');
                   
                    title('Spectrum'); xlabel('\bfFrequency (MHz)'); ylabel('\bfIntensity Norm. (dB)');
                    axis([0 25 -70 0]);
                    grid on;
                    set(gca,'FontSize',font);
                    
                    % figure(106)
                    % spmax = max(tempSp((1:NFFT/2+1)));
                    % set(gca,'FontSize',font);                                
                    % plot(f_MHz,10*log10(tempSp((1:NFFT/2+1))/spmax),'k');
                    % %title('Spectrum'); 
                    % xlabel('\bfFrequency (MHz)'); ylabel('\bfIntensity Norm. (dB)');
                    % axis([0 25 -70 0]);
                    % set(gca,'FontSize',font); 
                    

                elseif ii==round(m/2)
                    figure(103)
                    hold on;
                    set(gca,'FontSize',font);
                    % plot(f_MHz,10*log10(tempSp((1:NFFT/2+1))/max(tempSp((1:NFFT/2+1)))),'r');
                    plot(f_MHz,pow2db(tempSp((1:NFFT/2+1))),'r');
                    title('Spectrum'); xlabel('\bfFrequency (MHz)'); ylabel('\bfIntensity Norm. (dB)');
                    axis([0 25 -70 0]);
                    set(gca,'FontSize',font);

                    % figure(106)
                    % hold on;
                    % set(gca,'FontSize',font);
                    % plot(f_MHz,10*log10(tempSp((1:NFFT/2+1))/spmax),'r');
                    % %title('Spectrum');
                    % xlabel('\bfFrequency (MHz)'); ylabel('\bfIntensity Norm. (dB)');
                    % axis([0 25 -70 0]);
                    % set(gca,'FontSize',font);
                                                    
                elseif ii==floor(4.5*m/6)
                    figure(103)
                    hold on;
                    set(gca,'FontSize',font);
                    % plot(f_MHz,10*log10(tempSp((1:NFFT/2+1))/max(tempSp((1:NFFT/2+1)))),'b');   %
                    plot(f_MHz,pow2db(tempSp((1:NFFT/2+1))),'b');
                    title('Spectrum'); xlabel('\bfFrequency (MHz)'); ylabel('\bfIntensity Norm. (dB)');
                    axis([0 25 -70 0]);
                    set(gca,'FontSize',font), 
                    grid on;
                    %pause
                    legend('Top','Half','Bottom');

                    % figure(106)
                    % hold on;
                    % set(gca,'FontSize',font);
                    % plot(f_MHz,10*log10(tempSp((1:NFFT/2+1))/spmax),'b');   %
                    % %title('Spectrum'); 
                    % xlabel('\bfFrequency (MHz)'); ylabel('\bfIntensity normalized (dB)');
                    % axis([0 25 -70 0]);
                    % set(gca,'FontSize',font); grid minor;
                    % %pause
                    % legend('1 cm above focal dist.','At focal distance','1 cm below focal dist.');
                    % 
                end
        end
    end
end

%% Replicate Fig. 1a Gong et al. 2019

% Choose half lateral  
x_half = ceil(n/2);

% Find freqs MHz
freq1 = 3.5; freq2 = 6; freq3 = 8.5;
idx_f1 = find(band >= freq1, 1, 'first');
idx_f2 = find(band >= freq2, 1, 'first');
idx_f3 = find(band >= freq3, 1, 'first');

figure, 
sgtitle('Spectrum 1a')
hold on, grid on;
% plot(z_ACS*1E2, log(Sp(:, x_half, idx_f1)), 'DisplayName', sprintf('S_{%.2fMHz}', band(idx_f1)) )
% plot(z_ACS*1E2, log(Sp(:, x_half, idx_f2)), 'DisplayName', sprintf('S_{%.2fMHz}', band(idx_f2)) )
% plot(z_ACS*1E2, log(Sp(:, x_half, idx_f3)), 'DisplayName', sprintf('S_{%.2fMHz}', band(idx_f3)) )

plot(z_ACS*1E2, pow2db(Sp(:, x_half, idx_f1)/max(Sp(:, x_half, idx_f1))), 'DisplayName', sprintf('S_{%.2fMHz}', band(idx_f1)) )
plot(z_ACS*1E2, pow2db(Sp(:, x_half, idx_f2)/max(Sp(:, x_half, idx_f2))), 'DisplayName', sprintf('S_{%.2fMHz}', band(idx_f2)) )
plot(z_ACS*1E2, pow2db(Sp(:, x_half, idx_f3)/max(Sp(:, x_half, idx_f3))), 'DisplayName', sprintf('S_{%.2fMHz}', band(idx_f3)) )

xlabel('Depth [cm]')
ylabel('Spectrum ln [dB]')
legend('Location','Best')
set(gca, 'FontSize', 14)

%
RSp = zeros(m,n,p);
RSp(:, :, 2:end) = Sp(:, :, 2:end) ./ Sp(:, :, 1:end-1);
% For the last slice, keep the ratio the same as the previous slice
RSp(:, :, 1) = RSp(:, :, 2); % Assuming the first slice ratios are 1 as there is no "i-1"

figure, 
sgtitle('Spectrum 1b')
plot(z_ACS*1E2, log(RSp(:, x_half, idx_f1)), 'DisplayName', sprintf('S_{%.2fMHz}/S_{%.2fMHz}', band(idx_f1), band(idx_f1-1)) ), hold on, grid on;
plot(z_ACS*1E2, log(RSp(:, x_half, idx_f2)), 'DisplayName', sprintf('S_{%.2fMHz}/S_{%.2fMHz}', band(idx_f2), band(idx_f2-1)) ), hold on
plot(z_ACS*1E2, log(RSp(:, x_half, idx_f3)), 'DisplayName', sprintf('S_{%.2fMHz}/S_{%.2fMHz}', band(idx_f3), band(idx_f3-1)) ), hold on

xlabel('Depth [cm]')
ylabel('Spectrum ln [dB]')
legend('Location','Best')
set(gca, 'FontSize', 14)

%%
% figure, 1b HALF LINE ONLY
df_MHz = min(diff(band));
Rs_div_f1 = log(RSp(:, x_half, idx_f1)) / (4*df_MHz/Np2dB);
Rs_div_f2 = log(RSp(:, x_half, idx_f2)) / (4*df_MHz/Np2dB);
Rs_div_f3 = log(RSp(:, x_half, idx_f3)) / (4*df_MHz/Np2dB);

figure, 
sgtitle('Spectrum 1b (MEDIUM LINE)')
plot(z_ACS*1E2, Rs_div_f1, 'DisplayName', sprintf('S_{%.2fMHz}/S_{%.2fMHz}', band(idx_f1), band(idx_f1-1)) ), hold on, grid on;
plot(z_ACS*1E2, Rs_div_f2, 'DisplayName', sprintf('S_{%.2fMHz}/S_{%.2fMHz}', band(idx_f2), band(idx_f2-1)) ), hold on
plot(z_ACS*1E2, Rs_div_f3, 'DisplayName', sprintf('S_{%.2fMHz}/S_{%.2fMHz}', band(idx_f3), band(idx_f3-1)) ), hold on

xlabel('Depth [cm]')
ylabel('Spectrum ln [dB/MHz]')
legend('Location','Best')
set(gca, 'FontSize', 14)

%%%%%%%%%%%%% ALL LINES %%%%%%%%%%%%%
df_MHz = min(diff(band));
Rs_div_f1 = mean ( log(RSp(:, :, idx_f1)) / (4*df_MHz/Np2dB) , 2);
Rs_div_f2 = mean ( log(RSp(:, :, idx_f2)) / (4*df_MHz/Np2dB) , 2);
Rs_div_f3 = mean ( log(RSp(:, :, idx_f3)) / (4*df_MHz/Np2dB) , 2);



figure, 
sgtitle('Spectrum 1b (ALL LINES)')
plot(z_ACS*1E2, Rs_div_f1, 'DisplayName', sprintf('S_{%.2fMHz}/S_{%.2fMHz}', band(idx_f1), band(idx_f1-1)) ), hold on, grid on;
plot(z_ACS*1E2, Rs_div_f2, 'DisplayName', sprintf('S_{%.2fMHz}/S_{%.2fMHz}', band(idx_f2), band(idx_f2-1)) ), hold on
plot(z_ACS*1E2, Rs_div_f3, 'DisplayName', sprintf('S_{%.2fMHz}/S_{%.2fMHz}', band(idx_f3), band(idx_f3-1)) ), hold on

xlabel('Depth [cm]')
ylabel('Spectrum ln [dB/MHz]')
legend('Location','Best')
set(gca, 'FontSize', 14)

%% Make full ACS all freq
z_ACS_cm = z_ACS*1E2;
df_MHz = min(diff(band));

ACS = zeros(p,n);
for ff = 1:p 
    for jj = 1:n
        SFD = log(squeeze(squeeze(RSp(:,jj,ff))));
        [lin_slope , lin_intercept, lin_y, R_fit] = fit_linear(z_ACS_cm, SFD, 2); 
        clear SFD
        ACS(ff, jj) = -lin_slope/(4*df_MHz) *Np2dB; % -8.68*m/(4*delta_f) [dB/MHz/cm]
    end
end

figure, 
imagesc(x_ACS*1E2, band, ACS), colorbar

xlabel('Lateral [cm]')
ylabel('Freq [MHz]')

%%
gt_acs = 0.4;
figure;
boxplot(ACS(:)); grid on
yline(gt_acs, 'k--')
ylim([-5 5])
title('Box Plot of ACS');
xlabel('ACS');
ylabel('Values');
%% HISTOGRAM

figure;
histogram(ACS(:,16), 20, 'Normalization', 'probability'); % 20 bins
xlabel('Value');
ylabel('Frequency');
title('Histogram of ACS Elements');
grid on;
%%
% % Step 2: Calculate Spectral Ratio
% RSp = Sp(:,:, 2:end) ./ Sp(:,:, 1:end-1);
% RSd = Sd(:,:, 2:end) ./ Sd(:,:, 1:end-1);
% RSp = cat(3, RSp(:, :, 1), RSp); % Prepend the first element
% RSd = cat(3, RSd(:, :, 1), RSd); % Prepend the first element 

RSp = zeros(m,n,p);
RSd = zeros(m,n,p);

% Compute the ratio along the third dimension
RSp(:, :, 2:end) = Sp(:, :, 2:end) ./ Sp(:, :, 1:end-1);
% For the last slice, keep the ratio the same as the previous slice
RSp(:, :, 1) = 1; % Assuming the first slice ratios are 1 as there is no "i-1"

% Compute the ratio along the third dimension
RSd(:, :, 2:end) = Sd(:, :, 2:end) ./ Sd(:, :, 1:end-1);
% For the last slice, keep the ratio the same as the previous slice
RSd(:, :, 1) = 1; % Assuming the first slice ratios are 1 as there is no "i-1"


b = log(RSp) - log(RSd);
b_vec = b(:);

df = fs/NFFT; % 0.26MHz

zd_zp = (nz/2)*dz;   % zd_zp = 2*\delta_z = nz/2 =  %  [m]
zp_zd = -zd_zp;

% System of eq
A1 = kron( -4*zp_zd*df*ones(p) , speye(m*n) );
A2 = kron( ones(p) , speye(m*n) );

A = [A1 A2];


