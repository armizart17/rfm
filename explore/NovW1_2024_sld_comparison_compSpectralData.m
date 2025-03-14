% Test of calc_power spectra function (simplified Chahuara code)
% Nov 2024 EMZ

%%
clear all, clc, close all;
Np2dB = 20*log(exp(1));
dB2Np = 1/Np2dB;
addpath(genpath(pwd));
%% LOAD DATA
pathData = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\dataACS_kwave';

SAM = load(fullfile(pathData, 'sam0.mat'));
REF = load(fullfile(pathData, 'ref0.mat'));

% B-MODE CHECK
bmode_sam = db(hilbert(SAM.rf));
bmode_sam = bmode_sam - max(bmode_sam(:));

bmode_ref = db(hilbert(REF.rf));
bmode_ref = bmode_ref - max(bmode_ref(:));

figure,

subplot(121), 
imagesc(SAM.x*1E3, SAM.z*1E3, bmode_sam), axis("image");
xlabel('Lateral [mm]'), ylabel('Depth [mm]');
title('SAM')
colormap('gray')

subplot(122), 
imagesc(REF.x*1E3, REF.z*1E3, bmode_ref), axis("image");
xlabel('Lateral [mm]'), ylabel('Depth [mm]');
title('REF')
colormap('gray')


%% SLD PARAMETERS METHOD AC

clear parsSLD;

parsSLD.z_roi       = [5 45]*1E-3; % [m] 
parsSLD.x_roi       = [-17 17]*1E-3; % [m] 

parsSLD.cs = 1540; % speed of sound SoS in [m/s]
parsSLD.P = 2^10; % (number of points FFT)
parsSLD.bw = [3 9];   % BANDWITH i.e [ 4 10 ] % in [MHz]


% LINEAR DEPENDENCY
parsSLD.REF_num = 111; % (see function attenuation_phantoms_Np)
parsSLD.SLOPE = 0.4;

parsSLD.blocksize = 20; % datablock in wavelengths
parsSLD.overlap = 0.8;

parsSLD.window_type = 3;

REF.rf1 = REF.rf;
REF.rf2 = REF.rf;
REF.rf3 = REF.rf;
REF.rf4 = REF.rf;

% vCoila
% SLD_AC = SLD_function(SAM, REF, parsSLD);
pars.P = 1024;
parsSLD.nb_lambda_lateral = 20; 
% vRouyer
SLD_AC = SLD_function_v2(SAM, REF, parsSLD);
SLogRatio_AC = SLD_AC.SLD_term;



%% SLD MA METHOD PARAMETERS

% PARAMETERS
pars.bw          = [3 9]; % [MHz]
pars.overlap     = 0.8;
pars.blocksize   = 20; % wavelengths
pars.z_roi       = [0 45]*1E-3; % [m] 
pars.x_roi       = [-17 17]*1E-3; % [m] 
pars.window_type = 3; %  (1) Hanning, (2) Tuckey, (3) Hamming, (4) Tchebychev
pars.saran_layer = true;

% SAMPLE
spectralData_sam = calc_powerSpectra_prox_dis(SAM, pars);

Sp_sam = spectralData_sam.Sp;
Sd_sam = spectralData_sam.Sd;

% REFERENCE
num_ref = 1;
spectralData_ref = calc_powerSpectra_prox_dis(REF, pars);

band   = spectralData_sam.band; % [MHz]
zd_zp  = spectralData_sam.zd_zp; % [m]
Sp_ref = spectralData_ref.Sp;
Sd_ref = spectralData_ref.Sd;

% if num_ref > 1
%     S_ref_aux1 = zeros([size(S_sam), num_ref]);
%     S_ref_aux1 = zeros([size(S_sam), num_ref]);
% 
%     for iRef = 1:num_ref
%         spectralData_ref = calc_power_spectra(REF, pars);
%         S_ref_aux1(:,:,:,iRef) = spectralData_ref.Sp; 
%         S_ref_aux2(:,:,:,iRef) = spectralData_ref.Sd;     
%     end
%     Sp_ref = mean(S_ref_aux1, 4); clear(S_ref_aux1);
%     Sd_ref = mean(S_ref_aux2, 4); clear(S_ref_aux2);
% end

% COMPENSATION ATTENUATION
acs_ref     = 0.4; % [dB/cm/MHz]
att_ref     = acs_ref*band/Np2dB; % vector [Np/cm]
att_ref_map = ones(size(Sp_sam)) .* reshape(att_ref, 1, 1, []); % 3D array as SLogRatio

SLogRatio = log(Sp_sam./Sd_sam) - ( log(Sp_ref./Sd_ref) - 4*zd_zp*1E2*att_ref_map ); % (rows, cols, freqChannels)
clear('Sp_sam','Sd_sam', 'Sp_ref', 'Sd_ref', 'att_ref_map', 'att_ref');

%% OPTIMIZATION

% Optimization constants
tol = 1e-3;
mu1 = 10^4.5; 
mu2 = 10^4.5;

% V1 
[m, n, p] = size(SLogRatio_AC);
mask = ones(m,n,p);

[Bn,~] = AlterOpti_ADMM(SLD_AC.A1,SLD_AC.A2,SLogRatio_AC(:),mu1,mu2,m,n,tol,mask(:));
ACS_v1 = reshape(Bn*Np2dB,m,n);

bestTau = 0.001;
mu_wtnv = 2.7696; % best JASA
tol = 0.5e-4; % tolerance error
stableIter = 200;
weightEstimators = rescale(SLD_AC.muX, 1, max(SLD_AC.muX));
maxIter = 1000;
SLD_tnv = pdo_den_wtnv(SLogRatio_AC, mu_wtnv, bestTau, maxIter, tol, stableIter, weightEstimators);       
[ACS_v1, n_tnv] = cgs_ACS(SLD_AC.A, SLD_tnv);

%
bestTau = 0.001;
mu_wtnv = 2.7696; % best JASA
tol = 0.5e-4; % tolerance error
stableIter = 200;
weightEstimators = rescale(SLD_AC.muX, 1, max(SLD_AC.muX));
maxIter = 1000;
SLD_tnv = pdo_den_wtnv(SLogRatio, mu_wtnv, bestTau, maxIter, tol, stableIter, weightEstimators);       
[ACS_v1, n_tnv] = cgs_ACS([A1 A2], SLD_tnv);


% V2
% System of eq
mu1 = 10^3.5;
mu2 = 10^3.5;
[m, n, p] = size(SLogRatio);
mask = ones(m,n,p);
A1 = kron( 4*zd_zp*1E2*band , speye(m*n) );
A2 = kron( ones(size(band)) , speye(m*n) );
tol = 1e-3;
[Bn,~] = AlterOpti_ADMM(A1,A2,SLogRatio(:),mu1,mu2,m,n,tol,mask(:));
ACS_v2 = reshape(Bn*Np2dB,m,n);

%% Display 
figure, 
subplot(121), 
imagesc(ACS_v1), colormap("jet"), colorbar
title('ACS v1')

subplot(122), 
imagesc(ACS_v2, [0 1.1]), colormap("jet"), colorbar
title('ACS v2')

%%
%% INTERP OVERLAY BMODE + COLORIMAGE EXAMPLE

units = 1E3;

% V1
bmodeFull       = bmode_sam;
colorImg        = ACS_v1;
range_bmode     = [-60 0];
range_img       = [0 3];
transparency    = 0.7;
x_img           = SLD_AC.x_v2*units;
z_img           = SLD_AC.z_v2*units;
xFull           = SAM.x*units;
zFull           = SAM.z*units;

[X, Z] = meshgrid(xFull, zFull);
roi = and(X >= x_img(1), X <= x_img(end) ) & and( Z >= z_img(1), Z <= z_img(end));

figure, 

[~,hB,hColor] = imOverlayInterp(bmodeFull, colorImg, range_bmode, range_img, ...
                    transparency, x_img, z_img, roi, xFull, zFull);

hold on;
contour(xFull, zFull, roi, 1,'r--', 'LineWidth', 2)
hold off;


xlabel('Lateral'), ylabel('Axial');
hColor.Label.String = 'dB\cdotcm^{-1}\cdotMHz^{-1}';

title('ACS+Bmode v1')

%% V2
bmodeFull       = bmode_sam;
colorImg        = ACS_v2;
range_bmode     = [-60 0];
range_img       = [0 2];
transparency    = 0.7;
x_img           = spectralData_sam.lateral*units;
z_img           = spectralData_sam.depth*units;
xFull           = SAM.x*units;
zFull           = SAM.z*units;

[X, Z] = meshgrid(xFull, zFull);
roi = and(X >= x_img(1), X <= x_img(end) ) & and( Z >= z_img(1), Z <= z_img(end));


figure, 

[~,hB,hColor] = imOverlayInterp(bmodeFull, colorImg, range_bmode, range_img, ...
                    transparency, x_img, z_img, roi, xFull, zFull);

hold on;
contour(xFull, zFull, roi, 1,'r--', 'LineWidth', 2)
hold off;


xlabel('Lateral'), ylabel('Axial');
hColor.Label.String = 'dB\cdotcm^{-1}\cdotMHz^{-1}';

title('ACS+Bmode v2')