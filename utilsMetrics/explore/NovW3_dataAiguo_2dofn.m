% Test of calc_powerSpectra function (simplified Chahuara code)
% W3 Nov 2024 Data Aiguo 2DOF n

% clear all; 
% clc;
% close all;

Np2dB = 20*log10(exp(1));
dB2Np = 1/Np2dB;
range_bmode = [-60 0];
% Bmode = @(RF) 20*log10(abs(hilbert(RF))) - max (20*log10(abs(hilbert(RF))));
plotBSCdB = false;

methodEstim = 'normal';

reguType = 'high'; % normal, high, low

%%

pathDataUIUC = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\Data UIUC\';
pathEMZdata = [pathDataUIUC,'SiemensProcData_EMZ\'];

%%
operator = 'ZT'; % CHANGEABLE: AH or ZT
freqValue = 5.5; % 2.5, 4, 5.5 [MHZ]

numPhantomSam = 20;

numPhantomRef = 22;

% SAMPLE SPECS
if numPhantomSam == 20 % P20 ACS = 0.6851 dB/cm/MHz
    b_sam_the       = 5.705978E-06;
    n_sam_the       = 3.864606;       
    samPhan_alpha   = 0.6851; % [dB/cm/MHz]

elseif numPhantomSam == 22 % P22 ACS = 0.3007 dB/cm/MHz
    b_sam_the       = 5.125957E-07;
    n_sam_the       = 3.670544;
    samPhan_alpha   = 0.3007; % [dB/cm/MHz]
end

% REFERENCE SPECS
if numPhantomRef == 20 % P20 ACS = 0.6851 dB/cm/MHz    
    b_ref_the       = 5.705978E-06;
    n_ref_the       = 3.864606;    
    refPhan_alpha   = 0.6851; % [dB/cm/MHz]

elseif numPhantomRef == 22 % P22 ACS = 0.3007 dB/cm/MHz
    b_ref_the       = 5.125957E-07;
    n_ref_the       = 3.670544;
    refPhan_alpha   = 0.3007; % [dB/cm/MHz]
end

delta_n_the = n_sam_the - n_ref_the;
b_ratio_the = b_sam_the / b_ref_the;

delta_n     = delta_n_the;

% OLD 
% if numPhantomRef == 22  % P22 ACS = 0.3007 dB/cm/MHz
%     refPhan_alpha = 0.3007; % [dB/cm/MHz]
%     % delta_n       = 1.1;
%     delta_n       = 0.2; % rreal I think
% end
% 
% if numPhantomRef == 20  % P20 ACS = 0.6851 dB/cm/MHz
%     refPhan_alpha = 0.6851; % [dB/cm/MHz]
%     delta_n       = -1.1;
% end



phantomSamInfo = selectPhantom(operator, numPhantomSam, freqValue);
phantomRefInfo =  selectPhantom(operator, numPhantomRef, freqValue);

%%%%%%%%%%%%%%%%%%%%% RESUMEN %%%%%%%%%%%%%%%%%%%%%
% fprintf('Data resumen\n');

sam = {phantomSamInfo.operator, ['P', num2str(phantomSamInfo.numPhantom)], ...
    phantomSamInfo.freqChar, phantomSamInfo.numCase};

ref = {phantomRefInfo.operator, ['P', num2str(phantomRefInfo.numPhantom)], ...
    phantomRefInfo.freqChar, phantomRefInfo.numCase};

T = table(sam',ref', ...
    'VariableNames', {'SAMPLE','REFERENCE'}, ...
    'RowName',{'Operator','N° Phantom', 'Freq [MHz]', 'N° Case'}); 

disp(T)
clear sam ref

%%%%%%%%%%%%%%%%%%%%% FOLDER NAME %%%%%%%%%%%%%%%%%%%%%

folderNameSam =  ['P', num2str(phantomSamInfo.numPhantom), '_', ...
    'F', phantomSamInfo.freqChar, 'MHz_', ...
     num2str(phantomSamInfo.numCase)]

folderNameRef =  ['P', num2str(phantomRefInfo.numPhantom), '_', ...
    'F', phantomRefInfo.freqChar, 'MHz_', ...
     num2str(phantomRefInfo.numCase)]

%%%%%%%%%%%%%%%%%%%%% READ FOLDER %%%%%%%%%%%%%%%%%%%%%

pathData = [pathEMZdata, 'Operator_', phantomSamInfo.operator, ...
    '/', folderNameSam];
pathRef = [pathEMZdata, 'Operator_', phantomRefInfo.operator, ...
    '/', folderNameRef];

SAM = load([pathData, '/dataSiemens.mat']);
REF = load([pathRef, '/dataSiemens.mat']);

%% SET DATA
SAM.rf = SAM.RF_data(:,:,1);      
SAM.x  = SAM.lateral*1E-3; % was in [mm] set to [m]
SAM.z  = SAM.axial*1E-3;   % was in [mm] set to [m]
SAM.fs = 40E6; % [Hz]

REF.rf = REF.RF_data(:,:,1);      
REF.x  = REF.lateral*1E-3;  % was in [mm] set to [m]
REF.z  = SAM.axial*1E-3;  % was in [mm] set to [m]
REF.fs = 40E6; % [Hz]

% B-MODE CHECK
bmode_sam = db(hilbert(SAM.rf));
bmode_sam = bmode_sam - max(bmode_sam(:));

bmode_ref = db(hilbert(REF.rf));
bmode_ref = bmode_ref - max(bmode_ref(:));

% 
% figure,
% 
% subplot(121), 
% imagesc(SAM.x*1E2, SAM.z*1E2, bmode_sam, range_bmode), axis("image");
% xlabel('Lateral [cm]'), ylabel('Axial [cm]');
% title('SAM')
% colormap('gray')
% 
% subplot(122), 
% imagesc(REF.x*1E2, REF.z*1E2, bmode_ref,range_bmode), axis("image");
% xlabel('Lateral [cm]'), ylabel('Axial [cm]');
% title('REF')
% colormap('gray')

%% POWER LAW METHOD PARAMETERS

% pars.bw          = [1.6 4.7]; % [MHz]  % update depending the Power Spectra  original
% pars.bw          = [3 6]; % delta_n extra high
% pars.bw          = [1.9, 4]; % [MHz]
pars.bw          = [1.6 5.75]; % [MHz] % TUFFC
pars.bw          = [1.6 5.65];
pars.bw          = [1.6 5.5];
pars.overlap     = 0.8;
pars.blocksize   = 8; % wavelengths
pars.z_roi       = [3 10]*1E-2; % [m] 
pars.x_roi       = [-5 5]*1E-2; % [m] 

pars.saran_layer = false;

pars.window_type = 3; %  (1) Hanning, (2) Tuckey, (3) Hamming, (4) Tchebychev

%% REGULARIZATION PARAMETERS

% Implementation parameters
par_rpl.ini_method = 1; % METHOD LEAST SQUARES INITIALIZATION 
par_rpl.ini_tol = 1e-16;
par_rpl.df_op = 0; %(1) 1/2 f(x+1)-f(x-1), (0) f(x+1)-f(x)

par_rpl.tol   = 1e-16;
par_rpl.kmax  = 100;
par_rpl.eps_f = 1e-16;

% Robust
par_rpl.m_est   = 1;
par_rpl.c       = [1, 100, NaN, 100]; % c_fid, c_b, c_n, c_a
par_rpl.sigma   = [0.1, 0.1, NaN, 0.1]; % sigma_fid, sigma_b, sigma_n, sigma_a
% par_rpl.sigma   = [1.44, 1.44, 1.44, 1.44]; % sigma_fid, sigma_b, sigma_n, sigma_a

% Parameters for RPL-TV
%     mu_b  = 1E0-1E1; % ORIGINALS
%     mu_n  = 1E3-1E5; % ORIGINALS
%     mu_a  = 1E3-1E5; % ORIGINALS

% mu_rpl_tv    = [10^1.5, NaN, 1E4]; % [mu_b, mu_n, mu_a]
% mu_rpl_tv    = [10^1.75, NaN, 1E4]; % [mu_b, mu_n, mu_a]

if strcmp(reguType, 'high')
    mu_rpl_tv    = [1E3; 1E3; 1E5]; % [mu_b, mu_n, mu_a] % Test hyperregu
elseif strcmp(reguType, 'low')
    mu_rpl_tv    = [0.001; 0.001; 0.001]; % [mu_b, mu_n, mu_a] % Test no regu
elseif strcmp(reguType, 'normal')
    mu_rpl_tv    = [10^1.75, 1E4, 1E4]; % [mu_b, mu_n, mu_a] % Test regular    
end

%% POWER SPECTRA ESTIMATION

% spectralData_sam = calc_powerSpectra(SAM, pars);
spectralData_sam = calc_powerSpectra_vSimple(SAM, pars); % @

S_sam = spectralData_sam.powerSpectra;

num_ref = 1;
% spectralData_ref = calc_powerSpectra(REF, pars);
spectralData_ref = calc_powerSpectra_vSimple(REF, pars); % @
S_ref = spectralData_ref.powerSpectra;

% refFiles = dir([pathRef,'\*.mat']);
% num_ref = length(refFiles);
% 
% if num_ref > 1
%     S_ref_aux = zeros([size(S_sam), num_ref]);
% 
%     for iRef = 1:num_ref
%           REF = load(fullfile(pathRef, refFiles(iRef).name));
%         % spectralData_ref = calc_powerSpectra(REF, pars);
%         spectralData_ref = calc_powerSpectra_vSimple(REF, pars); % @
%         S_ref_aux(:,:,:,iRef) = spectralData_ref.powerSpectra;    
%     end
%     S_ref = mean(S_ref_aux, 4);
% end

%% SR and matrix creation (3rd dimension freq channels) Non Diagonal matrixes

% SR = S_sam ./ S_ref;
% 
% [r,p,q] = size(SR);
% 
% % Frequency vector (column vector of size [r, 1])
% band = spectralData_sam.band;
% f = band(:); % Ensure f is a column vector [r,1]
% 
% X = kron( ones(size(f)), speye(p*q) ); % Size: [p*q*r, p*q] 
% Z = kron( log(f) , speye(p*q) );       % Size: [p*q*r, p*q] 
% W = kron( -4*f, speye(p*q) );          % Size: [p*q*r, p*q] 
% 
% % 2DOF prior n
% delta_n = 0;   
% comp_ref = comp_ref_n_bsc_3rd(delta_n,band,p,q);
% 
% SR_comp = SR.*comp_ref;
% 
% % log-spectrum Ratio
% Y = log(SR_comp);

%% SR and matrix creation (1st dimension freq channels) Diagonal matrixes

SR_emz = S_sam ./ S_ref;

SR = permute(SR_emz,[3,1,2]);

[r,p,q] = size(SR);

% Frequency vector (column vector of size [r, 1])
band = spectralData_sam.band;
f = band(:); % Ensure f is a column vector [r,1]

X = kron( speye(p*q), ones(size(f)) ); % Size: [p*q*r, p*q] 
Z = kron( speye(p*q), log(f) );        % Size: [p*q*r, p*q] 
W = kron( speye(p*q), -4*f );          % Size: [p*q*r, p*q] 

%% 2DOF prior n COMPENSATION
 
comp_ref = comp_ref_n_bsc(delta_n,band,p,q);

SR_comp = SR.*comp_ref;

% log-spectrum Ratio Y = X.b + Z.n + W.a
Y = log(SR_comp);
%% WEIGHTS
% aSNR = 5; bSNR = 0.09;
% desvMin = 15;
% desvSNR = spectralData_sam.delta_snr;
% 
% wSNR = aSNR./(1 + exp(bSNR.*(desvSNR - desvMin)));

% low = 5e-10;
% high = 1;
% aux = desvSNR < 15;
% wSNR = high*aux + ~aux*low;
% 
% wSNR = ones(size(wSNR));
% 
% figure, 
% imagesc(wSNR), colorbar, title('Weights');

%% RPL

% initialization for RPL-based methods
u_0 = initialize_rpl_n_prior(Y,X,W,mu_rpl_tv,par_rpl);

if strcmp(methodEstim, 'normal')
    [u_opt,~] = rpl_tv_n_prior(Y,X,W,mu_rpl_tv,u_0,par_rpl);

elseif strcmp(methodEstim, 'robust')
    % Robust RPL estimation
    par_rpl.df_op = 0;
    [u_opt,~] = rpl_robtv_n_prior(Y,X,W,mu_rpl_tv,u_0,par_rpl);

elseif strcmp(methodEstim, 'weightedFid')
    par_rpl.weight = wSNR;
    [u_opt,~] = rpl_wFidtv_n_prior(Y,X,W,mu_rpl_tv,u_0,par_rpl);

elseif strcmp(methodEstim, 'weightedReg')
    par_rpl.weight = wSNR;
    [u_opt,~] = rpl_wRegtv_n_prior(Y,X,W,mu_rpl_tv,u_0,par_rpl);

elseif strcmp(methodEstim, 'weighted')
    par_rpl.weight = wSNR;
    [u_opt,~] = rpl_wtv_n_prior(Y,X,W,mu_rpl_tv,u_0,par_rpl);

elseif strcmp(methodEstim, 'rob-weighted')
    par_rpl.weight = wSNR;
    [u_opt,~] = rpl_robwtv_n_prior(Y,X,W,mu_rpl_tv,u_0,par_rpl);

elseif strcmp(methodEstim, 'spatial-admm')
    % computation
    tol = 1E-4;
    mask = ones(p,q,r);
    mu_b = 50; 
    mu_a = 5000;

    [u_opt,~] = ADMM_2dofn_WeightedTV(X,W,Y(:), mu_b, mu_a, ...
                    p, q, tol, mask(:), wSNR);
end
    
% A = [X, W];
% u_cgs = cgs(A'*A, A'*Y(:));

b = u_opt(1:p*q);
a = u_opt(p*q+1:2*p*q);

%% RESHAPE PARAMETERS

dy = 0.5*(diag(ones(p-1,1),1) - diag(ones(p-1,1),-1));
dy(1,1) = -1; dy(1,2) = 1; dy(end,end) = 1; dy(end,end-1) = -1;
Dy = sparse(kron(speye(q),dy));

z = 1E2*repmat(spectralData_sam.depth, 1, q);
% z = 1E2 * spectralData_sam.depth * ones(1, q); % 2d array
dz = reshape(Dy*z(:),p,q);
dz(end,:) = dz(end-1,:);

b_ratio     = reshape(exp(b),p,q);
b_ratio_dB  = 10*log10(b_ratio);
alpha_ratio = reshape(Np2dB*Dy*a./dz(:),p,q);
n_ratio     = delta_n*ones(size(b_ratio));

%%
%% METRICS
acs_sam = alpha_ratio + refPhan_alpha;

mean2d = @(x) mean(x(:));
std2d = @(x) std(x(:));
cv2d = @(x) 100*std(x(:))/mean(x(:));

m_n = mean2d(n_ratio);
s_n = std2d(n_ratio);
cv_n = cv2d(n_ratio);

m_b = mean2d(b_ratio);
s_b = std2d(b_ratio);
cv_b = cv2d(b_ratio);

m_a = mean2d(acs_sam);
s_a = std2d(acs_sam);
cv_a = cv2d(acs_sam);

disp(['\alpha : ', num2str(round(m_a, 3)), ' +/- ', ...
    num2str(round(s_a, 4)), ', %CV = ', num2str(round(cv_a, 4))]); 

disp(['\Delta b : ', num2str(round(m_b, 3)), ' +/- ', ...
    num2str(round(s_b, 4)), ', %CV = ', num2str(round(cv_b, 4))]); 

disp(['\Delta n : ', num2str(round(m_n, 4)), ' +/- ', ...
    num2str(round(s_n, 4)), ', %CV = ', num2str(round(cv_n, 4))]);
  
disp('--------')

metric_a = [m_a, s_a, cv_a];
metric_b = [m_b, s_b, cv_b];
metric_n = [m_n, s_n, cv_n];

% IMAGESC PLOTS
Xaxis = spectralData_ref.lateral;
Zaxis = spectralData_ref.depth;
cm = 1e2;

axis_n = [0 1.2];
axis_a = [0 3];
axis_b = [-60 0]; % dB
fontSize = 16;

figure, 
set(gcf,'units','normalized','outerposition',[0 0.25 1 0.65]); box on;
sgtitle('\bf2DOF n', 'FontSize', fontSize+2)
subplot(1,3,1)
% imagesc(Xaxis*cm, Zaxis*cm, Est_alphaz_sample, axis_a), colorbar
imagesc(Xaxis*cm, Zaxis*cm, acs_sam), colorbar
axis("image");
xlabel('Lateral [cm]'), ylabel('Depth [cm]'), colormap turbo;
% title('\Delta \alpha ');
title(['ACS : ', num2str(round(m_a, 3)), ' +/- ', num2str(round(s_a, 2)), ', %CV = ', num2str(round(cv_a, 3))])

h2 = colorbar; 
ylabel(h2,'dB\cdotcm^{-1}\cdotMHz^{-1}','FontSize', fontSize);
set(gca,'fontsize',fontSize)

subplot(1,3,2)
% imagesc(Xaxis*cm, Zaxis*cm, b_ratio_dB, axis_b), colorbar
% imagesc(Xaxis*cm, Zaxis*cm, b_ratio_dB), colorbar
imagesc(Xaxis*cm, Zaxis*cm, b_ratio)
h2 = colorbar; 
if plotBSCdB 
   imagesc(Xaxis*cm, Zaxis*cm, b_ratio_dB)
   h2 = colorbar;
   ylabel(h2,'dB','FontSize', fontSize);
end
axis("image");
xlabel('Lateral [cm]'), colormap turbo;
% title('\Delta b');
title(['\Delta b : ', num2str(round(m_b, 3)), ' +/- ', num2str(round(s_b, 2)), ', %CV = ', num2str(round(cv_b, 3))])
set(gca,'fontsize',fontSize)

subplot(1,3,3)
% imagesc(Xaxis*cm, Zaxis*cm, Est_n_sample, axis_n), colorbar
imagesc(Xaxis*cm, Zaxis*cm, n_ratio), colorbar
axis("image");
xlabel('Lateral [cm]'), colormap turbo;
% title('\Delta n  ');
title(['\Delta n : ', num2str(round(m_n, 3)), ' +/- ', num2str(round(s_n, 3)), ', %CV = ', num2str(round(cv_n, 3))])
set(gca,'fontsize',16)

%% FIGURE COLOR IMAGE DELTA ALPHA, DELTA B, DELTA N

% Xaxis = spectralData_ref.lateral;
% Zaxis = spectralData_ref.depth;
% cm = 1e2;
% 
% axis_n = [0 1.2];
% axis_a = [0 3];
% axis_b = [-60 0]; % dB
% fontSize = 16;
% 
% figure, 
% set(gcf,'units','normalized','outerposition',[0 0.25 1 0.65]); box on;
% 
% subplot(1,3,1)
% % imagesc(Xaxis*cm, Zaxis*cm, Est_alphaz_sample, axis_a), colorbar
% imagesc(Xaxis*cm, Zaxis*cm, alpha_ratio), colorbar
% axis("image");
% xlabel('Lateral [cm]'), ylabel('Depth [cm]'), colormap jet;
% title('\Delta \alpha ');
% h2 = colorbar; 
% ylabel(h2,'dB\cdotcm^{-1}\cdotMHz^{-1}','FontSize', fontSize);
% set(gca,'fontsize',fontSize)
% 
% subplot(1,3,2)
% % imagesc(Xaxis*cm, Zaxis*cm, b_ratio_dB, axis_b), colorbar
% % imagesc(Xaxis*cm, Zaxis*cm, b_ratio_dB), colorbar
% imagesc(Xaxis*cm, Zaxis*cm, b_ratio)
% h2 = colorbar; 
% if plotBSCdB 
%    imagesc(Xaxis*cm, Zaxis*cm, b_ratio_dB)
%    h2 = colorbar;
%    ylabel(h2,'dB','FontSize', fontSize);
% end
% axis("image");
% xlabel('Lateral [cm]'), colormap jet;
% title('\Delta b');
% set(gca,'fontsize',fontSize)
% 
% subplot(1,3,3)
% % imagesc(Xaxis*cm, Zaxis*cm, Est_n_sample, axis_n), colorbar
% imagesc(Xaxis*cm, Zaxis*cm, n_ratio), colorbar
% axis("image");
% xlabel('Lateral [cm]'), colormap jet;
% title('\Delta n  ');
% set(gca,'fontsize',16)

%% FIGURE INTERP OVERLAY BMODE, DELTA SNR, DELTA alpha, DELTA BSC, DELTA N

% fontSize = 14;
% 
% figure,
% set(gcf,'units','normalized','outerposition',[0 0 1 1])
% tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%% Bmode %%%%%%%%%%%%%%%%%%%%%%%%%%
% nexttile;
% axis off;
% 
% units           = 1E2;
% xFull           = SAM.x*units;
% zFull           = SAM.z*units;
% 
% tBmode = nexttile;
%     imagesc(xFull, zFull, bmode_sam, range_bmode), 
%     title(['Bmode ', T.SAMPLE{2}])
%     axis("image");
%     xlabel('Lateral'), ylabel('Depth');
%     c = colorbar;
%     c.Label.String = 'dB';
% %%%%%%%%%%%%%%%%%%%%%%%%%% Bmode %%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%% Delta SNR %%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% units           = 1E2;
% bmodeFull       = bmode_sam;
% colorImg        = spectralData_sam.delta_snr;
% range_bmode     = [-60 0];
% range_img       = [];
% transparency    = 0.65;
% x_img           = spectralData_sam.lateral*units;
% z_img           = spectralData_sam.depth*units;
% xFull           = SAM.x*units;
% zFull           = SAM.z*units;
% [X, Z] = meshgrid(xFull, zFull);
% roi = and(X >= x_img(1), X <= x_img(end)) & ...
%       and(Z >= z_img(1), Z <= z_img(end));
% 
% 
% t = nexttile;
% [~,hB,hColor] = imOverlayInterp(bmodeFull, colorImg, range_bmode, range_img, ...
%                     transparency, x_img, z_img, roi, xFull, zFull);
% hold on;
% contour(xFull, zFull, roi, 1,'r--', 'LineWidth', 2)
% hold off;
% xlabel('Lateral'), ylabel('Depth');
% hColor.Label.String = '';
% title('\Delta SNR')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%% Delta SNR %%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%% (ACS) %%%%%%%%%%%%%%%%%%%%%%%%%%
% acs_sam = alpha_ratio + refPhan_alpha;
% 
% units           = 1E2;
% bmodeFull       = bmode_sam;
% colorImg        = acs_sam;
% range_bmode     = [-60 0];
% range_img       = [];
% transparency    = 0.65;
% x_img           = spectralData_sam.lateral*units;
% z_img           = spectralData_sam.depth*units;
% xFull           = SAM.x*units;
% zFull           = SAM.z*units;
% [X, Z] = meshgrid(xFull, zFull);
% roi = and(X >= x_img(1), X <= x_img(end)) & ...
%       and(Z >= z_img(1), Z <= z_img(end));
% 
% t = nexttile;
%     [~,hB,hColor] = imOverlayInterp(bmodeFull, colorImg, range_bmode, range_img, ...
%                         transparency, x_img, z_img, roi, xFull, zFull);   
%     hold on;
%     contour(xFull, zFull, roi, 1,'r--', 'LineWidth', 2)
%     hold off;
%     xlabel('Lateral'), ylabel('Depth');
%     hColor.Label.String = 'dB\cdotcm^{-1}\cdotMHz^{-1}';
%     title('ACS')
% %%%%%%%%%%%%%%%%%%%%%%%%%% (ACS) %%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%% Delta alpha (BSC) %%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% units           = 1E2;
% bmodeFull       = bmode_sam;
% colorImg        = b_ratio;
% range_bmode     = [-60 0];
% range_img       = [];
% transparency    = 0.65;
% x_img           = spectralData_sam.lateral*units;
% z_img           = spectralData_sam.depth*units;
% xFull           = SAM.x*units;
% zFull           = SAM.z*units;
% [X, Z] = meshgrid(xFull, zFull);
% roi = and(X >= x_img(1), X <= x_img(end)) & ...
%       and(Z >= z_img(1), Z <= z_img(end));
% 
%     if plotBSCdB 
%        colorImg = b_ratio_dB;
%     end
% 
% t = nexttile;
% [~,hB,hColor] = imOverlayInterp(bmodeFull, colorImg, range_bmode, range_img, ...
%                     transparency, x_img, z_img, roi, xFull, zFull);
% hold on;
% contour(xFull, zFull, roi, 1,'r--', 'LineWidth', 2)
% hold off;
% xlabel('Lateral'), ylabel('Depth');
% hColor.Label.String = '';
%     if plotBSCdB 
%         hColor.Label.String ='dB';
%     end
% title('\Delta b')
% %%%%%%%%%%%%%%%%%%%%%%%%%% Delta alpha (BSC) in dB %%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%% Delta n %%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% units           = 1E2;
% bmodeFull       = bmode_sam;
% colorImg        = n_ratio;
% range_bmode     = [-60 0];
% range_img       = [];
% transparency    = 0.65;
% x_img           = spectralData_sam.lateral*units;
% z_img           = spectralData_sam.depth*units;
% xFull           = SAM.x*units;
% zFull           = SAM.z*units;
% [X, Z] = meshgrid(xFull, zFull);
% roi = and(X >= x_img(1), X <= x_img(end)) & ...
%       and(Z >= z_img(1), Z <= z_img(end));
% 
% 
% t = nexttile;
% [~,hB,hColor] = imOverlayInterp(bmodeFull, colorImg, range_bmode, range_img, ...
%                     transparency, x_img, z_img, roi, xFull, zFull);
% hold on;
% contour(xFull, zFull, roi, 1,'r--', 'LineWidth', 2)
% hold off;
% xlabel('Lateral'), ylabel('Depth');
% hColor.Label.String = '';
% title('\Delta n')
% %%%%%%%%%%%%%%%%%%%%%%%%%% Delta n %%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% % Apply font size to all axes labels, titles, and colorbars in one line
% colormap(tBmode,'gray')
% set(findall(gcf,'Type','text'),'FontSize',fontSize); % For all text elements
% set(findall(gcf,'Type','axes'),'FontSize',fontSize); % For axis labels and titles

%%

%% Litle HELP
% Delta_b*b_ref*f^(Delta_n +n_ref)

freq_bsc = spectralData_sam.band; % Given choosen BW

% OPTION A
% b_sam     = b_ratio * b_ref_the;
% n_sam     = n_ratio + n_ref_the;
% bsc_sam   = b_sam *(freq_bsc.^n_sam);
% med_bsc   = median(bsc_sam(:));

% OPTION B
b_sam     = median(b_ratio * b_ref_the, 'all');
n_sam     = median(n_ratio + n_ref_the, 'all');
med_bsc   = b_sam *(freq_bsc.^n_sam);

% THEORETICAL SAMPLE

%
fileNameBSC_P20 = 'BSC_vs_freq_P20.txt';
fileNameBSC_P22 = 'BSC_vs_freq_P22.txt';

pathBSC = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\Data UIUC\Siemens_S3000_6C1HD_THI_OFF';

% freqs | bsc

file_BSC_P20 = readmatrix(fullfile(pathBSC,fileNameBSC_P20)); % Freq, BSC
file_BSC_P22 = readmatrix(fullfile(pathBSC,fileNameBSC_P22)); % Freq, BSC

freqs_P20 = file_BSC_P20(:,1);
bsc_P20   = file_BSC_P20(:,2);

freqs_P22 = file_BSC_P22(:,1);
bsc_P22   = file_BSC_P22(:,2);


if numPhantomSam == 20
    % interp measurements for visualization (only see in choosen BW)
    BSC_gt = interp1(freqs_P20, bsc_P20, freq_bsc); 

elseif numPhantomSam == 22
    % interp measurements for visualization (only see in choosen BW)
    BSC_gt = interp1(freqs_P22, bsc_P22, freq_bsc); 
end

BSC_gt_dB  = 10*log10(BSC_gt);
med_bsc_dB = 10*log10(med_bsc);

abs_err = abs((med_bsc_dB - BSC_gt_dB)./BSC_gt_dB);

metrics.maedB       = mean(abs_err);
metrics.diff_fit_dB = mean ( abs(med_bsc_dB - BSC_gt_dB) ); 

fprintf('--------Metrics--------\n')
fprintf('Phantom       : P%d\n', numPhantomSam);
fprintf('Bandwidth (BW): %.2f - %.2f MHz\n', pars.bw(1), pars.bw(2));
fprintf('MAE dB        : %.4f\n', metrics.maedB);
fprintf('Diff Abs dB   : %.4f\n', metrics.diff_fit_dB)

% BSC RESULTS
linewidth = 2;

figure,
% set(gcf,'units','normalized','outerposition',[0 0 1 1]);
semilogy(band, med_bsc, '-', 'LineWidth', linewidth+1), hold on
semilogy(band, BSC_gt, 'k-.', 'LineWidth', linewidth), grid minor;

% Add a textbox to display GoF_dB value
gof_text = ['GoF_{dB} = ', num2str(metrics.diff_fit_dB, '%.2f')];
annotation('textbox', [0.15, 0.7, 0.3, 0.1], 'String', gof_text, ...
    'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', ...
    'FontSize', 10);

xlabel('Frequency [MHz]'), 
ylabel('BSC [cm^{-1}\cdot sr^{-1}]')


title(['BSC P', num2str(numPhantomSam)])
legend ('BSC', 'Ground truth','Location', 'Best')
fontsize(gcf, 10, 'points')


% %% PLOT 3DOF and 2DOF n same time
% 
% linewidth = 2;
% 
% figure,
% % set(gcf,'units','normalized','outerposition',[0 0 1 1]);
% loglog(band, med_bsc_2dof, '-', 'LineWidth', linewidth+1, 'DisplayName', '2-DoF'), hold on
% loglog(band, med_bsc_3dof, '-', 'LineWidth', linewidth+1, 'DisplayName', '3-DoF'), hold on
% loglog(band, BSC_gt, 'k-.', 'LineWidth', linewidth, 'DisplayName','Ground truth'), grid minor;
% 
% 
% xlabel('Frequency [MHz]'), 
% ylabel('BSC [cm^{-1}\cdot sr^{-1}]')
% 
% 
% title(['BSC P', num2str(numPhantomSam)])
% legend ('Location', 'Best')
% fontsize(gcf, 10, 'points')