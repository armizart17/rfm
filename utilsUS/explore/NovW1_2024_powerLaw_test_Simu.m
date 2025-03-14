% Test of calc_powerSpectra function (simplified Chahuara code)
% Nov 2024 EMZ

%%
clear all, clc, close all;

addpath(genpath(pwd))
%% LOAD DATA
pathData = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\dataACS_kwave';

SAM = load(fullfile(pathData, 'sam0.mat'));
REF = load(fullfile(pathData, 'ref0.mat'));

refPhan_alpha = 0.4; % [dB/cm/MHz]

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


%% POWER LAW METHOD PARAMETERS

pars.bw          = [3 9]; % [MHz]
pars.overlap     = 0.8;
pars.blocksize   = 20; % wavelengths
pars.z_roi       = [5 45]*1E-3; % [m] 
pars.x_roi       = [-15 15]*1E-3; % [m] 
pars.window_type = 3; %  (1) Hanning, (2) Tuckey, (3) Hamming, (4) Tchebychev
pars.saran_layer = true;

%% REGULARIZATION PARAMETERS

% Implementation parameters
par_rpl.tol   = 1e-16;
par_rpl.kmax  = 100;
par_rpl.eps_f = 1e-16;
par_rpl.m_est = 0; %Robust

% Parameters for RPL-TV
%     mu_b  = 1E0-1E1; % ORIGINALS
%     mu_n  = 1E3-1E5; % ORIGINALS
%     mu_a  = 1E3-1E5; % ORIGINALS

mu_rpl_tv    = [1E3; 1E3; 1E5]; % [mu_b, mu_n, mu_a]
par_rpl.ini_method = 1; % METHOD LEAST SQUARES INITIALIZATION 
par_rpl.ini_tol = 1e-16;
par_rpl.df_op = 1;


%% POWER SPECTRA ESTIMATION
% spectralData_sam = calc_powerSpectra(SAM, pars);
spectralData_sam = calc_powerSpectra_vSimple(SAM, pars);
S_sam = spectralData_sam.powerSpectra;

num_ref = 1;

% spectralData_ref = calc_powerSpectra(REF, pars);
spectralData_ref = calc_powerSpectra_vSimple(REF, pars);
S_ref = spectralData_ref.powerSpectra;


% if num_ref > 1
%     S_ref_aux = zeros([size(S_sam), num_ref]);
% 
%     for iRef = 1:num_ref
%         spectralData_ref = calc_powerSpectra(REF, pars);
%         S_ref_aux(:,:,:,iRef) = spectralData_ref.powerSpectra;    
%     end
%     S_ref = mean(S_ref_aux, 4);
% end



%% SR CALCULATION AND PERMUTATION
SR_emz = S_sam ./ S_ref;

SR = permute(SR_emz,[3,1,2]);

%% 2DOF prior n
delta_n = 0;
band = spectralData_sam.band;
[r,p,q] = size(SR);
    
comp_ref = comp_ref_n_bsc(delta_n,band,p,q);

SR_comp = SR.*comp_ref;

% indices initialization
f = band(:); % always column vector

% log-spectrum Ratio
Y = log(SR_comp);

% matrices for RPL-based algorithms
X = kron(speye(p*q), ones(r,1));
W = -4*kron(speye(p*q), f);

% initialization for RPL-based methods
u_0 = initialize_rpl_n_prior(Y,X,W,mu_rpl_tv,par_rpl);
    
% RPL estimation
[u_opt,~] = rpl_tv_n_prior(Y,X,W,mu_rpl_tv,u_0,par_rpl);

% A = [X, W];
% u_cgs = cgs(A'*A, A'*Y(:));

b = u_opt(1:p*q);
a = u_opt(p*q+1:2*p*q);

%% HELP

% check later 02 Nov
   % if df_op == 1
   %      dy = 0.5*(diag(ones(p-1,1),1) - diag(ones(p-1,1),-1));
   %      dy(1,1) = -1; dy(1,2) = 1; dy(end,end) = 1; dy(end,end-1) = -1;
   %  else
   %      dy = diag(ones(p-1,1),1) - diag([ones(p-1,1);0]);
   %  end
   % 
   %  Dy = sparse(kron(speye(q),dy));
   % 
   %  z = 1e+2*repmat(depth,1,q);
   %  dz = reshape(Dy*z(:),p,q);
   %  dz(end,:) = dz(end-1,:);
   % 
   %  b = u(1:p*q); n = u(1+p*q:2*p*q); a = u(1+2*p*q:end);
   % 
   %  b_ref = refvals{1}; n_ref = refvals{2}; a_ref = refvals{3};
   % 
   %  est_maps = cell(3+2*(src == 3),1);
   % 
   %  %b_param = exp(b)*b_ref;
   %  %n_param = n + n_ref;
   %  %a_param = 8.686*Dy*a./dz(:) + a_ref;
   %  b_param = exp(b);
   %  n_param = n;
   %  a_param = 8.686*Dy*a./dz(:);

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

Est_b_sample = b_ratio;
Est_n_sample = n_ratio;
Est_alphaz_sample = alpha_ratio + refPhan_alpha;

% Est_b_sample = 10*log10(Est_b_sample);

mean2d = @(x) mean(x(:));
std2d = @(x) std(x(:));
cv2d = @(x) 100*std(x(:))/mean(x(:));

m_n = mean2d(Est_n_sample);
s_n = std2d(Est_n_sample);
cv_n = cv2d(Est_n_sample);

m_b = mean2d(Est_b_sample);
s_b = std2d(Est_b_sample);
cv_b = cv2d(Est_b_sample);

m_a = mean2d(Est_alphaz_sample);
s_a = std2d(Est_alphaz_sample);
cv_a = cv2d(Est_alphaz_sample);


metric_a = [m_a, s_a, cv_a];
metric_b = [m_b, s_b, cv_b];
metric_n = [m_n, s_n, cv_n];
%%

Xaxis = spectralData_ref.lateral;
Zaxis = spectralData_ref.depth;
cm = 1e2;

axis_n = [0 1.2];

axis_a = [0 3];

axis_b = [-60 0]; % dB


figure, 
set(gcf,'units','normalized','outerposition',[0 0.25 1 0.65]); box on;

subplot(1,3,1)
% imagesc(Xaxis*cm, Zaxis*cm, Est_n_sample, axis_n), colorbar
imagesc(Xaxis*cm, Zaxis*cm, Est_n_sample), colorbar
xlabel('Lateral [cm]'), ylabel('Depth [cm]'), colormap jet;
title(['\Delta n : ', num2str(round(m_n, 3)), ' +/- ', num2str(round(s_n, 3)), ', %CV = ', num2str(round(cv_n, 3))])
set(gca,'fontsize',16)
subplot(1,3,2)
% imagesc(Xaxis*cm, Zaxis*cm, Est_b_sample, axis_b), colorbar
imagesc(Xaxis*cm, Zaxis*cm, Est_b_sample), colorbar
xlabel('Lateral [cm]'), ylabel('Depth [cm]'), colormap jet;
title(['\Delta b : ', num2str(round(m_b, 3)), ' +/- ', num2str(round(s_b, 2)), ', %CV = ', num2str(round(cv_b, 3))])
set(gca,'fontsize',16)

subplot(1,3,3)
% imagesc(Xaxis*cm, Zaxis*cm, Est_alphaz_sample, axis_a), colorbar
imagesc(Xaxis*cm, Zaxis*cm, Est_alphaz_sample), colorbar
xlabel('Lateral [cm]'), ylabel('Depth [cm]'), colormap jet;
title(['\Delta \alpha : ', num2str(round(m_a, 3)), ' +/- ', num2str(round(s_a, 2)), ', %CV = ', num2str(round(cv_a, 3))])
set(gca,'fontsize',16)


%% INTERP OVERLAY BMODE + COLORIMAGE EXAMPLE

units = 1E3;

bmodeFull       = bmode_sam;
colorImg        = Est_alphaz_sample;
range_bmode     = [-60 0];
range_img       = [0 1];
transparency    = 0.7;
x_img           = spectralData_sam.lateral*units;
z_img           = spectralData_sam.depth*units;
xFull           = SAM.x*units;
zFull           = SAM.z*units;

[X, Z] = meshgrid(xFull, zFull);
roi = and(X >= x_img(1), X <= x_img(end) ) & and( Z >= z_img(1), Z <= z_img(end));


figure,
% subplot(121),

imagesc(xFull, zFull, bmode_sam), axis("image");
xlabel('Lateral'), ylabel('Axial');
title('Bmode')
colormap('gray')

figure, 
% subplot(122),

[~,hB,hColor] = imOverlayInterp(bmodeFull, colorImg, range_bmode, range_img, ...
                    transparency, x_img, z_img, roi, xFull, zFull);

hold on;
contour(xFull, zFull, roi, 1,'r--', 'LineWidth', 2)
hold off;


xlabel('Lateral'), ylabel('Axial');
hColor.Label.String = '[dB\cdotcm^{-1}\cdotMHz^{-1}]';

title('ColorImage+Bmode')
% 
%%
%% FIGURE INTERP OVERLAY BMODE, DELTA SNR, ACS, DELTA BSC, DELTA N

fontSize = 14;
plotBSCdB = false;

figure,
tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

%%%%%%%%%%%%%%%%%%%%%%%%%% Bmode %%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile;
axis off;

units           = 1E3;
xFull           = SAM.x*units;
zFull           = SAM.z*units;

tBmode = nexttile;
    imagesc(xFull, zFull, bmode_sam, range_bmode), 
    title('Bmode')
    axis("image");
    xlabel('Lateral'), ylabel('Depth');
%%%%%%%%%%%%%%%%%%%%%%%%%% Bmode %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% Delta SNR %%%%%%%%%%%%%%%%%%%%%%%%%%

units           = 1E3;
bmodeFull       = bmode_sam;
colorImg        = spectralData_sam.delta_snr;
range_bmode     = [-60 0];
range_img       = [];
transparency    = 0.65;
x_img           = spectralData_sam.lateral*units;
z_img           = spectralData_sam.depth*units;
xFull           = SAM.x*units;
zFull           = SAM.z*units;
[X, Z] = meshgrid(xFull, zFull);
roi = and(X >= x_img(1), X <= x_img(end)) & ...
      and(Z >= z_img(1), Z <= z_img(end));


t = nexttile;
[~,hB,hColor] = imOverlayInterp(bmodeFull, colorImg, range_bmode, range_img, ...
                    transparency, x_img, z_img, roi, xFull, zFull);
hold on;
contour(xFull, zFull, roi, 1,'r--', 'LineWidth', 2)
hold off;
xlabel('Lateral'), ylabel('Depth');
hColor.Label.String = '';
title('\Delta SNR')

%%%%%%%%%%%%%%%%%%%%%%%%%% Delta SNR %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% Delta alpha (ACS) %%%%%%%%%%%%%%%%%%%%%%%%%%
acs_sam = alpha_ratio + refrefPhan_alpha;

units           = 1E3;
bmodeFull       = bmode_sam;
colorImg        = acs_sam;
range_bmode     = [-60 0];
range_img       = [0 0.85];
transparency    = 0.65;
x_img           = spectralData_sam.lateral*units;
z_img           = spectralData_sam.depth*units;
xFull           = SAM.x*units;
zFull           = SAM.z*units;
[X, Z] = meshgrid(xFull, zFull);
roi = and(X >= x_img(1), X <= x_img(end)) & ...
      and(Z >= z_img(1), Z <= z_img(end));

t = nexttile;
    [~,hB,hColor] = imOverlayInterp(bmodeFull, colorImg, range_bmode, range_img, ...
                        transparency, x_img, z_img, roi, xFull, zFull);   
    hold on;
    contour(xFull, zFull, roi, 1,'r--', 'LineWidth', 2)
    hold off;
    xlabel('Lateral'), ylabel('Depth');
    hColor.Label.String = 'dB\cdotcm^{-1}\cdotMHz^{-1}';
    title('ACS Sam')
%%%%%%%%%%%%%%%%%%%%%%%%%% Delta alpha (ACS) %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% Delta alpha (BSC) %%%%%%%%%%%%%%%%%%%%%%%%%%

units           = 1E3;
bmodeFull       = bmode_sam;
colorImg        = b_ratio;
range_bmode     = [-60 0];
range_img       = [];
transparency    = 0.65;
x_img           = spectralData_sam.lateral*units;
z_img           = spectralData_sam.depth*units;
xFull           = SAM.x*units;
zFull           = SAM.z*units;
[X, Z] = meshgrid(xFull, zFull);
roi = and(X >= x_img(1), X <= x_img(end)) & ...
      and(Z >= z_img(1), Z <= z_img(end));

    if plotBSCdB 
       colorImg = b_ratio_dB;
    end

t = nexttile;
[~,hB,hColor] = imOverlayInterp(bmodeFull, colorImg, range_bmode, range_img, ...
                    transparency, x_img, z_img, roi, xFull, zFull);
hold on;
contour(xFull, zFull, roi, 1,'r--', 'LineWidth', 2)
hold off;
xlabel('Lateral'), ylabel('Depth');
hColor.Label.String = '';
    if plotBSCdB 
        hColor.Label.String ='dB';
    end
title('\Delta BSC')
%%%%%%%%%%%%%%%%%%%%%%%%%% Delta alpha (BSC) in dB %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% Delta n %%%%%%%%%%%%%%%%%%%%%%%%%%

units           = 1E3;
bmodeFull       = bmode_sam;
colorImg        = n_ratio;
range_bmode     = [-60 0];
range_img       = [];
transparency    = 0.65;
x_img           = spectralData_sam.lateral*units;
z_img           = spectralData_sam.depth*units;
xFull           = SAM.x*units;
zFull           = SAM.z*units;
[X, Z] = meshgrid(xFull, zFull);
roi = and(X >= x_img(1), X <= x_img(end)) & ...
      and(Z >= z_img(1), Z <= z_img(end));


t = nexttile;
[~,hB,hColor] = imOverlayInterp(bmodeFull, colorImg, range_bmode, range_img, ...
                    transparency, x_img, z_img, roi, xFull, zFull);
hold on;
contour(xFull, zFull, roi, 1,'r--', 'LineWidth', 2)
hold off;
xlabel('Lateral'), ylabel('Depth');
hColor.Label.String = '';
title('\Delta n')
%%%%%%%%%%%%%%%%%%%%%%%%%% Delta n %%%%%%%%%%%%%%%%%%%%%%%%%%


% Apply font size to all axes labels, titles, and colorbars in one line
colormap(tBmode,'gray')
set(findall(gcf,'Type','text'),'FontSize',fontSize); % For all text elements
set(findall(gcf,'Type','axes'),'FontSize',fontSize); % For axis labels and titles



