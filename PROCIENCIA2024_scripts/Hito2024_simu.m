% Test of calc_powerSpectra function (simplified Chahuara code)
% Dec 2024 Data Aiguo 2DoF n
% By EAMZ

% clear all; 
clc;
warning('off');

Np2dB = 20*log10(exp(1));
dB2Np = 1/Np2dB;
range_bmode = [-60 0];
% Bmode = @(RF) 20*log10(abs(hilbert(RF))) - max (20*log10(abs(hilbert(RF))));
plotBSCdB = true;

% MODEL a b y n
% (1) case (variation a)
% SAM
% a0.6 b0.01 n1.5
% REF
% a0.4 b0.01 n1.5

% (2) case (variation a , b)
% SAM
% a0.6 b0.01 n1.5
% REF
% a0.4 b1 n1.5

% (2) case (variation a , b , n)
% SAM
% a0.6 b0.01 n1.5
% REF
% a0.4 b1 n0

%% PATH DATA
pathData = '.\NewPL_data\data_ok\';

%%

a_list = [0.4, 0.6];
b_list = [0.01 1];
n_list = [0 1.5];

numCase = 3;
switch numCase
    case 1 % ONLY VARIATION a
        a_sam_the = 0.6; a_ref_the = 0.4;
        b_sam_the = 0.01; b_ref_the = 0.01;
        n_sam_the = 1.5; n_ref_the = 1.5;
    case 2  % ONLY VARIATION a, b
        a_sam_the = 0.6; a_ref_the = 0.4;
        b_sam_the = 0.01; b_ref_the = 1;
        n_sam_the = 1.5; n_ref_the = 1.5;
    case 3  % ONLY VARIATION a, b, n
        a_sam_the = 0.6; a_ref_the = 0.4;
        b_sam_the = 0.01; b_ref_the = 1;
        n_sam_the = 1.5; n_ref_the = 0;
end

AUX = load([pathData,'ref0.mat']); % JUST AUXILIAR TO GIVE AXIS
name_sam = ['a',num2str(a_sam_the),'b',num2str(b_sam_the),'n',num2str(n_sam_the)];
name_ref = ['a',num2str(a_ref_the),'b',num2str(b_ref_the),'n',num2str(n_ref_the)];

SAM     = load([pathData, 'ref_', name_sam]);
SAM.rf  = SAM.rfnew;
SAM.x   = AUX.x;
SAM.z   = AUX.z;
SAM.fs  = AUX.fs;

REF     = load([pathData, 'ref_', name_ref]);
REF.rf  = REF.rfnew;
REF.x   = AUX.x;
REF.z   = AUX.z;
REF.fs  = AUX.fs;

clear AUX

delta_n_the = n_sam_the - n_ref_the;
b_ratio_the = b_sam_the / b_ref_the;
alpha_ref_the = a_ref_the;

delta_n     = delta_n_the;

%% SET DATA


% B-MODE CHECK
bmode_sam = db(abs(hilbert(SAM.rf)));
bmode_sam = bmode_sam - max(bmode_sam(:));

bmode_ref = db(abs(hilbert(REF.rf)));
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

pars.bw          = [3 9]; % [MHz]
pars.overlap     = 0.8;
pars.blocksize   = 10; % wavelengths

pars.z_roi       = [1.0 4.5]*1E-2; % [m] 
pars.x_roi       = [-1.5 1.5]*1E-2; % [m] 

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

% Robust (optional)
par_rpl.m_est   = 1;
par_rpl.c       = [1, 100, NaN, 100]; % c_fid, c_b, c_n, c_a
par_rpl.sigma   = [0.1, 0.1, NaN, 0.1]; % sigma_fid, sigma_b, sigma_n, sigma_a

% Parameters for RPL-TV
%     mu_b  = 1E0-1E1; % ORIGINALS
%     mu_n  = 1E3-1E5; % ORIGINALS
%     mu_a  = 1E3-1E5; % ORIGINALS

% mu_rpl_tv    = [10^1.75, 1E4, 1E4]; % [mu_b, mu_n, mu_a] % Test regular    
mu_rpl_tv    = [1E3, 1E4, 1E4]; % [mu_b, mu_n, mu_a] % Test regular    

%% POWER SPECTRA ESTIMATION

% spectralData_sam_the = calc_powerSpectra(SAM, pars);
spectralData_sam_the = calc_powerSpectra_vSimple(SAM, pars); % @
S_sam = spectralData_sam_the.powerSpectra;

% spectralData_ref_the = calc_powerSpectra(REF, pars);
spectralData_ref_the = calc_powerSpectra_vSimple(REF, pars); % @
S_ref = spectralData_ref_the.powerSpectra;

%% 2DoF n 
% SR and matrix creation (1st dimension freq channels) Diagonal matrixes

SR_emz = S_sam ./ S_ref;

SR = permute(SR_emz,[3,1,2]);

[r,p,q] = size(SR);

% Frequency vector (column vector of size [r, 1])
band = spectralData_sam_the.band;
f = band(:); % Ensure f is a column vector [r,1]

X = kron( speye(p*q), ones(size(f)) ); % Size: [p*q*r, p*q] 
Z = kron( speye(p*q), log(f) );        % Size: [p*q*r, p*q] 
W = kron( speye(p*q), -4*f );          % Size: [p*q*r, p*q] 

% 2DoF prior n COMPENSATION

comp_ref = comp_ref_n_bsc(delta_n,band,p,q);

SR_comp = SR.*comp_ref;

% log-spectrum Ratio Y = X.b + Z.n + W.a
Y = log(SR_comp);

% RPL 2DoF prior n

% initialization for RPL-based methods
u_0 = initialize_rpl_n_prior(Y,X,W,mu_rpl_tv,par_rpl);

[u_opt,~] = rpl_tv_n_prior(Y,X,W,mu_rpl_tv,u_0,par_rpl);

% \Deltas  
b = u_opt(1:p*q);
a = u_opt(p*q+1:2*p*q);

% RESHAPE PARAMETERS
dy = 0.5*(diag(ones(p-1,1),1) - diag(ones(p-1,1),-1));
dy(1,1) = -1; dy(1,2) = 1; dy(end,end) = 1; dy(end,end-1) = -1;
Dy = sparse(kron(speye(q),dy));

z = 1E2*repmat(spectralData_sam_the.depth, 1, q);
dz = reshape(Dy*z(:),p,q);
dz(end,:) = dz(end-1,:);

b_ratio     = reshape(exp(b),p,q);
b_ratio_dB  = 10*log10(b_ratio);
alpha_ratio = reshape(Np2dB*Dy*a./dz(:),p,q);
n_ratio     = delta_n*ones(size(b_ratio));

% METRICS 2DoF - n
acs_sam = alpha_ratio + alpha_ref_the;

calc2dStats = {@(x) mean(x(:)), @(x) std(x(:)), @(x) 100 * std(x(:)) / mean(x(:))};

mean2d = @(x) mean(x(:));
std2d = @(x) std(x(:));
cv2d = @(x) 100*std(x(:))/mean(x(:));

m_n = mean2d(n_ratio);    s_n = std2d(n_ratio);    cv_n = cv2d(n_ratio);
m_b = mean2d(b_ratio_dB); s_b = std2d(b_ratio_dB); cv_b = cv2d(b_ratio_dB);
m_a = mean2d(acs_sam);    s_a = std2d(acs_sam);    cv_a = cv2d(acs_sam);

fprintf('----2DoF n----\n');
fprintf('α_s [dB/cm/MHz]: %.3f ± %.4f, %%CV = %.4f\n', round(m_a, 3), round(s_a, 4), round(cv_a, 4));
fprintf('b_s/b_r [dB]   : %.3f ± %.4f, %%CV = %.4f\n', round(m_b, 3), round(s_b, 4), round(cv_b, 4));
fprintf('Δn [a.u.]      : %.4f ± %.4f, %%CV = %.4f\n', round(m_n, 4), round(s_n, 4), round(cv_n, 4));
fprintf('--------------\n');


metric_a = [m_a, s_a, cv_a];
metric_b = [m_b, s_b, cv_b];
metric_n = [m_n, s_n, cv_n];

% IMAGESC PLOTS
Xaxis = spectralData_ref_the.lateral;
Zaxis = spectralData_ref_the.depth;
cm = 1e2;

axis_n = [0 1.2];
axis_a = [0 3];
axis_b = [-60 0]; % dB

fontSize = 16;
%% IMAGESC 3DoF
figure, 
set(gcf,'units','normalized','outerposition',[0 0.15 1 0.75]); box on;
sgtitle('\bf2DoF n', 'FontSize', fontSize+2)

subplot(1,3,1)
imagesc(Xaxis*cm, Zaxis*cm, acs_sam), colorbar
axis("image");
xlabel('Lateral [cm]'), ylabel('Depth [cm]'), colormap turbo;
% title(['\alpha : ', num2str(round(m_a, 3)), ' +/- ', num2str(round(s_a, 2)), ', %CV = ', num2str(round(cv_a, 3))])
title({'$\alpha_s$', ...
       [num2str(round(m_a, 3)), ' $\pm$ ', num2str(round(s_a, 3)), ', CV = ', num2str(round(cv_a, 3))]}, ...
      'Interpreter', 'latex');
h2 = colorbar; 
ylabel(h2,'[dB\cdotcm^{-1}\cdotMHz^{-1}]','FontSize', fontSize);
set(gca,'fontsize',fontSize)

subplot(1,3,2)
imagesc(Xaxis*cm, Zaxis*cm, b_ratio)
h2 = colorbar; 
if plotBSCdB 
   imagesc(Xaxis*cm, Zaxis*cm, b_ratio_dB)
   h2 = colorbar;
   ylabel(h2,'[dB]','FontSize', fontSize);
end
axis("image");
xlabel('Lateral [cm]'), colormap turbo;
% title(['\Delta b: ', num2str(round(m_b, 3)), ' +/- ', num2str(round(s_b, 2)), ', %CV = ', num2str(round(cv_b, 3))])
title({'$\Delta b = b_s / b_r$', ...
       [num2str(round(m_b, 3)), ' $\pm$ ', num2str(round(s_b, 3)), ', CV = ', num2str(round(cv_b, 3))]}, ...
      'Interpreter', 'latex');
set(gca,'fontsize',fontSize)

subplot(1,3,3)
imagesc(Xaxis*cm, Zaxis*cm, n_ratio), colorbar
h2 = colorbar;
ylabel(h2,'[a.u.]','FontSize', fontSize);
axis("image");
xlabel('Lateral [cm]'), colormap turbo;
title(['\Delta n: ', num2str(round(m_n, 3)), ' +/- ', num2str(round(s_n, 3)), ', %CV = ', num2str(round(cv_n, 3))])
title({'$\Delta n = n_s - n_r$', ...
       [num2str(round(m_n, 3)), ' $\pm$ ', num2str(round(s_n, 3)), ', CV = ', num2str(round(cv_n, 3))]}, ...
      'Interpreter', 'latex');
set(gca,'fontsize',16)

%%%%%%%%%%%%%% FORM BSC %%%%%%%%%%%%%%
freq_bsc        = spectralData_sam_the.band; % Given choosen BW

b_sam           = median(b_ratio * b_ref_the, 'all');
n_sam           = median(n_ratio + n_ref_the, 'all');
BSC_2dof        = b_sam *(freq_bsc.^n_sam);

%% 2DoF n prior delta imoverlay plots

% FIGURE INTERP OVERLAY BMODE, DELTA SNR, ACS, DELTA BSC, DELTA N
fontSize = 14;

%%%%%%%%%%%%%%%%%%%%%%%%%% Delta alpha (ACS) %%%%%%%%%%%%%%%%%%%%%%%%%%

acs_sam = alpha_ratio + alpha_ref_the;

units           = 1E2;
bmodeFull       = bmode_sam;
colorImg        = acs_sam;
range_bmode     = [-60 0];
range_img       = [0 0.7];
transparency    = 0.65;
x_img           = spectralData_sam_the.lateral*units;
z_img           = spectralData_sam_the.depth*units;
xFull           = SAM.x*units;
zFull           = SAM.z*units;
[X, Z] = meshgrid(xFull, zFull);
roi = and(X >= x_img(1), X <= x_img(end)) & ...
      and(Z >= z_img(1), Z <= z_img(end));

figure,
    [~,hB,hColor] = imOverlayInterp(bmodeFull, colorImg, range_bmode, range_img, ...
                        transparency, x_img, z_img, roi, xFull, zFull);   
    hold on;
    axis("image")
    contour(xFull, zFull, roi, 1,'w--', 'LineWidth', 2)
    hold off;
    xlabel('Lateral [cm]'), ylabel('Axial [cm]');
    hColor.Label.String = '[dB\cdotcm^{-1}\cdotMHz^{-1}]';
    title('$\alpha_s$', 'Interpreter', 'latex');
    set(gca,'fontsize',fontSize)
%%%%%%%%%%%%%%%%%%%%%%%%%% Delta alpha (ACS) %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% b ratio (dB) %%%%%%%%%%%%%%%%%%%%%%%%%%

plotBSCdB = true;
units           = 1E2;
bmodeFull       = bmode_sam;
colorImg        = b_ratio;
range_bmode     = [-60 0];
range_img       = [-21 -19];
transparency    = 0.65;
x_img           = spectralData_sam_the.lateral*units;
z_img           = spectralData_sam_the.depth*units;
xFull           = SAM.x*units;
zFull           = SAM.z*units;
[X, Z] = meshgrid(xFull, zFull);
roi = and(X >= x_img(1), X <= x_img(end)) & ...
      and(Z >= z_img(1), Z <= z_img(end));

    if plotBSCdB 
       colorImg = b_ratio_dB;
    end

figure, 
    [~,hB,hColor] = imOverlayInterp(bmodeFull, colorImg, range_bmode, range_img, ...
                        transparency, x_img, z_img, roi, xFull, zFull);
    hold on;
    contour(xFull, zFull, roi, 1,'w--', 'LineWidth', 2)
    hold off;
    xlabel('Lateral [cm]'), ylabel('Axial [cm]');
    hColor.Label.String = '';
        if plotBSCdB 
            hColor.Label.String ='[dB]';
        end
    % title('\Deltab')
    title('$\Delta b = {b_s}/{b_r}$', 'Interpreter', 'latex')
    set(gca,'fontsize',fontSize)
%%%%%%%%%%%%%%%%%%%%%%%%%% Delta b %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% Delta n %%%%%%%%%%%%%%%%%%%%%%%%%%

units           = 1E2;
bmodeFull       = bmode_sam;
colorImg        = n_ratio;
range_bmode     = [-60 0];
range_img       = [1.3 1.8];
transparency    = 0.65;
x_img           = spectralData_sam_the.lateral*units;
z_img           = spectralData_sam_the.depth*units;
xFull           = SAM.x*units;
zFull           = SAM.z*units;
[X, Z] = meshgrid(xFull, zFull);
roi = and(X >= x_img(1), X <= x_img(end)) & ...
      and(Z >= z_img(1), Z <= z_img(end));

figure,
    [~,hB,hColor] = imOverlayInterp(bmodeFull, colorImg, range_bmode, range_img, ...
                        transparency, x_img, z_img, roi, xFull, zFull);
    hold on;
    contour(xFull, zFull, roi, 1,'w--', 'LineWidth', 2)
    title('$\Delta n = n_s - n_r$', 'Interpreter', 'latex')
    hold off;
    xlabel('Lateral [cm]'), ylabel('Axial [cm]');
    hColor.Label.String = '[a.u.]';
    set(gca,'fontsize',fontSize)
    
%%%%%%%%%%%%%%%%%%%%%%%%%% Delta n %%%%%%%%%%%%%%%%%%%%%%%%%%


%% 3DoF 
% SR and matrix creation (1st dimension freq channels) Diagonal matrixes

SR_emz = S_sam ./ S_ref;

SR = permute(SR_emz,[3,1,2]);

[r,p,q] = size(SR);

% Frequency vector (column vector of size [r, 1])
band = spectralData_sam_the.band;
f = band(:); % Ensure f is a column vector [r,1]

X = kron( speye(p*q), ones(size(f)) ); % Size: [p*q*r, p*q] 
Z = kron( speye(p*q), log(f) );        % Size: [p*q*r, p*q] 
W = kron( speye(p*q), -4*f );          % Size: [p*q*r, p*q] 

%
% log-spectrum Ratio Y = X.b + Z.n + W.a
Y = log(SR);

% RPL
% initialization for RPL-based methods
u_0 = initialize_rpl(Y,X,Z,W,mu_rpl_tv,par_rpl);

[u_opt,~] = rpl_tv(Y,X,Z,W,mu_rpl_tv,u_0,par_rpl);

b = u_opt(1:p*q);
n = u_opt(p*q+1:2*p*q);
a = u_opt(2*p*q+1:end);

% RESHAPE PARAMETERS
dy = 0.5*(diag(ones(p-1,1),1) - diag(ones(p-1,1),-1));
dy(1,1) = -1; dy(1,2) = 1; dy(end,end) = 1; dy(end,end-1) = -1;
Dy = sparse(kron(speye(q),dy));

z = 1E2*repmat(spectralData_sam_the.depth, 1, q);
dz = reshape(Dy*z(:),p,q);
dz(end,:) = dz(end-1,:);

b_ratio     = reshape(exp(b),p,q);
b_ratio_dB  = 10*log10(b_ratio);
alpha_ratio = reshape(Np2dB*Dy*a./dz(:),p,q);
n_ratio     = reshape(n,p,q); % 3dof

% METRICS 3 DoF
acs_sam = alpha_ratio + alpha_ref_the;

mean2d = @(x) mean(x(:));
std2d = @(x) std(x(:));
cv2d = @(x) 100*std(x(:))/mean(x(:));

m_n = mean2d(n_ratio);    s_n = std2d(n_ratio);    cv_n = cv2d(n_ratio);
m_b = mean2d(b_ratio_dB); s_b = std2d(b_ratio_dB); cv_b = cv2d(b_ratio_dB);
m_a = mean2d(acs_sam);    s_a = std2d(acs_sam);    cv_a = cv2d(acs_sam);

fprintf('-----3DoF-----\n');
fprintf('α_s [dB/cm/MHz]: %.3f ± %.4f, %%CV = %.4f\n', round(m_a, 3), round(s_a, 4), round(cv_a, 4));
fprintf('b_s/b_r [dB]   : %.3f ± %.4f, %%CV = %.4f\n', round(m_b, 3), round(s_b, 4), round(cv_b, 4));
fprintf('Δn [a.u.]      : %.4f ± %.4f, %%CV = %.4f\n', round(m_n, 4), round(s_n, 4), round(cv_n, 4));
fprintf('--------------\n');

metric_a = [m_a, s_a, cv_a];
metric_b = [m_b, s_b, cv_b];
metric_n = [m_n, s_n, cv_n];

% IMAGESC PLOTS
Xaxis = spectralData_ref_the.lateral;
Zaxis = spectralData_ref_the.depth;
cm = 1e2;

axis_n = [0 1.2];
axis_a = [0 3];
axis_b = [-60 0]; % dB
fontSize = 16;

%% IMAGESC 2DoF n
figure, 
set(gcf,'units','normalized','outerposition',[0 0.15 1 0.75]); box on;
sgtitle('\bf3DoF', 'FontSize', fontSize+2)

subplot(1,3,1)
imagesc(Xaxis*cm, Zaxis*cm, acs_sam), colorbar
axis("image");
xlabel('Lateral [cm]'), ylabel('Depth [cm]'), colormap turbo;
% title(['\alpha : ', num2str(round(m_a, 3)), ' +/- ', num2str(round(s_a, 2)), ', %CV = ', num2str(round(cv_a, 3))])
title({'$\alpha_s$', ...
       [num2str(round(m_a, 3)), ' $\pm$ ', num2str(round(s_a, 3)), ', CV = ', num2str(round(cv_a, 3))]}, ...
      'Interpreter', 'latex');
h2 = colorbar; 
ylabel(h2,'[dB\cdotcm^{-1}\cdotMHz^{-1}]','FontSize', fontSize);
set(gca,'fontsize',fontSize)

subplot(1,3,2)
imagesc(Xaxis*cm, Zaxis*cm, b_ratio)
h2 = colorbar; 
if plotBSCdB 
   imagesc(Xaxis*cm, Zaxis*cm, b_ratio_dB)
   h2 = colorbar;
   ylabel(h2,'[dB]','FontSize', fontSize);
end
axis("image");
xlabel('Lateral [cm]'), colormap turbo;
% title(['\Delta b: ', num2str(round(m_b, 3)), ' +/- ', num2str(round(s_b, 2)), ', %CV = ', num2str(round(cv_b, 3))])
title({'$\Delta b = b_s / b_r$', ...
       [num2str(round(m_b, 3)), ' $\pm$ ', num2str(round(s_b, 3)), ', CV = ', num2str(round(cv_b, 3))]}, ...
      'Interpreter', 'latex');
set(gca,'fontsize',fontSize)

subplot(1,3,3)
imagesc(Xaxis*cm, Zaxis*cm, n_ratio), colorbar
h2 = colorbar;
ylabel(h2,'[a.u.]','FontSize', fontSize);
axis("image");
xlabel('Lateral [cm]'), colormap turbo;
title(['\Delta n: ', num2str(round(m_n, 3)), ' +/- ', num2str(round(s_n, 3)), ', %CV = ', num2str(round(cv_n, 3))])
title({'$\Delta n = n_s - n_r$', ...
       [num2str(round(m_n, 3)), ' $\pm$ ', num2str(round(s_n, 3)), ', CV = ', num2str(round(cv_n, 3))]}, ...
      'Interpreter', 'latex');
set(gca,'fontsize',16)

%%%%%%%%%%%%%% FORM BSC %%%%%%%%%%%%%%
freq_bsc        = spectralData_sam_the.band; % Given choosen BW

b_sam           = median(b_ratio * b_ref_the, 'all');
n_sam           = median(n_ratio + n_ref_the, 'all');
BSC_3dof        = b_sam *(freq_bsc.^n_sam);

%% 3DoF imoverlay plots

% FIGURE INTERP OVERLAY BMODE, DELTA SNR, ACS, DELTA BSC, DELTA N
fontSize = 14;

%%%%%%%%%%%%%%%%%%%%%%%%%% Delta alpha (ACS) %%%%%%%%%%%%%%%%%%%%%%%%%%

acs_sam = alpha_ratio + alpha_ref_the;

units           = 1E2;
bmodeFull       = bmode_sam;
colorImg        = acs_sam;
range_bmode     = [-60 0];
range_img       = [0 0.7];
transparency    = 0.65;
x_img           = spectralData_sam_the.lateral*units;
z_img           = spectralData_sam_the.depth*units;
xFull           = SAM.x*units;
zFull           = SAM.z*units;
[X, Z] = meshgrid(xFull, zFull);
roi = and(X >= x_img(1), X <= x_img(end)) & ...
      and(Z >= z_img(1), Z <= z_img(end));

figure,
    [~,hB,hColor] = imOverlayInterp(bmodeFull, colorImg, range_bmode, range_img, ...
                        transparency, x_img, z_img, roi, xFull, zFull);   
    hold on;
    axis("image")
    contour(xFull, zFull, roi, 1,'w--', 'LineWidth', 2)
    hold off;
    xlabel('Lateral [cm]'), ylabel('Axial [cm]');
    hColor.Label.String = '[dB\cdotcm^{-1}\cdotMHz^{-1}]';
    title('$\alpha_s$', 'Interpreter', 'latex');
    set(gca,'fontsize',fontSize)
%%%%%%%%%%%%%%%%%%%%%%%%%% Delta alpha (ACS) %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% b ratio (dB) %%%%%%%%%%%%%%%%%%%%%%%%%%

plotBSCdB = true;
units           = 1E2;
bmodeFull       = bmode_sam;
colorImg        = b_ratio;
range_bmode     = [-60 0];
range_img       = [-21 -19];
transparency    = 0.65;
x_img           = spectralData_sam_the.lateral*units;
z_img           = spectralData_sam_the.depth*units;
xFull           = SAM.x*units;
zFull           = SAM.z*units;
[X, Z] = meshgrid(xFull, zFull);
roi = and(X >= x_img(1), X <= x_img(end)) & ...
      and(Z >= z_img(1), Z <= z_img(end));

    if plotBSCdB 
       colorImg = b_ratio_dB;
    end

figure, 
    [~,hB,hColor] = imOverlayInterp(bmodeFull, colorImg, range_bmode, range_img, ...
                        transparency, x_img, z_img, roi, xFull, zFull);
    hold on;
    contour(xFull, zFull, roi, 1,'w--', 'LineWidth', 2)
    hold off;
    xlabel('Lateral [cm]'), ylabel('Axial [cm]');
    hColor.Label.String = '';
        if plotBSCdB 
            hColor.Label.String ='[dB]';
        end
    % title('\Deltab')
    title('$\Delta b = {b_s}/{b_r}$', 'Interpreter', 'latex')
    set(gca,'fontsize',fontSize)
%%%%%%%%%%%%%%%%%%%%%%%%%% Delta b %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% Delta n %%%%%%%%%%%%%%%%%%%%%%%%%%

units           = 1E2;
bmodeFull       = bmode_sam;
colorImg        = n_ratio;
range_bmode     = [-60 0];
range_img       = [1.3 1.8];
transparency    = 0.65;
x_img           = spectralData_sam_the.lateral*units;
z_img           = spectralData_sam_the.depth*units;
xFull           = SAM.x*units;
zFull           = SAM.z*units;
[X, Z] = meshgrid(xFull, zFull);
roi = and(X >= x_img(1), X <= x_img(end)) & ...
      and(Z >= z_img(1), Z <= z_img(end));

figure,
    [~,hB,hColor] = imOverlayInterp(bmodeFull, colorImg, range_bmode, range_img, ...
                        transparency, x_img, z_img, roi, xFull, zFull);
    hold on;
    contour(xFull, zFull, roi, 1,'w--', 'LineWidth', 2)
    title('$\Delta n = n_s - n_r$', 'Interpreter', 'latex')
    hold off;
    xlabel('Lateral [cm]'), ylabel('Axial [cm]');
    hColor.Label.String = '[a.u.]';
    set(gca,'fontsize',fontSize)
    
%%%%%%%%%%%%%%%%%%%%%%%%%% Delta n %%%%%%%%%%%%%%%%%%%%%%%%%%
%% Litle HELP THERETICAL

BSC_gt = b_sam_the *(freq_bsc.^n_sam_the);

% PLOT 3DoF and 2DoF n same time

% Metrics
BSC_gt_dB   = 10*log10(BSC_gt);
BSC_2dof_dB = 10*log10(BSC_2dof);
BSC_3dof_dB = 10*log10(BSC_3dof);

%%%%%%%%%% BSC 2DoF n METRICS %%%%%%%%%%
metrics2dof.rmse    = rmse(BSC_2dof, BSC_gt);
metrics2dof.rmse_dB = rmse(BSC_2dof, BSC_gt);
metrics2dof.diff_dB = mean( abs(BSC_2dof_dB - BSC_gt_dB) ); 

err     = (BSC_2dof - BSC_gt)/BSC_gt;
abs_err = abs(BSC_2dof - BSC_gt)/BSC_gt;
err_sq  = (BSC_2dof - BSC_gt).^2;

metrics2dof.mpe     = mean(err);
metrics2dof.mae     = mean(abs_err);
metrics2dof.rmse2   = sqrt(mean(err_sq));
metrics2dof.nrmse2  = metrics2dof.rmse / mean(BSC_2dof); 

%%%%%%%%%% BSC 3DoF METRICS %%%%%%%%%%
metrics3dof.rmse    = rmse(BSC_3dof, BSC_gt);
metrics3dof.rmse_dB = rmse(BSC_3dof, BSC_gt);
metrics3dof.diff_dB = mean( abs(BSC_3dof_dB - BSC_gt_dB) ); 

err     = (BSC_3dof - BSC_gt)/BSC_gt;
abs_err = abs(BSC_3dof - BSC_gt)/BSC_gt;
err_sq  = (BSC_3dof - BSC_gt).^2;

metrics3dof.mpe     = mean(err);
metrics3dof.mae     = mean(abs_err);
metrics3dof.rmse2   = sqrt(mean(err_sq));
metrics3dof.nrmse2  = metrics3dof.rmse2 / mean(BSC_2dof); 


lw = 2;
figure,
% set(gcf,'units','normalized','outerposition',[0 0 1 1]);
h1=loglog(band, BSC_3dof, '-', 'LineWidth', lw+1, 'DisplayName', '3-DoF'); hold on
h2=loglog(band, BSC_2dof, '-', 'LineWidth', lw+1, 'DisplayName', '2-DoF');
h3=loglog(band, BSC_gt, 'k-.', 'LineWidth', lw, 'DisplayName','Ground truth');
grid on;

% Add a textbox to display RMSE (3DoF)
gof_text_3dof = ['\bfRMSE (3DoF) = ', num2str(metrics3dof.rmse, '%.2f')];
annotation('textbox', [0.15, 0.8, 0.3, 0.1], 'String', gof_text_3dof, ...
    'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', ...
    'FontSize', 10, 'Color', h1.Color); % Match font color to line color

% Add a textbox to display RMSE (2DoF)
gof_text_2dof = ['\bfRMSE (2DoF) = ', num2str(metrics2dof.rmse, '%.2f')];
annotation('textbox', [0.15, 0.7, 0.3, 0.1], 'String', gof_text_2dof, ...
    'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', ...
    'FontSize', 10, 'Color', h2.Color); % Match font color to line color

xlabel('Frequency [MHz]'), 
ylabel('BSC [cm^{-1}\cdot sr^{-1}]')
title(['BSC: ', name_sam])
legend ('Location', 'Best')
set(gca,'fontsize',fontSize)

%%
% dir_figures = './Hito2_PROCIENCA2024/figures';
% dir_figures = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\Hito2_PROCIENCA2024\figures';
% title_figures = ['Simu_', ( name_sam ), '_fig'];
% if ~exist(dir_figures,"dir"); mkdir(dir_figures); end
% save_all_figures_to_directory(dir_figures, title_figures);
