% ======================================================================
% ======================================================================
%% ALL PHANTOMS
clear, clc;
close all

Np2dB = 20*log10(exp(1));
dB2Np = 1/Np2dB;
range_bmode = [-60 0];
plotBSCdB = true;

addpath(genpath(pwd))

%%

addpath(genpath(pwd));

pathData = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\phantomTargets\SonixTouch\ID316';
pathRef = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\phantomTargets\SonixTouch\refPh544\';
refPhan_alpha = 0.53; % [dB/cm/MHz]

% resultsDir = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\phantomTargets\SonixTouch\ID316_out';

targetFiles = dir([pathData,'\*.mat']);

% targetFiles = targetFiles(end-2:end); 
num_targets = length(targetFiles);

groundTruthTargets    = [0.52, 0.55, 0.74, 0.81, 0.75, 0.97, 0.95, 0.95]; % [dB/cm/MHz]
groundTruthBackground = 0.55; % [dB/cm/MHz]

% if ~exist("resultsDir","dir"); mkdir(resultsDir); end
tableName = 'targetsPL.xlsx';

%% Constants

pars.bw          = [2.5 7.5]; % [MHz]
pars.overlap     = 0.8;
pars.blocksize   = 8; % wavelengths
pars.z_roi       = [2 35]*1E-3; % [m] 
pars.x_roi       = [1 38]*1E-3; % [m] 
pars.saran_layer = false;

pars.window_type = 3; %  (1) Hanning, (2) Tuckey, (3) Hamming, (4) Tchebychev

%% REGULARIZATION PARAMETERS

% Implementation parameters
par_rpl.ini_method = 1; % METHOD LEAST SQUARES INITIALIZATION 
par_rpl.ini_tol = 1e-16;
par_rpl.df_op = 1;

par_rpl.tol   = 1e-16;
par_rpl.kmax  = 100;
par_rpl.eps_f = 1e-16;
par_rpl.m_est = 0; %Robust

% Parameters for RPL-TV
%     mu_b  = 1E0-1E1; % ORIGINALS
%     mu_n  = 1E3-1E5; % ORIGINALS
%     mu_a  = 1E3-1E5; % ORIGINALS

mu_rpl_tv    = [10^1.5, NaN, 1E4]; % [mu_b, mu_n, mu_a]
mu_rpl_tv    = [10^1.75, NaN, 1E4]; % [mu_b, mu_n, mu_a]
% mu_rpl_tv    = [1E3; 1E3; 1E5]; % [mu_b, mu_n, mu_a] % before


%% Weight parameters
ratioCutOff = 10;
reject = 0.1;
extension = 3;

% SWTV
aSNR = 5; bSNR = 0.09;
desvMin = 15;


% Plotting constants
range_acs = [0.2,1.2];

tol = 1e-3;

c1x = 1.95; c1z = 1.93;
roiL = 1; roiD = 0.6;
roiLz = 1.5;
%% For looping each phantom

for iAcq = 1:num_targets

%% LOAD SAMPLE
fprintf("-----------------------------------")
fprintf("Phantom no. %i, %s\n",iAcq,targetFiles(iAcq).name);
SAM = load(fullfile(pathData,targetFiles(iAcq).name));

%% SWITCH EMZ

switch iAcq
    case 6
        c1x = 1.8; c1z = 1.9;
    case 7
        c1x = 1.95; c1z = 1.95;
    case 8
        c1x = 1.85; c1z = 1.9;
    otherwise
        c1x = 1.85; c1z = 1.9;
end


%% POWER SPECTRA ESTIMATION

% spectralData_sam = calc_powerSpectra(SAM, pars);
spectralData_sam = calc_powerSpectra_vSimple(SAM, pars); % @

S_sam = spectralData_sam.powerSpectra;

% num_ref = 1;
% spectralData_ref = calc_powerSpectra(REF, pars);
% S_ref = spectralData_ref.powerSpectra;

refFiles = dir([pathRef,'\*.mat']);
num_ref = length(refFiles);

if num_ref > 1
    S_ref_aux = zeros([size(S_sam), num_ref]);

    for iRef = 1:num_ref
          REF = load(fullfile(pathRef, refFiles(iRef).name));
          if ~isfield(REF, 'rf')
              REF.rf = REF.RF(:,:,1);      
          end  
        % spectralData_ref = calc_powerSpectra(REF, pars);
        spectralData_ref = calc_powerSpectra_vSimple(REF, pars); % @
        S_ref_aux(:,:,:,iRef) = spectralData_ref.powerSpectra;    
    end
    S_ref = mean(S_ref_aux, 4);
end

%% SR CALCULATION AND PERMUTATION
SR_emz = S_sam ./ S_ref;

SR = permute(SR_emz,[3,1,2]);

%% 2DOF prior n
delta_n = 0;
band = spectralData_sam.band;
[r,p,q] = size(SR);
    
%     delta_n = -1.5; %2.8; %-6; % 2.5 4 nsam-nref
% bsc_band -> n = delta_n del tejido
% bsc_band -> n = n + delta_n del tejido

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

%% SETTING FOR PLOT

rInc = 0.95; % [cm]

x = SAM.x*1E2; % [cm]
z = SAM.z*1E2; % [cm]
x_ACS = spectralData_sam.lateral*1E2; % [cm]
z_ACS = spectralData_sam.depth*1E2; % [cm]

acs_sam = alpha_ratio + refPhan_alpha;

if isfield(SAM, 'rf')
    bmode_sam = db(hilbert(SAM.rf(:,:,1)));
    bmode_sam = bmode_sam - max(bmode_sam(:));    
end

if isfield(SAM, 'RF')
    bmode_sam = db(hilbert(SAM.RF(:,:,1)));
    bmode_sam = bmode_sam - max(bmode_sam(:));    
end

%%
figure('Units','centimeters', 'Position',[5 2 6 16]);
tl = tiledlayout(4,1, "Padding","tight");

t1 = nexttile;
imagesc(x, z, bmode_sam, range_bmode)
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
% xlabel('Lateral [cm]'),
ylabel('Axial [cm]')
axis image
colormap(t1,gray)
title('B-mode')
%subtitle(' ')
c = colorbar;
c.Label.String = 'dB';
% hold on 
% rectangle('Position',[c1x-rInc c1z-rInc 2*rInc 2*rInc], 'LineStyle','--', ...
%     'LineWidth',1, 'Curvature',1)
% hold off

% fontsize(gcf,8,'points')

t2 = nexttile;
imagesc(x_ACS, z_ACS, acs_sam, range_acs)
% xlabel('Lateral [cm]'),
ylabel('Axial [cm]')
colormap(t2,turbo)
axis image
title('ACS')
c = colorbar;
c.Label.String = 'dB/cm/MHz';
% fontsize(gcf,8,'points')
hold on 
rectangle('Position',[c1x-rInc c1z-rInc 2*rInc 2*rInc], 'LineStyle','--', ...
    'LineWidth',1, 'Curvature',1, 'EdgeColor','w')
hold off

t3 = nexttile;
if plotBSCdB
    imagesc(x_ACS, z_ACS, b_ratio_dB);
    c = colorbar;
    c.Label.String = 'dB';
else
    imagesc(x_ACS, z_ACS, b_ratio);
    c = colorbar;
    c.Label.String = '';
end
% xlabel('Lateral [cm]'),
ylabel('Axial [cm]')
colormap(t3,turbo)
axis image
title('\Delta BSC')

% fontsize(gcf,8,'points')
hold on 
rectangle('Position',[c1x-rInc c1z-rInc 2*rInc 2*rInc], 'LineStyle','--', ...
    'LineWidth',1, 'Curvature',1, 'EdgeColor','w')
hold off

t5 = nexttile;
imagesc(x_ACS, z_ACS, n_ratio)
xlabel('Lateral [cm]'),
ylabel('Axial [cm]')
colormap(t5,turbo)
axis image
title('\Delta n')
c = colorbar;
c.Label.String = '';

% fontsize(gcf,8,'points')
hold on 
rectangle('Position',[c1x-rInc c1z-rInc 2*rInc 2*rInc], 'LineStyle','--', ...
    'LineWidth',1, 'Curvature',1, 'EdgeColor','w')
hold off
end

%%

% pathOut = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\phantomTargets\SonixTouch\ID316_outPL';
% if ~exist("pathOut","dir"); mkdir(pathOut); end
% name = 'PL_target';
% if plotBSCdB
%     name = 'PL_b_dB_target';
% end
% save_all_figures_to_directory(pathOut, name, 'png'); % @
%%
% close all
% %%
% results1 = struct2table(MetricsTV);
% results2 = struct2table(MetricsSWTV);
% results4 = struct2table(MetricsSWIFT);
% 
% %%  % @
% T = [results1;results2;results4];
% writetable(T,fullfile(resultsDir,tableName),...
%      'WriteRowNames',true);
%%