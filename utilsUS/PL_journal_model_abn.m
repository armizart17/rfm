%%% GENERAL TEST FOR JOURNAL PL 3 DOF vs 2DOF with simulations%%%
% 26/09/2023
%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%
% clear all; 
% close all

% close all;
addpath('.\RPL_TV reduced\')
pathData = '.\NewPL_data\data_ok\';

a_list = [0.4, 0.6];
b_list = [0.01 1];
n_list = [0 1.5];

% numCase = input('Enter type of case: ');
numCase = 3;
% estim_method = input('Enter type of method (1) for 3DOF, (2) *a* prior 2DOF, (3) *b* prior 2DOF, (4) *n* prior: ');
% estim_method = 2;
% blocksize = input('Enter the blocksize in wavelenghts: ');
% regu = input('Regularization (1) True (0) False [0.01 0.01 0.01] ');

l_estim_method = [1 2 3 4];
for lll = 1:length(l_estim_method)
    estim_method = l_estim_method(lll);
% estim_method = 2;
blocksize = 10; 
regu = 1;

switch numCase
    case 1 % ONLY VARIATION a
        a_sam = 0.6; a_ref = 0.4;
        b_sam = 0.01; b_ref = 0.01;
        n_sam = 1.5; n_ref = 1.5;
    case 2  % ONLY VARIATION a, b
        a_sam = 0.6; a_ref = 0.4;
        b_sam = 0.01; b_ref = 1;
        n_sam = 1.5; n_ref = 1.5;
    case 3  % ONLY VARIATION a, b, n
        a_sam = 0.6; a_ref = 0.4;
        b_sam = 0.01; b_ref = 1;
        n_sam = 1.5; n_ref = 0;
end

switch regu

    case 0
        struct.mu = [0.001 0.001 0.001];
    case 1
        struct.mu = [1e3, 1e3, 1e5] ; % mu_b, mu_n mu_a
        % struct.mu = [1e3, 1e5, 1e5] ; % mu_b, mu_n mu_a

        % for TUFFC
        % struct.mu = [1e3, 1e3, 1e4];
        % struct.mu = [1e3, 1e3, 1e4];
%     mu_b  = 1e+0;%1e+0; %1e+0; % ORIGINALS
%     mu_n  = 1e+3;%1e+3; %1e+2; % ORIGINALS
%     mu_a  = 1e+3;%1e+3; %1e+3; % ORIGINALS
end

%

AUX = load([pathData,'ref0.mat']); % JUST AUXILIAR TO GIVE AXIS
cm = 1e2;

name_sam = ['a',num2str(a_sam),'b',num2str(b_sam),'n',num2str(n_sam)];
name_ref = ['a',num2str(a_ref),'b',num2str(b_ref),'n',num2str(n_ref)];

SAM = load([pathData, 'ref_', name_sam]);
SAM.rf = SAM.rfnew;
SAM.x = AUX.x;
SAM.z = AUX.z;
SAM.fs = AUX.fs;
clear AUX

REF = load([pathData, 'ref_', name_ref]);

bmode_caxix = [-60 0];

% figure, 
% set(gcf,'units','normalized','outerposition',[0 0.25 0.5 0.5]); box on;
% subplot (1,2,1)
% imagesc( SAM.x*cm, SAM.z*cm, my_RF2Bmode(SAM.rfnew), bmode_caxix), colormap gray
% title(['SAM: ', name_sam]), xlabel('Lateral [cm]'), ylabel('Depth [cm]');
% set(gca,'fontsize',18)
% 
% subplot (1,2,2)
% imagesc( SAM.x*cm, SAM.z*cm, my_RF2Bmode(REF.rfnew), bmode_caxix), colormap gray
% title(['REF: ', name_ref]), xlabel('Lateral [cm]'), ylabel('Depth [cm]');
% set(gca,'fontsize',18)

% ADJUST DATA FOR HECTOR CODE

Fs = SAM.fs;
BmodeWidth = 0.3*128; % [mm]
phan_alpha = 0; % unknown
bsc_band = 0; % unknown
cm = 1e2;
y1 = 10;
y2 = 45;
x1 = -15;
x2 = 15;
c0 = 1540;
BW = [3 9];
tissue = SAM.rf;
SAM.RF = SAM.rf;
phantom = REF.rfnew;
line_ROI_lateral = blocksize;
lambda_axial = blocksize;
delta_flag = 1;

delta_n = nan;
delta_b = nan;
delta_alpha = nan;

switch estim_method
    case 1 % 3DOF
        
    case 2 % 2DOF a prior
        delta_alpha = a_sam - a_ref;    
    case 3 % 2DOF b prior
        delta_b = b_sam / b_ref;
    case 4 % 2DOF n prior
        delta_n = n_sam - n_ref;
end

 [Est_n_sample, Est_b_sample, Est_alphaz_sample,delta_SNR,SNR,est_CF,DATA,Xaxis,Zaxis,SR,z,...
    nb_ROI_axial,nb_ROI_lateral,mu_rpl_tv,mu_rpl_robust_tv] = qUS_estimation_v2(Fs, BmodeWidth, phan_alpha, bsc_band,...
    x1, x2, y1,y2, c0, BW, tissue, phantom, line_ROI_lateral, lambda_axial, estim_method, delta_flag, SAM, delta_n, struct, ...
    delta_b, delta_alpha, a_ref);

%  [Est_n_sample, Est_b_sample, Est_alphaz_sample,delta_SNR,SNR,est_CF,DATA,Xaxis,Zaxis,SR,z,...
%     nb_ROI_axial,nb_ROI_lateral,mu_rpl_tv,mu_rpl_robust_tv] = qUS_estimation(Fs, BmodeWidth, phan_alpha, bsc_band,...
%     x1, x2, y1,y2, c0, BW, tissue, phantom, line_ROI_lateral, lambda_axial, estim_method, delta_flag, SAM, delta_n, struct);


% Est_b_sample = 10*log10(Est_b_sample); % if want in dB


mean2d = @(x) mean(x(:));
std2d = @(x) std(x(:));
cv2d = @(x) 100*std(x(:))/mean(x(:));

%% DELTA METRICS
% m_n = mean2d(Est_n_sample);
% s_n = std2d(Est_n_sample);
% cv_n = cv2d(Est_n_sample);
% 
% m_b = mean2d(Est_b_sample);
% s_b = std2d(Est_b_sample);
% cv_b = cv2d(Est_b_sample);
% 
% m_a = mean2d(Est_alphaz_sample);
% s_a = std2d(Est_alphaz_sample);
% cv_a = cv2d(Est_alphaz_sample);
% 
% disp(['Case: ' , num2str(numCase),' | Wavelength: ', num2str(blocksize), ...
%      ' | Method: ',  num2str(estim_method), ' | Regu: ', num2str(regu)]);
% disp(['\Delta alpha : ', num2str(round(m_a, 3)), ' +/- ', ...
%     num2str(round(s_a, 4)), ', %CV = ', num2str(round(cv_a, 4))]); 
% 
% disp(['\Delta b : ', num2str(round(m_b, 3)), ' +/- ', ...
%     num2str(round(s_b, 4)), ', %CV = ', num2str(round(cv_b, 4))]); 
% 
% disp(['\Delta n : ', num2str(round(m_n, 4)), ' +/- ', ...
%     num2str(round(s_n, 4)), ', %CV = ', num2str(round(cv_n, 4))]);
% 
% disp('--------')
% 
% metric_a = [m_a, s_a, cv_a];
% metric_b = [m_b, s_b, cv_b];
% metric_n = [m_n, s_n, cv_n];
% 
% %
% axis_n = [0 1.7];
% axis_b = [0 0.11];
% 
% axis_a = [0 0.25];
% 
% axis_b = [-60 0]; % dB
% 
% 
% figure, 
% set(gcf,'units','normalized','outerposition',[0 0.25 1 0.65]); box on;
% 
% subplot(1,3,1)
% % imagesc(Xaxis*cm, Zaxis*cm, Est_alphaz_sample, axis_a), colorbar
% imagesc(Xaxis*cm, Zaxis*cm, Est_alphaz_sample), colorbar
% xlabel('Lateral [cm]'), ylabel('Depth [cm]'), colormap jet;
% title(['\Delta \alpha : ', num2str(round(m_a, 3)), ' +/- ', num2str(round(s_a, 2)), ', %CV = ', num2str(round(cv_a, 3))])
% set(gca,'fontsize',16)
% 
% subplot(1,3,2)
% % imagesc(Xaxis*cm, Zaxis*cm, Est_b_sample, axis_b), colorbar
% imagesc(Xaxis*cm, Zaxis*cm, Est_b_sample), colorbar
% xlabel('Lateral [cm]'), ylabel('Depth [cm]'), colormap jet;
% title(['\Delta b : ', num2str(round(m_b, 3)), ' +/- ', num2str(round(s_b, 2)), ', %CV = ', num2str(round(cv_b, 3))])
% set(gca,'fontsize',16)
% 
% subplot(1,3,3)
% % imagesc(Xaxis*cm, Zaxis*cm, Est_n_sample, axis_n), colorbar
% imagesc(Xaxis*cm, Zaxis*cm, Est_n_sample), colorbar
% xlabel('Lateral [cm]'), ylabel('Depth [cm]'), colormap jet;
% title(['\Delta n : ', num2str(round(m_n, 3)), ' +/- ', num2str(round(s_n, 3)), ', %CV = ', num2str(round(cv_n, 3))])
% set(gca,'fontsize',16)





%% NOT ALPHA, REF KNOWN (Feb 2025)
Est_alphaz_sample_sam = Est_alphaz_sample + a_ref;
Est_b_sample = 10*log10(Est_b_sample);

m_n = mean2d(Est_n_sample);
s_n = std2d(Est_n_sample);
cv_n = cv2d(Est_n_sample);

m_b = mean2d(Est_b_sample);
s_b = std2d(Est_b_sample);
cv_b = cv2d(Est_b_sample);

m_a = mean2d(Est_alphaz_sample_sam);
s_a = std2d(Est_alphaz_sample_sam);
cv_a = cv2d(Est_alphaz_sample_sam);

disp(['Case: ' , num2str(numCase),' | Wavelength: ', num2str(blocksize), ...
     ' | Method: ',  num2str(estim_method), ' | Regu: ', num2str(regu)]);
disp(['\Delta alpha : ', num2str(round(m_a, 3)), ' +/- ', ...
    num2str(round(s_a, 4)), ', %CV = ', num2str(round(cv_a, 4))]); 

disp(['\Delta b : ', num2str(round(m_b, 3)), ' +/- ', ...
    num2str(round(s_b, 4)), ', %CV = ', num2str(round(cv_b, 4))]); 

disp(['\Delta n : ', num2str(round(m_n, 4)), ' +/- ', ...
    num2str(round(s_n, 4)), ', %CV = ', num2str(round(cv_n, 4))]);

disp('--------')

metric_a = [m_a, s_a, cv_a];
metric_b = [m_b, s_b, cv_b];
metric_n = [m_n, s_n, cv_n];

%
axis_n = [0 1.7];
axis_b = [0 0.11];

axis_a = [0 0.25];

axis_b = [-60 0]; % dB

%
figure, 
set(gcf,'units','normalized','outerposition',[0 0.25 1 0.65]); box on;

subplot(1,3,1)
% imagesc(Xaxis*cm, Zaxis*cm, Est_alphaz_sample, axis_a), colorbar
imagesc(Xaxis*cm, Zaxis*cm, Est_alphaz_sample_sam), colorbar
xlabel('Lateral [cm]'), ylabel('Depth [cm]'), colormap jet;
title(['\alpha_s : ', num2str(round(m_a, 3)), ' +/- ', num2str(round(s_a, 2)), ', %CV = ', num2str(round(cv_a, 3))])
set(gca,'fontsize',16)

subplot(1,3,2)
% imagesc(Xaxis*cm, Zaxis*cm, Est_b_sample, axis_b), colorbar
imagesc(Xaxis*cm, Zaxis*cm, Est_b_sample), colorbar
xlabel('Lateral [cm]'), ylabel('Depth [cm]'), colormap jet;
title(['\Delta b : ', num2str(round(m_b, 3)), ' +/- ', num2str(round(s_b, 2)), ', %CV = ', num2str(round(cv_b, 3))])
set(gca,'fontsize',16)

subplot(1,3,3)
% imagesc(Xaxis*cm, Zaxis*cm, Est_n_sample, axis_n), colorbar
imagesc(Xaxis*cm, Zaxis*cm, Est_n_sample), colorbar
xlabel('Lateral [cm]'), ylabel('Depth [cm]'), colormap jet;
title(['\Delta n : ', num2str(round(m_n, 3)), ' +/- ', num2str(round(s_n, 3)), ', %CV = ', num2str(round(cv_n, 3))])
set(gca,'fontsize',16)


end