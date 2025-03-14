%%% GENERAL TEST FOR JOURNAL PL 3 DOF vs 2DOF with simulations%%%
% 26/09/2023
%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clear all;
addpath(genpath('.\QUS_codes\'))
pathData = '.\QUS_codes\NewPL_data\';
pathData = '.\QUS_codes\NewPL_data\data_ok';

a_values = {'0.4', '0.01', ''}

a_list = [0.4, 0.6];
b_list = [0.01 1];
n_list = [0 1.5];

%% MODEL a y b
% (1) case (varia a)
% SAM
% a0.6 b0.01 n0
% REF
% a0.4 b0.01 n0

% (2) case (varia a , b)
% SAM
% a0.6 b0.01 n0
% REF
% a0.4 b1 n0

% (3) case (varia a , b , n)
% SAM
% a0.6 b0.01 n0
% REF
% a0.4 b1 n1.5

%% MODEL a b y n
% () Primer caso (varia a)
% SAM
% a0.6 b0.01 n1.5
% REF
% a0.4 b0.01 n1.5

% (2) Primer caso (varia a , b)
% SAM
% a0.6 b0.01 n1.5
% REF
% a0.4 b1 n1.5

% (2) Primer caso (varia a , b , n)
% SAM
% a0.6 b0.01 n1.5
% REF
% a0.4 b1 n0
%%

AUX = load([pathData,'ref0.mat']);

% REF = load(['.\QUS_codes\NewPL_data\ref_a0.6b0.01n1.5.mat']);
SAM = load([pathData, 'ref_a0.4b1n0.mat']);
SAM.rf = SAM.rfnew;
SAM.x = AUX.x;
SAM.z = AUX.z;
SAM.fs = AUX.fs;


REF = load([pathData, 'ref_a0.6b0.01n1.5.mat']);

bmode_caxix = [-60 0];
%%
figure, 
subplot (1,2,1)
% imagesc(SAM.x,SAM.z, my_RF2Bmode(SAM.rf), bmode_caxix), colormap gray

imagesc( my_RF2Bmode(SAM.rfnew), bmode_caxix), colormap gray
title(['Bmode SAM: ', ])

subplot (1,2,2)
imagesc( my_RF2Bmode(REF.rfnew), bmode_caxix), colormap gray
title('Bmode REF')
%% ADJUST DATA FOR HECTOR CODE



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
line_ROI_lateral = 20;
lambda_axial = 20;
delta_flag = 1;

estim_method = 1; % 1----- 3 DOF

% estim_method = 2; % 2 ---- 2 DOF
% delta_n = -1.5;

 [Est_n_sample, Est_b_sample, Est_alphaz_sample,delta_SNR,SNR,est_CF,DATA,Xaxis,Zaxis,SR,z,...
    nb_ROI_axial,nb_ROI_lateral,mu_rpl_tv,mu_rpl_robust_tv] = qUS_estimation(Fs, BmodeWidth, phan_alpha, bsc_band,...
    x1, x2, y1,y2, c0, BW, tissue, phantom, line_ROI_lateral, lambda_axial, estim_method, delta_flag, SAM, delta_n);

%%
axis_b = [];
axis_c = [];
axis_a = [-0.18 0];

figure, 
set(gcf,'units','normalized','outerposition',[0 0.25 1 0.65]); box on;
subplot(1,3,1)
imagesc(Xaxis*cm, Zaxis*cm, Est_n_sample), colorbar
m_n=round( mean(Est_n_sample(:)), 2);
std_n=round(std(Est_n_sample(:)), 2);
title(['n: ', num2str(m_n), '+/-', num2str(std_n)]), xlabel('Lateral'), ylabel('Depth'), colormap jet;

subplot(1,3,2)
imagesc(Xaxis*cm, Zaxis*cm, Est_b_sample), colorbar
m_b=round(mean(Est_b_sample(:)), 2);
std_b=round(std(Est_b_sample(:)), 2);
title(['b: ', num2str(m_b), '+/-', num2str(std_b)]), xlabel('Lateral'), ylabel('Depth')

subplot(1,3,3)
imagesc(Xaxis*cm, Zaxis*cm, Est_alphaz_sample, axis_a), colorbar
m_alpha = round( mean(Est_alphaz_sample(:)), 2);
std_alpha= round( std(Est_alphaz_sample(:)), 2);
title(['\alpha: ', num2str(m_alpha), '+/-', num2str(std_alpha)]), xlabel('Lateral'), ylabel('Depth')

