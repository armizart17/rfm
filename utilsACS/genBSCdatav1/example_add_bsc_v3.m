% add bsc_v3 geometrical form

% CREATE THIS GAUSSIAN BSCs

% SAM
% a0.6 g0.01 s0.3
% REF
% a0.4 g1 s0.1

%%
clear all, clc, close all;
addpath(genpath(pwd))
%% LOAD DATA ref_a0.4
% pathData = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\dataACS_kwave';
pathData = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\qus-lim\NewPL_data\data_ok';

name_in = 'ref_a0.4';
SAM = load(fullfile(pathData, name_in+".mat"));

% B-MODE CHECK
bmode_sam = db(hilbert(SAM.rf));
bmode_sam = bmode_sam - max(bmode_sam(:));

figure,
imagesc(SAM.x*1E3, SAM.z*1E3, bmode_sam), axis("image");
xlabel('Lateral [mm]'), ylabel('Depth [mm]');
title('SAM')
colormap('gray')

% g           = 1;
% s           = 0.001;
g           = 1;
s           = -0.2;
geometry    = 'homo';
bsc_model   = 'gauss';
rf_input    = SAM.rf;
xnew        = SAM.x;
znew        = SAM.z;
fs          = SAM.fs;

rfnew = add_bsc_v3(rf_input, fs, geometry, bsc_model, g, s);


%keyboard
c = 1540; dt = 1/fs;
zbsc = c*(1:length(rfnew))*dt/2;

% B-MODE CHECK
bmode_bsc = db(hilbert(rfnew));
bmode_bsc = bmode_bsc - max(bmode_bsc(:));

figure,

imagesc(SAM.x*1E3, zbsc*1E3, bmode_bsc), axis("image");
xlabel('Lateral [mm]'), ylabel('Depth [mm]');
title('BSC')
colormap('gray')

%% % SAVE MODIFIED
name_out = sprintf('%sg%.2fs%.4f', name_in, g, s);
rf = rfnew;
x = SAM.x;
z = SAM.z;
save(fullfile(pathData,name_out+".mat"), 'g', 's', 'rf', 'fs', 'x', 'z')


%% LOAD DATA ref_a0.6
% pathData = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\dataACS_kwave';
pathData = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\qus-lim\NewPL_data\data_ok';

name_in = 'ref_a0.6';
SAM = load(fullfile(pathData, name_in+".mat"));

% B-MODE CHECK
bmode_sam = db(hilbert(SAM.rf));
bmode_sam = bmode_sam - max(bmode_sam(:));

figure,
imagesc(SAM.x*1E3, SAM.z*1E3, bmode_sam), axis("image");
xlabel('Lateral [mm]'), ylabel('Depth [mm]');
title('SAM')
colormap('gray')

g           = 0.01;
s           = 0.;
% g           = 1;
% s           = -0.2;
geometry    = 'homo';
bsc_model   = 'gauss';
rf_input    = SAM.rf;
xnew        = SAM.x;
xnew        = SAM.z;
fs          = SAM.fs;

% g           = 0.01;
% s           = 1.5;
% geometry    = 'homo';
% bsc_model   = 'power-law';
% rf_input    = SAM.rf;
% xnew        = SAM.x;
% xnew        = SAM.z;
% fs          = SAM.fs;

rfnew = add_bsc_v3(rf_input, fs, geometry, bsc_model, g, s);


%keyboard
c = 1540; dt = 1/fs;
zbsc = c*(1:length(rfnew))*dt/2;

% B-MODE CHECK
bmode_bsc = db(hilbert(rfnew));
bmode_bsc = bmode_bsc - max(bmode_bsc(:));

figure,

imagesc(SAM.x*1E3, zbsc*1E3, bmode_bsc), axis("image");
xlabel('Lateral [mm]'), ylabel('Depth [mm]');
title('BSC')
colormap('gray')


%% % SAVE MODIFIED
name_out = sprintf('%sg%.2fs%.4f', name_in, g, s);
rf = rfnew;
x = SAM.x;
z = SAM.z;
save(fullfile(pathData,name_out+".mat"), 'g', 's', 'rf', 'fs', 'x', 'z')


