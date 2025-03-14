% add bsc_v2 geometrical form

%%
clear all, clc, close all;

addpath(genpath(pwd))
%% LOAD DATA SAM0
pathData = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\dataACS_kwave';

SAM = load(fullfile(pathData, 'sam0.mat'));

% B-MODE CHECK
bmode_sam = db(hilbert(SAM.rf));
bmode_sam = bmode_sam - max(bmode_sam(:));

figure,
imagesc(SAM.x*1E3, SAM.z*1E3, bmode_sam), axis("image");
xlabel('Lateral [mm]'), ylabel('Depth [mm]');
title('SAM')
colormap('gray')


b = 1e-2;
n = 1.5;
rf = SAM.rf;
fs = SAM.fs;
b_vals = [1, 1E-2];
b_vals = [1-100, 1E2];
n_vals = [1.5 1.5];
% n_vals = [1 1];
% geometry = 'two_layer';
% rfbsc = add_bsc_v2(rf, b_vals, n_vals, fs, geometry);

geometry = 'circle';
mask = SAM.attenuation_map;
mask(mask<1) = 0;

mask = imresize(mask, size(SAM.rf), 'nearest');
mask = mask>0.5;


rfbsc = add_bsc_v2(rf, b_vals, n_vals, fs, geometry, mask);


%keyboard

c = 1540; dt = 1/fs;
zbsc = c*(1:length(rfbsc))*dt/2;

% B-MODE CHECK
bmode_bsc = db(hilbert(rfbsc));
bmode_bsc = bmode_bsc - max(bmode_bsc(:));

figure,

imagesc(SAM.x*1E3, zbsc*1E3, bmode_bsc), axis("image");
xlabel('Lateral [mm]'), ylabel('Depth [mm]');
title('BSC')
colormap('gray')

%% SAV SAMO MODIFIED
name = 'sam0_vertLayer.mat';
name = 'sam0_vertLayer.mat';
rf = rfbsc;
x = SAM.x;
z = SAM.z;
save(name, 'b_vals', 'n_vals', 'rf', 'fs', 'x', 'z')
%%




%% LOAD REF REF
pathData = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\dataACS_kwave';

SAM = load(fullfile(pathData, 'ref0.mat'));

% B-MODE CHECK
bmode_sam = db(hilbert(SAM.rf));
bmode_sam = bmode_sam - max(bmode_sam(:));

figure,
imagesc(SAM.x*1E3, SAM.z*1E3, bmode_sam), axis("image");
xlabel('Lateral [mm]'), ylabel('Depth [mm]');
title('SAM')
colormap('gray')


b = 1e-2;
n = 1.5;
rf = SAM.rf;
fs = SAM.fs;
b_vals = [1E-2, 1E-2];
n_vals = [1.5 1.5];
% n_vals = [1.5 1];
geometry = 'two_layer';

rfbsc = add_bsc_v2(rf, b_vals, n_vals, fs, geometry);


%keyboard

c = 1540; dt = 1/fs;
zbsc = c*(1:length(rfbsc))*dt/2;

% B-MODE CHECK
bmode_bsc = db(hilbert(rfbsc));
bmode_bsc = bmode_bsc - max(bmode_bsc(:));

figure,

imagesc(SAM.x*1E3, zbsc*1E3, bmode_bsc), axis("image");
xlabel('Lateral [mm]'), ylabel('Depth [mm]');
title('BSC')
colormap('gray')

%% SAV SAMO MODIFIED
name = 'ref0_mod.mat';
rf = rfbsc;
x = SAM.x;
z = SAM.z;
save(name, 'b_vals', 'n_vals', 'rf', 'fs', 'x', 'z')