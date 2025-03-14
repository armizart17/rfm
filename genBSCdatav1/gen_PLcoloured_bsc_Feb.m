% add bsc_v3 geometrical form

% CREATE COLOURED BSC POWER LAW FEB2025

% SAM
% a0.6 g0.01 s0.3
% REF
% a0.4 g1 s0.1

%%
clear all, clc, close all;
addpath(genpath(pwd))
%% %%%%%%%%%%%% LOAD DATA ref_a0.6 (b=0.01, n=1.5)

pathData = ['C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\' ...
    'dataACS_kwave\simFeb2025\PowerLaw1p\a0p6_pow1p\'];

% Get list of all .mat files in the folder
matFiles = dir(fullfile(pathData, '*.mat'));

% BSC specs
b           = 0.01;
n           = 1.5;
geometry    = 'homo';
bsc_model   = 'power-law';

for i = 1:length(matFiles)
    fileNameIn = matFiles(i).name(1:end-4); % do not add .mat
    SAM = load( fullfile(pathData, fileNameIn ) );

    % B-MODE CHECK
    bmode_sam = db(hilbert(SAM.rf));
    bmode_sam = bmode_sam - max(bmode_sam(:));
    figure,
    imagesc(SAM.x*1E3, SAM.z*1E3, bmode_sam), axis("image");
    xlabel('Lateral [mm]'), ylabel('Depth [mm]');
    title('SAM') , colormap('gray')

    rf_input    = SAM.rf;
    xnew        = SAM.x;
    znew        = SAM.z;
    fs          = SAM.fs;

    rfnew = add_bsc_v3(rf_input, fs, geometry, bsc_model, b, n);

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

    % SAVE NEW RF
    rf = rfnew;
    x = SAM.x;
    z = SAM.z;
    
    fileNameOut = [fileNameIn,'b',num2str(b),'n',num2str(n)];

    save(fullfile(pathData, fileNameOut+".mat"), 'b', 'n', 'rf', 'fs', 'x', 'z')
end

%% %%%%%%%%%%%% LOAD DATA ref_a0.4 (b=0.01, n=1.5)


pathData = ['C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\' ...
    'dataACS_kwave\simFeb2025\PowerLaw1p\a0p4_pow1p\'];

% Get list of all .mat files in the folder
matFiles = dir(fullfile(pathData, '*.mat'));

% BSC specs
b           = 1;
n           = 0;
geometry    = 'homo';
bsc_model   = 'power-law';

for i = 1:length(matFiles)
    fileNameIn = matFiles(i).name(1:end-4); % do not add .mat
    SAM = load( fullfile(pathData, fileNameIn ) );

    % B-MODE CHECK
    bmode_sam = db(hilbert(SAM.rf));
    bmode_sam = bmode_sam - max(bmode_sam(:));
    figure,
    imagesc(SAM.x*1E3, SAM.z*1E3, bmode_sam), axis("image");
    xlabel('Lateral [mm]'), ylabel('Depth [mm]');
    title('SAM') , colormap('gray')

    rf_input    = SAM.rf;
    xnew        = SAM.x;
    znew        = SAM.z;
    fs          = SAM.fs;

    rfnew = add_bsc_v3(rf_input, fs, geometry, bsc_model, b, n);

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

    % SAVE NEW RF
    rf = rfnew;
    x = SAM.x;
    z = SAM.z;
    
    fileNameOut = [fileNameIn,'b',num2str(b),'n',num2str(n)];

    save(fullfile(pathData, fileNameOut+".mat"), 'b', 'n', 'rf', 'fs', 'x', 'z')

end

%%
% GENERATE ONE RF in each plane each case

