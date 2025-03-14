%% SIMPLE COD TO GENERATE .mat files from old SONIX TOUCH
% Sonix Touch acq RF -> .mat
%% HOMOGENEOUS ID261

baseDir = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\sonixtouch\ID261V2\06-08-2023-Generic';
resultsDir = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\sonixtouch\ID261V2\';

frame_ini = 10;
frame_end = 20;

% files = dir(baseDir);
files = dir(fullfile(baseDir,'*.rf'));

for iFiles = 1:length(files)

    filestr = files(iFiles).name(1:end-3); % omit .rf

    fileName = fullfile(baseDir, filestr);

    DATA = lectura_OK([fileName , '.rf']);

    rf    = DATA.RF(:,:,frame_ini:frame_end);
    Bmode = DATA.Bmode(:,:,frame_ini:frame_end);
    fs    = DATA.fs;
    fc    = DATA.fc;
    x     = DATA.x; z = DATA.z;

    figure,
    imagesc(x*1e2, z*1e2, Bmode(:,:,1));
    axis image
    colormap gray;
    clim([-60 0]);
    colorbar
    title(filestr)
    ylabel('Depth [cm]', 'FontSize', 10);
    xlabel('Lateral [cm]', 'FontSize', 10);

    % Save data
    save(fullfile(resultsDir, filestr + ".mat"), 'rf','x','z','fs','fc','Bmode')

    % Save referential bmode image
    saveas(gcf, fullfile(resultsDir, filestr +".png"))         
        
end

%% HOMOGENEOUS ID544

baseDir = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\sonixtouch\ID544V2\06-08-2023-Generic';
resultsDir = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\sonixtouch\ID544V2\';

frame_ini = 10;
frame_end = 20;

% files = dir(baseDir);
files = dir(fullfile(baseDir,'*.rf'));

for iFiles = 1:length(files)

    filestr = files(iFiles).name(1:end-3); % omit .rf

    fileName = fullfile(baseDir, filestr);

    DATA = lectura_OK([fileName , '.rf']);

    rf    = DATA.RF(:,:,frame_ini:frame_end);
    Bmode = DATA.Bmode(:,:,frame_ini:frame_end);
    fs    = DATA.fs;
    fc    = DATA.fc;
    x     = DATA.x; z = DATA.z;

    figure,
    imagesc(x*1e2, z*1e2, Bmode(:,:,1));
    axis image
    colormap gray;
    clim([-60 0]);
    colorbar
    title(filestr)
    ylabel('Depth [cm]', 'FontSize', 10);
    xlabel('Lateral [cm]', 'FontSize', 10);

    % Save data
    save(fullfile(resultsDir, filestr + ".mat"), 'rf','x','z','fs','fc','Bmode')

    % Save referential bmode image
    saveas(gcf, fullfile(resultsDir, filestr +".png"))         
        
end