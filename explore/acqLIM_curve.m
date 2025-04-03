
clear; clc; close all;
fileName = 'Q:\CAvendanoData\42460445\42460445_IHR_F\42460445_IHR_F.mat';
preSetName = 'Q:\CAvendanoData\42460445\42460445_IHR_F\42460445_IHR_F_preSet.mat';


ax = axes(figure);
% generateBMode(fileName, preSetName, ax)

loadAndBfCurve_C52v(fileName, preSetName, ax)

%%

inputDir = 'Q:\CAvendanoData\';
outputDir0 = 'D:\emirandaz\qus\data\liver';

folderOut = 'bf_M04_D02';

outputDir = fullfile(outputDir0, folderOut);

if ~exist(outputDir) mkdir(outputDir); end

folders = dir(inputDir);
folders = folders(3:6); 

for iFolder = 1:length(folders)

    folderStr = folders(iFolder).name;
    subfolders = dir(fullfile(inputDir,folderStr));
    subfolders = subfolders(3:end); % omit '.', '..'
    % isSubfolder = [subdirs.isdir] & ~ismember({subdirs.name}, {'.', '..'});
    
    for isubFolder = 1:length(subfolders)
        subFolderStr = subfolders(isubFolder).name;
        acqDir = dir(fullfile(inputDir,folderStr,subFolderStr));
        acqDir = acqDir(3:end); % omit '.', '..'
        numAcq = fix(length(acqDir)/2); % .mat + preset.mat so only half acq

        for iAcq = 1:numAcq
            samName = acqDir(iAcq).name;
            fileName = fullfile(inputDir, folderStr, subFolderStr, samName);
            presetName = fileName(1:end-4) + "_preSet.mat";
            
            if ~exist(presetName,'file'), continue; end
        
            [rf,bMode,xp,zp,z0p,xr,zr,r0,fs] = loadAndBfCurve_C52v(fileName, presetName);
            
            figure('Units','normalized','Position',[0.1 0.1 0.7 0.6]);

            caption = strrep(samName(1:end-4), '_', ' ');
            fontSize = 12;
            
            t = tiledlayout(1, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');
            
            % First plot (left side)
            nexttile
            imagesc(xr, zr*1e3, bMode);
            cbar = colorbar;
            ylabel(cbar, '[dB]', 'FontSize', 12);
            clim([-55 0]);    
            colormap gray;
            ylabel('Depth [mm]', 'FontSize', 14);
            xlabel('\theta [°]', 'FontSize', 14);
            axis image 
            % axis tight
            title(caption, 'FontSize', fontSize);
            
            % Second plot (right side)
            nexttile
            pcolor(xp*1e3, zp*1e3, bMode);
            shading interp
            cbar = colorbar;
            ylabel(cbar, '[dB]', 'FontSize', 12);
            clim([-55 0]);
            colormap gray               
            ylabel('Depth [mm]', 'FontSize', 14);
            xlabel('Lateral [mm]', 'FontSize', 14);
            axis equal ij tight
            title(caption, 'FontSize', fontSize);
            set(gca, 'Color', 'k');  % axes background

            % Save Figure
            saveas(gcf, fullfile(outputDir,samName(1:end-4)+".png"));
            pause(0.5), close
            
            % Save bf data
            save(fullfile(outputDir,samName(1:end-4)+".mat"), ...
                'rf','bMode','xp','zp','z0p','xr','zr','r0','fs');
 
        end
    end
end

