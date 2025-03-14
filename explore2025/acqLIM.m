baseDir = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\SavedDataQUSPhantom';
acqNames = ["","_2","_3"];
deadBand = 0.1e-2;

folders = dir(baseDir);
% folders = folders(3:end);
% folders = folders(10); % only homogeneous
folders = folders(5:12);

for iFolder = 1:length(folders)
    folderStr = folders(iFolder).name;
    subFolderStr = folderStr + "_F";

    for iAcq = 1:length(acqNames)
        samName = subFolderStr + acqNames(iAcq);
        fileName = fullfile(baseDir, folderStr, subFolderStr, samName);
        presetName = fileName + "_preSet.mat";

        if ~exist(presetName,'file'), continue; end
        [rf,x,z,fs] = loadAndBfLinear(fileName, presetName);

        bMode = db(hilbert(rf));
        bMode = bMode - max(bMode(z>deadBand,:),[],'all');
        figure,
        imagesc(x*1e2, z*1e2, bMode);
        axis image
        colormap gray;
        clim([-60 0]);
        colorbar
        ylabel('[mm]', 'FontSize', 10);
        xlabel('[mm]', 'FontSize', 10);
        ylim([deadBand*100 5])


        resultsDir = fullfile(baseDir, 'bf');
        mkdir(resultsDir)
        saveas(gcf, fullfile(resultsDir,samName+'.png'))
        pause(1), close,
        save(fullfile(resultsDir,samName),'rf','x','z','fs')
    end
end