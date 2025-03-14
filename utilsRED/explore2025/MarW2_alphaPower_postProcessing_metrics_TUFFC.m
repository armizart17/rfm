
clear;
clc;
close all;

myColormap = 'sky'; % hot, cool, summer, autumn, bone, copper, pink, sky
% myColormap = 'hot';
% myColormap = 'autumn';
font_size = 30;

%% AC METRICS



methods = {'3-DoF', '2-DoF-b', '2-DoF-n'}; % Excluding 2-DoF-a

% Load Data
pathData = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\qus-lim\out\alpha_power\';

fileNameRPLTV = 'results_1p1_1p9';
rpltv = load(fullfile(pathData, fileNameRPLTV));

fileNamefreq = 'freq';
load(fullfile(pathData, fileNamefreq));

% Extract necessary data
maps_results_all = rpltv.maps_results_all;
j_sam_values = rpltv.j_sam_values;
j_ref_values = rpltv.j_ref_values;
alpha_sam = rpltv.alpha_sam;

num_sam = length(j_sam_values);
num_ref = length(j_ref_values);

% Preallocate metrics storage
metrics_all = cell(num_sam, num_ref, length(methods));

% Compute metrics for selected methods
tic
for iMethod = 1:length(methods)
    for iSam = 1:num_sam
        for iRef = 1:num_ref
            % Extract the a parameter (attenuation coefficient map)
            method_index = find(strcmp({'3-DoF', '2-DoF-a', '2-DoF-b', '2-DoF-n'}, methods{iMethod}));
            if ~isempty(maps_results_all{iSam, iRef, method_index}) 
                a_map = maps_results_all{iSam, iRef, method_index}{1}; 

                % Define a mask (since all values are used, logical ones)
                mask_homo = logical(ones(size(a_map)));

                % Compute metrics
                metrics = get_metrics_homo_gt(a_map, mask_homo, alpha_sam, methods{iMethod});

                % Store metrics
                metrics_all{iSam, iRef, iMethod} = metrics;
            end
        end
    end
end
fprintf('AC Elapsed time: %.2f seconds\n', toc)

%% PLOT AC NRMSE HEATMAP FOR SELECTED METHODS (3-DoF, 2-DoF-b, 2-DoF-n)
j_sam_values = 1.1:0.1:1.5;
j_ref_values = j_sam_values;

methods = {sprintf('3-DoF'), sprintf('2-DoF "b"'),  sprintf('2-DoF "n"')};
for iMethod = 1:length(methods)
    metric_name = 'nrmse_homo';  
    metric_matrix = nan(length(j_sam_values), length(j_ref_values));

    % Extract metric values into a matrix
    for iSam = 1:length(j_sam_values)
        for iRef = 1:length(j_ref_values)
            if ~isempty(metrics_all{iSam, iRef, iMethod})
                metric_matrix(iSam, iRef) = metrics_all{iSam, iRef, iMethod}.(metric_name);
            end
        end
    end

    % Plot the heatmap
    % subplot(1,3,iMethod
    figure,
    set(gcf, 'Units', 'pixels', 'Position', [100, 100, 900, 800]); % [x, y, width, height] in pixels
    imagesc(j_sam_values, j_ref_values, metric_matrix .* 100); % Convert to percentage
    axis("image")
    colormap(myColormap);
    % clim([0 70])
    hb2 = colorbar; ylabel(hb2, 'NRMSE (%)', 'FontSize', font_size+2);
    
    title(sprintf('AC: %s', methods{iMethod}), 'FontWeight', 'bold');
    set(gca, 'FontSize', font_size);
    xlabel('$m_s$', 'Interpreter', 'latex', 'FontSize', font_size+15);
    ylabel('$m_r$', 'Interpreter', 'latex', 'FontSize', font_size+15);
    set(gca, 'YDir', 'normal'); % Ensures correct orientation
end


%% BSC METRICS

methods = {'3-DoF', '2-DoF-a', '2-DoF-b', '2-DoF-n'}; % Excluding 2-DoF-a

% Load Data
pathData = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\qus-lim\out\alpha_power\';

fileNameRPM = 'resultsRPM_1p1_1p9';
rpm = load(fullfile(pathData, fileNameRPM));

fileNameRPLTV = 'results_1p1_1p9';
rpltv = load(fullfile(pathData, fileNameRPLTV));

fileNamefreq = 'freq';
load(fullfile(pathData, fileNamefreq));

% Extract necessary data
bsc_results_all = rpltv.bsc_results_all;
j_sam_values = rpltv.j_sam_values;
j_ref_values = rpltv.j_ref_values;
alpha_sam = rpltv.alpha_sam;

bsc_results_all_rpm = rpm.bsc_results_all;

num_sam = length(j_sam_values);
num_ref = length(j_ref_values);

% Preallocate metrics storage
metrics_all_bsc = cell(num_sam, num_ref, length(methods));

% Compute metrics for selected methods
tic
for iMethod = 1:length(methods)
    for iSam = 1:num_sam
        for iRef = 1:num_ref
            % Extract the a parameter (attenuation coefficient map)
            method_index = find(strcmp({'3-DoF', '2-DoF-a', '2-DoF-b', '2-DoF-n'}, methods{iMethod}));
            if ~isempty(bsc_results_all{iSam, iRef, method_index}) 
                bsc_rpl = bsc_results_all{iSam, iRef, method_index}{1}; 

                bsc_rpm = bsc_results_all_rpm{iSam, iRef}{1}; 

                % Define a mask (since all values are used, logical ones)
                mask_homo = logical(ones(size(bsc_rpl)));

                % Compute metrics
                metrics = get_metrics_homo_gt(bsc_rpl, mask_homo, bsc_rpm, methods{iMethod});

                % Store metrics
                metrics_all_bsc{iSam, iRef, iMethod} = metrics;
            end
        end
    end
end
fprintf('BSC Elapsed time: %.2f seconds\n', toc)

%% PLOT NRMSE BSC HEATMAP FOR SELECTED METHODS (3-DoF, 2-DoF-b, 2-DoF-n)
j_sam_values = 1.1:0.1:1.5;
j_ref_values = j_sam_values;

methods = {sprintf('3-DoF'), sprintf('2-DoF "a"'), sprintf('2-DoF "b"'),  sprintf('2-DoF "n"')};

for iMethod = 1:length(methods)
    metric_name = 'nrmse_homo';  
    metric_matrix = nan(length(j_sam_values), length(j_ref_values));

    % Extract metric values into a matrix
    for iSam = 1:length(j_sam_values)
        for iRef = 1:length(j_ref_values)
            if ~isempty(metrics_all_bsc{iSam, iRef, iMethod})
                metric_matrix(iSam, iRef) = metrics_all_bsc{iSam, iRef, iMethod}.(metric_name);
            end
        end
    end

    % Plot the heatmap
    figure,
    set(gcf, 'Units', 'pixels', 'Position', [100, 100, 920, 800]); % [x, y, width, height] in pixels
    imagesc(j_sam_values, j_ref_values, metric_matrix .* 100); % Convert to percentage
    axis("image")
    colormap(myColormap);
    % clim([0 100])
    hb2 = colorbar; ylabel(hb2, 'NRMSE (%)', 'FontSize', font_size+2);
    
    title(sprintf('BSC: %s', methods{iMethod}), 'FontWeight', 'bold');
    set(gca, 'FontSize', font_size);
    xlabel('$m_s$', 'Interpreter', 'latex', 'FontSize', font_size+15);
    ylabel('$m_r$', 'Interpreter', 'latex', 'FontSize', font_size+15);
    set(gca, 'YDir', 'normal'); % Ensures correct orientation
    
end

%% SAVE FIGURES

keyboard
dirFigout = '.\TUFFC25\simuGauss\ac0p5\gaussianmodel\alpha_powerv2';
if (~exist(dirFigout)); mkdir (dirFigout); end
titleFigout = 'Fig';
save_all_figures_to_directory(dirFigout, titleFigout, 'svg')