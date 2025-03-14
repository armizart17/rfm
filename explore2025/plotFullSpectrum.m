%% PLOT FULL Spectrum through Depth

spectralData_sam = calc_powerSpectra_vSimple(SAM, pars);
% spectralData_sam = calc_powerSpectraFull_prox_dis(SAM, pars);

ratio_dB = -20;
ratio = db2mag(ratio_dB);

% Get size of Sp

Sfull = spectralData_sam.Sfull; % POWER LAW 3DOF

[m, n, ~] = size(spectralData_sam.Sfull);

nLines = 5;
lin_cen = round(n / 2); 
lat_range = max(1, lin_cen-fix(nLines/2)):min(n, lin_cen+fix(nLines/2)); 

S_2d = squeeze(mean(spectralData_sam.Sfull(:, lat_range, :), 2)); % Mean over 2nd dim (lateral)
S_2d_dB  = pow2db(S_2d ./ max(S_2d, [], 2));

figure;
set(gcf,'units','normalized','outerposition',[0 0.1 0.5 0.5]); box on;

subplot(1,2,1)
imagesc(spectralData_sam.bandFull, spectralData_sam.depth*1e3, S_2d_dB),
xline(pars.bw(1), 'w--', 'LineWidth', 2)
xline(pars.bw(2), 'w--', 'LineWidth', 2)
xlim([0 SAM.fs/2]*1e-6); % visualize only positive size
xlabel('Frequency [MHz]');
ylabel('Depth [mm]');
h2 = colorbar; 
ylabel(h2,'dB');
title('SAM Norm Power Spectrum by depth');

subplot(1,2,2)
plot(spectralData_sam.bandFull, S_2d_dB(1, :), 'DisplayName', 'Top')
hold on, grid on;
plot(spectralData_sam.bandFull, S_2d_dB(round(m/2), :), 'DisplayName', 'Half')
plot(spectralData_sam.bandFull, S_2d_dB(end, :), 'DisplayName', 'Bottom')
yline(ratio_dB, 'k--', 'DisplayName', '')
xlim([0 SAM.fs/2]*1e-6); % visualize only positive size
hold off;
xlabel('Frequency [MHz]');
ylabel('Norm Power Spectrum [dB]');
title('SAM Norm Power Spectrum');
legend('Location', 'Best');

% FIND PEAKS at all depths

mDepths = size(S_2d, 1);
f = spectralData_sam.bandFull;
arrayBands = zeros(mDepths,  2);
for dd = 1:mDepths
    y = S_2d(dd, :);
    [fLeft,fRight] = findFreqBand(f, y, ratio);
    arrayBands(dd,1) = fLeft;
    arrayBands(dd,2) = fRight;
end

figure,
imagesc(f, spectralData_sam.depth*1e3, S_2d_dB);
hold on
plot(arrayBands(:,1), spectralData_sam.depth * 1e3, 'r', 'LineWidth', 2); % Left boundary
plot(arrayBands(:,2), spectralData_sam.depth * 1e3, 'r', 'LineWidth', 2); % Right boundar
hold off
xlabel('Frequency [MHz]');
xlim([0 SAM.fs/2]*1e-6); % visualize only positive size
ylabel('Depth [mm]');
h2 = colorbar; 
ylabel(h2,'dB');
title('Norm Power Spectrum');

