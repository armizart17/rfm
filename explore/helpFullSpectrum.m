% Simple code to plot Full Spectrum in RFM
% Also identifies and prints BW according to a specif ratio_dB

ratio_dB = -20; % PW
% 40V
% ratio_dB = -30; % PWC 21Ang


spectralData_sam.Sfull      = Sp_k;
spectralData_sam.bandFull   = bandFull;
spectralData_sam.depth      = z_ACS;


[m, n, ~] = size(spectralData_sam.Sfull);

nLines = n;
lin_cen = round(n / 2); 
lat_range = max(1, lin_cen-fix(nLines/2)):min(n, lin_cen+fix(nLines/2)); 

S_2d = squeeze(mean(spectralData_sam.Sfull(:, lat_range, :), 2)); % Mean over 2nd dim (lateral)
S_2d_dB  = pow2db(S_2d ./ max(S_2d, [], 2));

[fL_ratio,fH_ratio] = findFreqBand(bandFull, S_2d(end, :), db2pow(ratio_dB));
fprintf('Usable Freq Range (%.fdB decay from max) is: [%.2f, %.2f] MHz\n', ratio_dB, fL_ratio, fH_ratio)

figure;
set(gcf,'units','normalized','outerposition',[0 0.1 0.5 0.5]); box on;
sgtitle(caption, 'FontWeight', 'Bold');

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
xline(pars.bw(1), 'k--', 'DisplayName', '')
xline(pars.bw(2), 'k--', 'DisplayName', '')
xlim([0 SAM.fs/2]*1e-6); % visualize only positive size
hold off;
xlabel('Frequency [MHz]');
ylabel('Norm Power Spectrum [dB]');
title('SAM Norm Power Spectrum');
legend('Location', 'Best');