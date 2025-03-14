% TEST OLD SLD
% ACQ with Sonix Touch, function calc_ProxDis
% Check targets
% TO SLD LINE
% regionMaskAcs = ones(m,n);
% sldLine = squeeze(sum(sum(b.*regionMaskAcs,2),1))/4/L*Np2dB/sum(regionMaskAcs(:));
% sldLine = squeeze(sum(sum(b.*regionMaskAcs,2),1)) / sum(regionMaskAcs(:)) / (4*L)*Np2dB;
% fit1 = f\sldLine; % slope = fit1
% fit2 = [f ones(length(f),1)]\sldLine; % slope=acs=fit2(1) intercept=fit2(1)
% cgs(f'*f, f'*sldLine)
%%
%% ALL PHANTOMS
% clear, clc;
% close all

Np2dB       = 20*log10(exp(1));
dB2Np       = 1/Np2dB;

% Weight parameters
ratioCutOff = 10;
reject = 0.1;
extension = 3;

% SWTV
aSNR = 5; bSNR = 0.09;
desvMin = 15;

% Plotting constants
dynRange = [-50,0];
attRange = [0.2,1.2];

tol = 1e-3;

c1x = 1.95; c1z = 1.93;
roiL = 1; roiD = 0.6;
roiLz = 1.5;

addpath(genpath(pwd));

dataDir = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\phantomTargets\SonixTouch\ID316';
refDir = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\phantomTargets\SonixTouch\refPh544\';

% resultsDir = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\phantomTargets\SonixTouch\ID316_out';

rawFiles = dir([dataDir,'\*.rf']);
targetFiles = dir([dataDir,'\*.mat']);

% targetFiles = targetFiles(end-2:end); 
num_targets = length(targetFiles);

groundTruthTargets    = [0.52, 0.55, 0.74, 0.81, 0.75, 0.97, 0.95, 0.95]; % [dB/cm/MHz]
groundTruthBackground = 0.55; % [dB/cm/MHz]

for iAcq = 1:num_targets
    %% SWITCH EMZ
switch iAcq
    % Optimal reg for BS 8x12, circular ROIs
    case 6
        muBtv = 10^3.5; muCtv = 10^3.5;
        muBswtv = 10^3; muCswtv = 10^2.5;
        muBtvl1 = 10^3.5; muCtvl1 = 10^2;
        muBswift = 10^3.5; muCswift = 10^2;
    case 7
        muBtv = 10^3; muCtv = 10^3;
        muBswtv = 10^3; muCswtv = 10^0;
        muBtvl1 = 10^3.5; muCtvl1 = 10^1;
        muBswift = 10^3.5; muCswift = 10^1;
    case 8
        muBtv = 10^3.5; muCtv = 10^3.5;
        muBswtv = 10^3; muCswtv = 10^0;
        muBtvl1 = 10^3.5; muCtvl1 = 10^1;
        muBswift = 10^3.5; muCswift = 10^1;
    otherwise
        muBtv = 10^3.5; muCtv = 10^3.5;
        muBswtv = 10^3; muCswtv = 10^0;
        muBtvl1 = 10^3.5; muCtvl1 = 10^1;
        muBswift = 10^3.5; muCswift = 10^1;
end

switch iAcq
    case 6
        c1x = 1.8; c1z = 1.9;
    case 7
        c1x = 1.95; c1z = 1.95;
    case 8
        c1x = 1.85; c1z = 1.9;
    otherwise
        c1x = 1.85; c1z = 1.9;
end
%%
fprintf("Phantom no. %i, %s\n",iAcq,targetFiles(iAcq).name);
SAM = load(fullfile(dataDir,targetFiles(iAcq).name));
SAM.rf = SAM.RF(:,:,1);
SAM.RF = SAM.RF(:,:,1);

pars.bw             = [2.5 7.5];
pars.overlap        = 0.8;
pars.blocksize      = 8;
pars.x_roi          = [0.1 3.8]*1e-2;
pars.z_roi          = [0.2 3.5]*1e-2;
pars.window_type    = 3;
pars.saran_layer    = 0;
pars.ratio_zx       = 1.5;


refFiles = dir([refDir,'\*.mat']);
numRefs = length(refFiles);

REF     = load([refDir,'\',refFiles(1).name]);
newrf  = nan([size(REF.RF), numRefs], 'like', REF.RF(:,:,1)); % Use 'like' for type consistency
for i = 1:numRefs
    newrf(:,:,i) = load([refDir,'\',refFiles(i).name],'RF').RF(:,:,1);
end
REF.rf = newrf;
alpha_ref = 0.53;
%% PLOT BMODE

% SAMPLE
spectralData_sam = calc_powerSpectra_prox_dis(SAM, pars);

Sp_sam = spectralData_sam.Sp;
Sd_sam = spectralData_sam.Sd;

% REFERENCE
num_ref = 1;
spectralData_ref = calc_powerSpectra_prox_dis(REF, pars);

band   = spectralData_sam.band; % [MHz]
zd_zp  = spectralData_sam.zd_zp; % [m]
Sp_ref = spectralData_ref.Sp;
Sd_ref = spectralData_ref.Sd;

% COMPENSATION ATTENUATION
acs_ref     = alpha_ref; % [dB/cm/MHz]
att_ref     = acs_ref*band/Np2dB; % vector [Np/cm]
att_ref_map = ones(size(Sp_sam)) .* reshape(att_ref, 1, 1, []); % 3D array as SLogRatio

SLogRatio = log(Sp_sam./Sd_sam) - ( log(Sp_ref./Sd_ref) - 4*zd_zp*1E2*att_ref_map ); % (rows, cols, freqChannels)

%% ROI selection

[X,Z] = meshgrid(spectralData_sam.lateral,spectralData_sam.depth);
[Xq,Zq] = meshgrid(SAM.x,SAM.z);
rInc = 0.95;
inc = ((Xq-c1x).^2 + (Zq-c1z).^2)<= (rInc-0.1)^2;
back = ((Xq-c1x).^2 + (Zq-c1z).^2) >= (rInc+0.1)^2;

x0mask = c1x - roiL/2; 
z0mask = c1z - roiLz/2;
% [back,inc] = getRegionMasks(x,z,c1x,c1z,roiL,roiD,roiLz);

% Setting up
[m, n, p] = size(SLogRatio);
b = SLogRatio;
A1 = kron( 4*zd_zp*1E2*band , speye(m*n) );
A2 = kron( ones(size(band)) , speye(m*n) );
mask = ones(m,n,p);

%% RSLD-TV

tic
[Bn,Cn,ite] = AlterOpti_ADMM(A1,A2,b(:),muBtv,muCtv,m,n,tol,mask(:));
exTime = toc;
BR = (reshape(Bn*Np2dB,m,n));
CR = (reshape(Cn,m,n));

AttInterp = interp2(X,Z,BR,Xq,Zq);
r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdInc = std(AttInterp(inc),"omitnan");
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.biasBack = mean( AttInterp(back) - groundTruthBackground,"omitnan");
r.biasInc = mean( AttInterp(inc) - groundTruthTargets(iAcq),"omitnan");
r.rmseBack = sqrt( mean( (AttInterp(back) - groundTruthBackground).^2,...
    "omitnan") );
r.rmseInc = sqrt( mean( (AttInterp(inc) - groundTruthTargets(iAcq)).^2,...
    "omitnan") );
r.cnr = abs(r.meanBack - r.meanInc)/sqrt(r.stdInc^2 + r.stdBack^2);
r.method = 'TV';
r.sample = iAcq;
r.ite = ite;
r.exTime = exTime;
MetricsTV(iAcq) = r;

%% British Columbia Approach

% Weights

desvSNR = spectralData_sam.delta_snr;
wSNR = aSNR./(1 + exp(bSNR.*(desvSNR - desvMin)));

% computation
tic
[Bn,Cn,ite] = AlterOptiAdmmAnisWeighted(A1,A2,b(:),muBswtv,muCswtv, ...
    m,n,tol,mask(:),wSNR);
exTime = toc;
BRBC = (reshape(Bn*Np2dB,m,n));
CRBC = (reshape(Cn,m,n));

AttInterp = interp2(X,Z,BRBC,Xq,Zq);
r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdInc = std(AttInterp(inc),"omitnan");
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.biasBack = mean( AttInterp(back) - groundTruthBackground,"omitnan");
r.biasInc = mean( AttInterp(inc) - groundTruthTargets(iAcq),"omitnan");
r.rmseBack = sqrt( mean( (AttInterp(back) - groundTruthBackground).^2,...
    "omitnan") );
r.rmseInc = sqrt( mean( (AttInterp(inc) - groundTruthTargets(iAcq)).^2,...
    "omitnan") );
r.cnr = abs(r.meanBack - r.meanInc)/sqrt(r.stdInc^2 + r.stdBack^2);
r.method = 'SWTV';
r.sample = iAcq;
r.ite = ite;
r.exTime = exTime;
MetricsSWTV(iAcq) = r;

%% SWIFT
% First iteration
[~,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muBswift,muCswift,m,n,tol,mask(:));
bscMap = reshape(Cn*Np2dB,m,n);

% Weight map
w = (1-reject)*(abs(bscMap)<ratioCutOff)+reject;
wExt = movmin(w,extension);

% Weight matrices and new system
W = repmat(wExt,[1 1 p]);
W = spdiags(W(:),0,m*n*p,m*n*p);
bw = W*b(:);        
A1w = W*A1;
A2w = W*A2;

% Second iteration
tic
[Bn,~,ite] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,muBswift,muCswift,m,n,tol,mask(:),w);
exTime = toc;
BSWIFT = reshape(Bn*Np2dB,m,n);

AttInterp = interp2(X,Z,BSWIFT,Xq,Zq);
r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdInc = std(AttInterp(inc),"omitnan");
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.biasBack = mean( AttInterp(back) - groundTruthBackground,"omitnan");
r.biasInc = mean( AttInterp(inc) - groundTruthTargets(iAcq),"omitnan");
r.rmseBack = sqrt( mean( (AttInterp(back) - groundTruthBackground).^2,...
    "omitnan") );
r.rmseInc = sqrt( mean( (AttInterp(inc) - groundTruthTargets(iAcq)).^2,...
    "omitnan") );
r.cnr = abs(r.meanBack - r.meanInc)/sqrt(r.stdInc^2 + r.stdBack^2);
r.method = 'SWIFT';
r.sample = iAcq;
r.ite = ite;
r.exTime = exTime;
MetricsSWIFT(iAcq) = r;

%%
x = SAM.x*1e2;
z = SAM.z*1e2;
x_ACS = spectralData_sam.lateral*1e2;
z_ACS = spectralData_sam.depth*1e2;
Bmode = db(hilbert(spectralData_sam.rf_roi));
Bmode = Bmode - max(Bmode(:));

figure('Units','centimeters', 'Position',[5 2 6 16]);
tl = tiledlayout(4,1, "Padding","tight");

t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
% xlabel('Lateral [cm]'),
ylabel('Axial [cm]')
axis image
colormap(t1,gray)
title('B-mode')
%subtitle(' ')
c = colorbar;
c.Label.String = 'dB';
% hold on 
% rectangle('Position',[c1x-rInc c1z-rInc 2*rInc 2*rInc], 'LineStyle','--', ...
%     'LineWidth',1, 'Curvature',1)
% hold off

% fontsize(gcf,8,'points')

t2 = nexttile;
imagesc(x_ACS,z_ACS,BR, attRange)
% xlabel('Lateral [cm]'),
ylabel('Axial [cm]')
colormap(t2,turbo)
axis image
title('RSLD')
c = colorbar;
c.Label.String = 'ACS [dB/cm/MHz]';
% fontsize(gcf,8,'points')
hold on 
rectangle('Position',[c1x-rInc c1z-rInc 2*rInc 2*rInc], 'LineStyle','--', ...
    'LineWidth',1, 'Curvature',1, 'EdgeColor','w')
hold off

t3 = nexttile;
imagesc(x_ACS,z_ACS,BRBC, attRange)
% xlabel('Lateral [cm]'),
ylabel('Axial [cm]')
colormap(t3,turbo)
axis image
title('SWTV-ACE')
c = colorbar;
c.Label.String = 'ACS [dB/cm/MHz]';
% fontsize(gcf,8,'points')
hold on 
rectangle('Position',[c1x-rInc c1z-rInc 2*rInc 2*rInc], 'LineStyle','--', ...
    'LineWidth',1, 'Curvature',1, 'EdgeColor','w')
hold off

t5 = nexttile;
imagesc(x_ACS,z_ACS,BSWIFT, attRange)
xlabel('Lateral [cm]'),
ylabel('Axial [cm]')
colormap(t5,turbo)
axis image
title('SWIFT')
c = colorbar;
c.Label.String = 'ACS [dB/cm/MHz]';

% fontsize(gcf,8,'points')
hold on 
rectangle('Position',[c1x-rInc c1z-rInc 2*rInc 2*rInc], 'LineStyle','--', ...
    'LineWidth',1, 'Curvature',1, 'EdgeColor','w')
hold off
end

%%

keyboard
%%
regionMaskAcs = ones(m,n);
sldLine = squeeze(sum(sum(b.*regionMaskAcs,2),1))/4/L*Np2dB/sum(regionMaskAcs(:));
sldLine = squeeze(sum(sum(b.*regionMaskAcs,2),1)) / sum(regionMaskAcs(:)) / (4*L)*Np2dB;

fit1 = f\sldLine;
fit2 = [f ones(length(f),1)]\sldLine;

cgs(f'*f, f'*sldLine)
figure('Units','centimeters', 'Position',[5 5 20 10]),
tiledlayout(1,2),
nexttile
plot(f,sldLine),
hold on,
plot(f,fit1*f, '--')
plot(f,fit2(1)*f + fit2(2), '--')
hold off,
grid on,
xlim([0,freq_H*1.1]/1e6),
ylim([-1 5]),
xlabel('Frequency [MHz]')
ylabel('Attenuation [dB/cm]')
title('Mean SLD')
legend('Liver',sprintf('ACS ZC = %.2f\n',fit1), ...
    sprintf('ACS NZC = %.2f\n',fit2(1)), 'Location','northwest')

sldLine = squeeze(sum(sum(b .* regionMaskAcs, 2), 1))  / sum(regionMaskAcs(:));


sldLine2 = squeeze(mean(mean(b .* regionMaskAcs, 2), 1)) / 4 / L * Np2dB;
sldLine3 = squeeze(mean(mean(b .* regionMaskAcs, 1), 2)) / 4 / L * Np2dB;

%% FIND PEAKS at all depths

mDepths = size(Sd_Pos, 1);
f = freq_Pos;
ratio = db2mag(-20);
arrayBands = zeros(mDepths,  2);
for dd = 1:mDepths
    y = Sd_Pos(dd, :);
    [fLeft,fRight] = findFreqBand(f, y, ratio);
    arrayBands(dd,1) = fLeft;
    arrayBands(dd,2) = fRight;
end
%% FIGURE OF BW at all depths with peaks identified
figure,
imagesc(freq_Pos, spectralData_sam.depth*1e3, Sd_Pos_dB);
hold on
plot(arrayBands(:,1), spectralData_sam.depth * 1e3, 'r', 'LineWidth', 2); % Left boundary
plot(arrayBands(:,2), spectralData_sam.depth * 1e3, 'r', 'LineWidth', 2); % Right boundar
hold off
xlabel('Frequency [MHz]');
ylabel('Depth [mm]');
h2 = colorbar; 
ylabel(h2,'dB');
title('Norm Power Spectrum');


%%
%%
% pathOut = pwd;
resultsDir = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\SavedDataQUSPhantom\bf\targetOut';
if (~exist(resultsDir)) mkdir(resultsDir); end
save_all_figures_to_directory(resultsDir,'targets','png'); % @
% close all
