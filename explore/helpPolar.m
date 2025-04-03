figure, 
pcolor(xp*1e3, zp*1e3, ...
    bModeImage(:, :,end))
colorbar
clim([-55, 0]);
colormap gray
%%
% ax = axes(figure);
figure, 
pcolor(xp*1e3, zp*1e3, bModeImage(:, :,end));
% axes.Color = [0 0 0];
% axes.FontSize = 13;
cbar = colorbar;
ylabel(cbar, '[dB]');
% cbar.Ticks = [];
clim([-55, 0]);
colormap gray
fontSize = 18;

% title(caption, 'FontSize', fontSize);
ylabel('[mm]', 'FontSize', 12);
xlabel('[mm]', 'FontSize', 12);
% shading("interp")
shading("flat")
axis equal ij tight

%%
z2 = (0:size(bModeImage, 1) - 1)/sampleFreq*soundSpeed/2; % [m]
    
[xp_cs, zp_cs] = toPolarGrid(size(bModeImage, [1, 2]), z2(end),param);

% function [xPolar, zPolar] = toPolarGrid(siz,zmax,param)                    
%     N = param.Nelements;
%     R = param.radius;
%     p = param.pitch;
%     L = 2*R*sin(asin(p/2/R)*(N-1)); % chord length
%     d = sqrt(R^2-L^2/4); % apothem
%     z0 = -d;
%     th      = -(linspace(atan2(L/2,d),atan2(-L/2,d),siz(2)))*180/pi;
%     r       = linspace(R+p,-z0+zmax,siz(1));
% 
%     [thth,rr] = meshgrid(...
%         linspace(atan2(L/2,d),atan2(-L/2,d),siz(2))+pi/2,...
%         linspace(R+p,-z0+zmax,siz(1)));
%     [xPolar,zPolar] = pol2cart(thth,rr);
%     zPolar = zPolar+z0;
% end    

%%


[xPolar,zPolar, z0Polar] = impolgrid(size(bModeImage), z(end),param);

xFull = th; % [deg]
r0    = r(1);
zFull = (r-r0)*1e2; % [cm]

%%

[TH_acs,R_acs]        = meshgrid(-x_ACS*pi/180 + pi/2,z_ACS/100 + r0);
[xPolarACS,zPolarACS] = pol2cart(TH_acs,R_acs);
zPolarACS = zPolarACS + z0Polar;

BmodeFull = db(hilbert(rf));
BmodeFull = BmodeFull - max(BmodeFull(:));

dynRange    = [-60 0];
% xPolar      = 
% zPolar      =
% xPolarACS   =
% zPolarACS   = 
figure('Units','centimeters', 'Position',[5 5 6 6]);


[ax1,~] = imOverlayPolar(BmodeFull,ones(m,n),dynRange,attRange,0, ...
    xPolar,zPolar,xPolarACS,zPolarACS);
yticks(ax1,[4 8 12 16])
title(ax1,'SWTV-ACE')
xlabel(ax1,'Lateral [cm]'), ylabel(ax1,'Axial [cm]')
xlim([-9 9]), ylim([0 18])
hold on
contour(xPolar*1e2, zPolar*1e2, roi,1,'w--')
hold off
% exportgraphics(gcf,fullfile(figsDir,"sample"+samName+"_swtv.png"), ...
%     'Resolution','300')

[TH_acs,R_acs] = meshgrid(-x_ACS*pi/180 + pi/2,z_ACS/100 + r0);
[xPolarACS,zPolarACS] = pol2cart(TH_acs,R_acs);
zPolarACS = zPolarACS + z0Polar;


%% FIGURES
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
axis equal 
axis tight
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

%%
caption = strrep(samName(1:end-4), '_', ' ');
fontSize = 12;

%% First Figure – RFM (e.g., angular view)
figure('Name', 'RFM View', 'Units','normalized','Position',[0.1 0.4 0.275 0.5]);

imagesc(xr, zr*1e3, bMode);
cbar = colorbar;
ylabel(cbar, '[dB]', 'FontSize', 12);
clim([-55 0]);    
colormap gray;
ylabel('Depth [mm]', 'FontSize', 14);
xlabel('\theta [°]', 'FontSize', 14);
% axis ("equal")
title(caption, 'FontSize', fontSize);

%% Second Figure – Polar View
figure('Name', 'Polar View', 'Units','normalized','Position',[0.55 0.4 0.4 0.5]);

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
set(gca, 'Color', 'k');  % Make polar plot background black

%%
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
