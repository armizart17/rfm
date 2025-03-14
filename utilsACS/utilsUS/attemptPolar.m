
% param = getparam('C5-2v');
siz = size(rf);
z = SAM.z;
zmax = z(end);
Rad = param.radius;
pit = param.pitch;
Nel = param.Nelements;

% SIEMENS 6C1 HD
Nel = 192;
L = 2*R*sin(asin(pit/2/Rad)*(Nel-1)); % chord length
d = sqrt(R^2-L^2/4); % apothem
z0 = -d;

th = -(linspace(atan2(L/2,d),atan2(-L/2,d),siz(2)))*180/pi;
r = linspace(R+p,-z0+zmax,siz(1));
r0 = r(1);

% To Polar Coordinates
[xPolar,zPolar, z0Polar] = impolgrid(size(b_mode), z(end),param);
%%

units           = 1E2;
bmodeFull       = bmode_sam;
colorImg        = alpha_ratio;
range_bmode     = [-60 0];
range_img       = [0 0.7];
transparency    = 0.65;
x_img           = spectralData_sam.lateral*units;
z_img           = spectralData_sam.depth*units;
xFull           = SAM.x*units;
zFull           = SAM.z*units;
[X, Z] = meshgrid(xFull, zFull);
roi = and(X >= x_img(1), X <= x_img(end)) & ...
      and(Z >= z_img(1), Z <= z_img(end));


%% POLAR

% 
units           = 1E0;
x_img           = spectralData_sam.lateral*units;
z_img           = spectralData_sam.depth*units;
bmodeFull       = bmode_sam;

% SIEMENS 6C1 HD
param.radius        = 49.91E-3;
param.Nelements     = 192; 
param.pitch         = 0.318E-3;
[xPolar,zPolar, z0Polar] = impolgrid(size(bmodeFull), SAM.z(end), param);

L = 2*param.radius*sin(asin(param.pitch /2/param.radius)*(param.Nelements-1)); % chord length
d = sqrt(param.radius^2 - L^2/4); % apothem
z0 = -d;
r = linspace(param.radius + param.pitch, -z0+SAM.z(end), size(SAM.rf,1) );
r0 = r(1);



[TH_acs,R_acs] = meshgrid(-x_img*pi/180 + pi/2, z_img + r0);
[xPolarACS, zPolarACS] = pol2cart(TH_acs, R_acs);
zPolarACS = zPolarACS + z0Polar;

figure('Units','centimeters', 'Position',[5 5 6 6]);
% [ax1,~] = imOverlayPolar(bmodeFull, ones(p,q) ,range_bmode, [8 11], 0.5, ...
[ax1,~] = imOverlayPolar(bmodeFull, b_ratio_dB, range_bmode, [8 11], 0.5, ...
    xPolar,zPolar,xPolarACS,zPolarACS);
yticks(ax1,[4 8 12 16])
title(ax1,'B-mode')
xlabel(ax1,'Lateral [cm]'), ylabel(ax1,'Axial [cm]')
xlim([-8 8]), ylim([0 18])
hold on
contour(xPolar*1e2, zPolar*1e2, roi,1,'w--')
hold off
%%







figure('Units','centimeters', 'Position',[5 5 6 6]);
[ax1,~] = imOverlayPolar(bmode_sam,ones(p,q),range_bmode,colorImg,0, ...
    xPolar,zPolar,xPolarACS,zPolarACS);
yticks(ax1,[4 8 12 16])
title(ax1,'B-mode')
xlabel(ax1,'Lateral [cm]'), ylabel(ax1,'Axial [cm]')
xlim([-8 8]), ylim([0 18])
hold on
contour(xPolar*1e2, zPolar*1e2, roi,1,'w--')
hold off