% BEAMFORMING (bf_planewave)
% of K-WAVE SIMULATIONS WITH PLANE WAVE AND STERRING ANGLE 
clear all, clc;
%% TEST DATA ONE CASE
% rf_prebf1_ang0p0_a_0p4_1p0
% data = load('Q:\emiranda\proj\dof\out\M05_W4_pwi_v1\rf_prebf1_ang0p0_a_0p4_1p0.mat');

data = load('Q:\emiranda\proj\dof\out\IUS_PWI_v1_0p3\rf_prebf1_ang0p0_a_0p3.mat');

% IUS_PWI_v1_0p3

raw_data    = data.rf_prebf;
fs          = data.fs;
angles_deg  = data.steering_angle;

fnumber     = 2;
pars.sos    = 1540; % [m/s]
pars.pitch  = 0.30e-3; % [m]

[rf, rf_ang] = bf_planewave_comp(raw_data, fs, fnumber, angles_deg, pars);

%%

range_bmode = [-60 0];
units       = 1e2;

% CROP
offset_plot = 110;

axAxis      = 0:size(rf,1)-1; axAxis = axAxis*1/fs*pars.sos/2;
latAxis     = 0:size(rf,2)-1; latAxis = latAxis-mean(latAxis); latAxis = latAxis*pars.pitch;

% FINAL 
rff         = rf(offset_plot:end,:);
z           = axAxis(offset_plot:end);
x           = latAxis;

bmode       = my_RF2Bmode(rff);

figure, 
imagesc(x*units, z*units, bmode, range_bmode)
xlabel('Lateral [cm]'), ylabel('Depth [cm]')
title('Bmode'),
colorbar, colormap("gray")
axis("image")


figure, 
imagesc(bmode, range_bmode)
title('Bmode'),
colorbar, colormap("gray")
% axis("image")

