% Test generador de pulsos 

% 52Vpp 30dB
% 18Vpp 20dB
% 5.6Vpp  10dB
% 3.6Vpp 5dB

% Therefore is 10^(delta_dB/20)

sos = 1.46;  % mm/us % Cluster.Vel_sound
t_v = 1E-6;
d_v = (t_v * 1e6 * sos) / 2  % Distance vector in mm





