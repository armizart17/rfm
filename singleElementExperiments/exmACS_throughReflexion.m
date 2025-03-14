%% EXAMPLE FOR THROUGH REFLEXION METHOD FOR ACS ESTIMATION
% EMZ Nov2024

% General specs
NFFT            =  ;
Fs              =  ; % [Hz]
dist            =  ; % width of the phantom [cm]


%% Signal plexy
signal_plexy    =  ; % [V]
time_vec        =  (0:length(signal_plexy)-1)/Fs; %[s])

Y_plexy = fftshift(fft(signal_plexy));

S_plexy = abs(Y_plexy) .^ 2; % Power Spectrum

S_plexy_dB = 10*log10(S_plexy);

axis_f = (0:(NFFT-1))'*fs/NFFT;  % [Hz] 
%% Transmission loss
rho     = 1.05;     % densisty [g/cm^3]
cs      = 1470.5;   % speed of sound [m/s]

rho_ref = 1;        % densisty [g/cm^3] reference water
cs_ref  = 1480.81;  % speed of sound [m/s] reference water

D1      = 2*(rho_ref*cs_ref)/((rho_ref*cr)+(rho*c)); % Tcw
D2      = 2*(rho*cs)/((rho_ref*cs_ref)+(rho*c));     % Twc
TL      = (D1*D2)^2; % Transmission_loss


%% Signal with phantom and plexy
signal_phantom  =  ; % [V]
time_vec        =  (0:length(signal_phantom)-1)/Fs; %[s]

Y_phantom = fftshift(fft(signal_phantom));

S_phantom = abs(Y_signal_vec) .^ 2; % Power Spectrum

S_phantom_dB = 10*log10(S_phantom /TL);

axis_f = (0:(NFFT-1))'*fs/NFFT;  % [Hz] 

%% Display Norm Power in dB
S_plexy_normdB   = S_plexy(1:end_pos) / max(S_plexy(1:end_pos));
S_phantom_normdB = S_phantom(1:end_pos) / max(S_phantom(1:end_pos));


figure, 
title('Power Norm (dB)')
plot(freq_pos, S_plexy_normdB, 'DisplayName', 'S_{plex}')
hold on;
plot(freq_pos, S_phantom_normdB, 'DisplayName', 'S_{sam}')
hold off;
xlabel('Freq [MHz]'), ylabel('Norm Power [dB]')


%% Attenuation calculation
att_full = (S_plexy_dB - S_phantom_dB)/(2*dist); 

% Display
end_pos = ceil(NFFT/2)+1;
freq_pos = axis_f(1:end_pos)*1E-6;

figure, 
plot(freq_pos, att_full(1:end_pos) )
xlabel('Freq [MHz]'), ylabel('Power Ratio [dB]')


% Select Bandwith considering -15dB decay 
f_ini   = ; % [MHz]
f_end   = ; % [MHz]

ind_f   = axis_f >= f_ini & axis_f <= f_end;

band_range = axis_f(ind_f)*1E-6; % [MHz]
att_range  = att_full(ind_f); % [dB/cm]

%% FIT 
band_range = band_range(:); % to col vector
att_range  = att_range(:); % to col vector
   
fit_func   = fittype('alpha0 * f^gamma', 'independent', 'f', 'coefficients', {'alpha0', 'gamma'});
    
[fit_result, gof] = fit(band_range, att_range, fit_func);
    
alpha0  = fit_result.alpha0;
gamma   = fit_result.gamma;
        
slope   = alpha0 * gamma*f.^(gamma-1); % derivate

% Display  
figure,
plot(band_range, att_range, '--');
hold on;
plot(band_range, fit_result(band_range), 'r-');
hold off;
xlabel('Frequency [MHz]');
ylabel('Attenuation [dB/cm]');
title('Attenuation Coefficient')
legend('Data', 'Fit');