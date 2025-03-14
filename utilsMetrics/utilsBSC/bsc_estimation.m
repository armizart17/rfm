function [BSC, band] = bsc_estimation(DATA, REF, pars)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS: 
%         - DATA : Sample data
%                 - DATA.RF
%                 - DATA.acs
%         - REF  : Reference phantom
%                 - REF.RF
%                 - REF.acs
%                 - REF.BSC_ref
%         - pars : 
%                 - pars.z_ini; ROI ini axial (depth)
%                 - pars.z_end; ROI end axial (depth)
%                 - pars.x_ini; ROI ini lateral
%                 - pars.x_end; ROI ini lateral
%                 - pars.fs
%                 - pars.c
%                 - pars.BW
%                 - pars.nb_lambda_axial
%                 - pars.overlap_axial;
%                 - pars.overlap_lateral;
%                 - pars.P (number of points FFT)
% OUTPUTS: 
%         - BSC: Object containing data BSC axial*lateral*band
% AUTHORs: José Timaná - Edmundo Arom Miranda & Andres Coila, based on LIM repository
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% add here the ACS map premade, make an algorithm to make that segmentation
% of layers in the tissue, draw and select an att -> use codes instead of
% values, then change the mask to values when computing, because that
% can change depending on the AC value you want to in it.
% x and z at the starting time

% Transmisson coefficient for Saran layer
T = 1; %transmission_membrane(band);

if isempty(DATA.x)
    DATA.x = 0:size(DATA.RF,2)-1;
    % lateral plot in m
    DATA.x = DATA.x*(pars.BmodeWidth*(10^(-3))/size(DATA.RF,2));
    % Depth resolution
    DATA.z = (0:size(DATA.RF,1)-1)*(pars.c/2)/pars.fs;

    last = min([length(DATA.RF), length(REF.RF)]);
    DATA.z = DATA.z(:,1:last);
end

% The bandwidth:
% lambdac = (pars.c/1e3)/mean(pars.BW); % The wavelength at the transducer center frequency

nb_lambda_axial = pars.nb_lambda_axial; % [wavelength]

overlap_axial = pars.overlap_axial;   
overlap_lateral = pars.overlap_lateral;

%Bande frequentiel de calcul
P = pars.P; %for number of points 
fct = 1e6; %for MHz
Xfreq = (0 : (P/2-1))*pars.fs/P/fct;    % [MHz]                           % [MHz]
lambda = pars.c/(mean([pars.BW(1) pars.BW(2)])*1e6);% [wavelength]
% fprintf('\nCalculation bandwidth:\nBW = %g - %g MHz, Fc = %g MHz\nwavelength (wl) = %g mm\n', pars.BW(1), pars.BW(2), mean([pars.BW(1) pars.BW(2)]), lambda*1e3);
Fi = find(Xfreq >= pars.BW(1));
Ff = find(Xfreq >= pars.BW(2));
band_ind = Fi(1) : Ff(1);               % [sample]
band = band_ind*pars.fs/fct/P;          % [MHz]
% band_k = 2*pi*band*1e6/pars.c;   % [rad]

% Choice of window type
window_type = 2; % 1. Hanning, 2. Tuckey-0.25, 3. Hamming, 4. Chebyshev

Zinitial = pars.z_ini;
Zfinal   = pars.z_end;

Xinitial = pars.x_ini;
Xfinal = pars.x_end;


ind_x = Xinitial <= DATA.x & DATA.x <= Xfinal;
ind_z = Zinitial <= DATA.z & DATA.z <= Zfinal;

Jmin = find(ind_x == 1,1);
Jmax = find(ind_x == 1,1,'last');

Imin = find(ind_z == 1,1);
Imax = find(ind_z == 1,1,'last');

% Jmin = find(DATA.x > Xinitial);  
% Jmax = find(DATA.x > Xfinal);
% 
% Imin = find(DATA.z > Zinitial);  
% Imax = find(DATA.z > Zfinal);

Zaxis = DATA.z(Imin(1):Imax(1));
Xaxis = DATA.x(Jmin(1):Jmax(1));

DATA_ROI_seg = DATA.RF(Imin(1):Imax(1), Jmin(1):Jmax(1), 1);
PHAN_ROI_seg = REF.RF(Imin(1):Imax(1), Jmin(1):Jmax(1), 1);

[I2,J2] = size(DATA_ROI_seg);

% Axial window lenght
axial_gate_length   = nb_lambda_axial*lambda;          % [m]
delta_z = (Zaxis(end)-Zaxis(1))/length(Zaxis);         % [m]
nb_sample_ROI_axial = floor(axial_gate_length/delta_z);% [sample]

% Axial overlap
nb_sample_overlap_axial = fix(nb_sample_ROI_axial*overlap_axial);% [sample]
nb_ROI_axial = fix((I2-nb_sample_overlap_axial)/(nb_sample_ROI_axial - nb_sample_overlap_axial));
% delta_interpixel_axial = (Zaxis(end)-Zaxis(1))/ nb_ROI_axial;
% fprintf('\nAxial gate length = %d wl = %g mm = %g samp.\nNumber of axial region = %g\n', ...
%     nb_lambda_axial, axial_gate_length*1e3, nb_sample_ROI_axial, nb_ROI_axial);

% Lateral windows length
%     nb_lambda_lateral = nb_line_ROI_lateral*delta_lat/lambdac;
nb_lambda_lateral = nb_lambda_axial;
lateral_gate_length  = nb_lambda_lateral*lambda;               % [m]
indice_lateral = find((Xaxis-Xaxis(1)) >= lateral_gate_length);
nb_line_ROI_lateral = indice_lateral(1); clear indice_lateral;

% Recouvrement inter-fenetre lateral
nb_sample_overlap_lateral = fix(nb_line_ROI_lateral*overlap_lateral);    % [sample]
nb_ROI_lateral = fix((J2-nb_sample_overlap_lateral)/(nb_line_ROI_lateral - nb_sample_overlap_lateral));
% delta_interpixel_lateral = (Xaxis(end)-Xaxis(1))/ nb_ROI_lateral;
% fprintf('\nLateral gate length = %2.1f wl = %g mm = %g lines\nNumber of lateral region = %g\n',...
% nb_lambda_lateral, lateral_gate_length*1e3, nb_line_ROI_lateral, nb_ROI_lateral );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X_ROI   = zeros(1, nb_ROI_lateral);
Z_ROI   = zeros(1, nb_ROI_axial);
% Z_ROI_ini = zeros(1, nb_ROI_axial);

SR = zeros(nb_ROI_axial, nb_ROI_lateral, length(band_ind));
BSC   = nan(nb_ROI_axial, nb_ROI_lateral, length(band));

% Main loop -> data block per data block
% where averaged power spectrum are calculated and employed for
% the parameter estimation
for i = 1 : nb_ROI_axial
    pixi = (i-1)*(nb_sample_ROI_axial - nb_sample_overlap_axial) + (1:nb_sample_ROI_axial-1);
    wind = window_choice(length(pixi),window_type);
    if i == nb_ROI_axial
        u = 0;
        while pixi(end) > I2
            u = u + 1; clear pixi mat_window wind;
            pixi = (i-1)*(nb_sample_ROI_axial - nb_sample_overlap_axial) ...
                + (1:(nb_sample_ROI_axial-u));
            warning('\nFor i = %d, ROI final sample (%d) is larger than the RF-signal size (%d)\n', ...
                i, pixi(end), I2);
        end
        wind = window_choice(length(pixi),window_type);
    end
    
    for j = 1 : nb_ROI_lateral
        pixj = (j-1)*(nb_line_ROI_lateral - nb_sample_overlap_lateral) + (1:nb_line_ROI_lateral);
        mat_window = (ones(length(pixj),1)*wind')';
        if j == nb_ROI_lateral
            w = 0;
            while pixj(end) > J2
                w = w + 1; clear pixj;
                pixj = (j-1)*(nb_line_ROI_lateral - nb_sample_overlap_lateral) + (1:nb_line_ROI_lateral-w);
                warning('\nFor j = %d, ROI final sample (%d) is larger than the ROW size (%d)\n', ...
                    j, pixj(end), J2);
            end
            mat_window = (ones(length(pixj),1)*wind')';
        end
        
        Z_ROI(i) = Zaxis( round((pixi(end) + pixi(1))/2 ) )*1e2; % [cm]
        X_ROI(j) = Xaxis( round((pixj(end) + pixj(1))/2 ) )*1e2; % [cm]
        
        % Selection of data in data-block
        phan_block = PHAN_ROI_seg(pixi, pixj,:);
        data_block = DATA_ROI_seg(pixi, pixj);
        
        there_is_nan = double(isnan(data_block));
        [z_find_nan,x_find_nan] = find(there_is_nan == 1);
        
        if  isempty(z_find_nan) == 0 && isempty(x_find_nan) == 0
            
        elseif (isempty(z_find_nan) == 1 && isempty(x_find_nan) == 1)
            
            % Remove the continuous component
            phan_block = remove_DC_in_pixel_ref(squeeze(mean(phan_block)), phan_block, mat_window);
            data_block = remove_DC_in_pixel(mean(data_block), data_block, mat_window);
            
            % Averaged power spectra
            tf_phan = power_spectrum_averaged_in_pixel_ref(phan_block,P,0,0);
            tf_data = power_spectrum_averaged_in_pixel(data_block,P,0,0);

            % Selection of the useful bandwidth
            
            tf_phan = tf_phan(band_ind)./T(:); % If both phantoms have Saran wrap, T = 1
%             tf_phan = tf_phan(band_ind);
            tf_data = tf_data(band_ind);
            
            % Spectral ratio for BSC estimation
            SR(i,j,:) = tf_data ./ tf_phan;
            
            %% Cumulative attenution for BSC compensation for uniform media
            % See Eq. (5) in 10.1016/j.ultras.2021.106376

            % Compensation of the attenuation
            X_ACS_ind = find(1e2*DATA.x >= X_ROI(j));
            Z_ACS_ind = find(1e2*DATA.z >= Z_ROI(i));

%             Z_ROI_ini(i) = Zaxis( pixi(1) ) *1e2; % [cm]
%             L = Zaxis( pixi(end) ) *1e2 - Zaxis( pixi(1) ) *1e2; % [cm]

            % Attenuation of the sample phantom
            ACS_sam = ones(size(DATA.RF)).*DATA.acs/8.686;
            
%             ATT_sam = DATA.acs*band/8.686;
%             cumul_att_data =  exp(-4 * Z_ROI_ini(i) * ATT_sam) .* ( (1 - exp(-ATT_sam * L))./(ATT_sam * L) ).^2; % [Np]  Uniform medium
            cumul_att_data =  exp(-4 * (delta_z*1e2)*sum(ACS_sam(1:Z_ACS_ind(1), X_ACS_ind(1)) ).* band); % [Np/cm]

            % Attenuation of the reference phantom
            ACS_ref = REF.acs/8.686;

%             ATT_ref = REF.acs* band/8.686;
%             cumul_att_phan =  exp(-4 * Z_ROI_ini(i) * ATT_ref) .* ( (1 - exp(-ATT_ref * L))./(ATT_ref * L) ).^2; % [Np]  Uniform medium   
            cumul_att_phan = exp(-4 * Z_ROI(i) * ACS_ref * band);

            % Backscatter coefficient estimation
            if isempty(REF.BSC_ref)
                BSC_data_w_uniform_map = squeeze(squeeze(SR(i,j,:)))'.* cumul_att_phan ./ cumul_att_data;
            else
                BSC_data_w_uniform_map = squeeze(squeeze(SR(i,j,:)))'.* REF.BSC_ref(band) .* cumul_att_phan ./ cumul_att_data;
            end

            BSC(i,j,:)        =  BSC_data_w_uniform_map;
        end
    end
end
end