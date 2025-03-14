function [DATA,X_ROI,Z_ROI,delta_SNR,SNR,est_CF,l_CF,r_CF] = dSNR(Fs, BmodeWidth, x1, x2, y1, y2, c0, BW, tissue, line_ROI_lateral, lambda_axial, dB)
% [DATA,Xaxis,Zaxis,delta_SNR,SNR,est_CF,l_CF,r_CF]
% ---------------------------------------------------------------
% DESCRIPTION
% Provides the b, n and a maps of the sample

% Wavelength at the transducer center frequency
lambda = c0/(mean([BW(1) BW(2)])*1e6);% [wavelength]
% Transducer footprint in mm
transducer_footprint = BmodeWidth;
%%
% Size of each ROI in lines
nb_line_ROI_lateral = line_ROI_lateral;
% Size of each ROI
nb_lambda_axial = lambda_axial;
% Overlapping
overlap_axial = 0.8;
overlap_lateral = 0.8;

% Heterogeneous tissue
DATA.RF    = tissue;

lateral_plot = 0:size(tissue,2)-1;
lateral_plot = lateral_plot*(transducer_footprint*(10^(-3))/size(tissue,2)); % lateral plot in m
DATA.x     = lateral_plot; % Always the same number of lines
DATA.z     = (0:size(DATA.RF,1)-1)*(c0/2)/Fs; % Depth resolution can change

last       = min([length(DATA.z), length(DATA.RF)]);

DATA.z     = DATA.z(:,1:last);
% DATA.RF    = DATA.RF(1:last,:,1); % Not neccesary as we are given frames
DATA.RF    = DATA.RF(1:last,:,:);
env = abs(hilbert(DATA.RF(1:last,:,:)));
%%
P     = 2^10;
fct   = 1e+6; % MHz Factor
Xfreq = (0:(P/2-1))*Fs/P/fct; % [MHz] Getting half of frequency

window_type = 2;

% % Number of lines in each Bmode
latlines  = size(DATA.RF,2);
% % Distance in between adjacent lines
delta_lat = transducer_footprint/latlines;

% Axis to mm, the data is plotted in cm
axis_x = 1000*DATA.x;
axis_z = 1000*DATA.z;

% ROI size in mm
% xs are lateral distance
% ys are axial distance

% counter clockwise rectangle coordinates
seg.location_lesion(1,2) = y1;
seg.location_lesion(4,2) = y1;
seg.location_lesion(2,2) = y2;
seg.location_lesion(3,2) = y2;
seg.location_lesion(1,1) = x1;
seg.location_lesion(4,1) = x2;
seg.location_lesion(2,1) = x1;
seg.location_lesion(3,1) = x2;

% Transforming from distance to pixel values:
seg.position_lesion = seg.location_lesion;
seg.position_lesion(:,1) = 1 + (length(axis_x)-1)/(axis_x(end) - axis_x(1))*(seg.position_lesion(:,1) - axis_x(1)); % Get the index for axis_x
seg.position_lesion(:,2) = 1 + (length(axis_z)-1)/(axis_z(end) - axis_z(1))*(seg.position_lesion(:,2) - axis_z(1)); % Get the index for axis_z
seg.mask_lesion = poly2mask(seg.position_lesion(:,1),seg.position_lesion(:,2),size(DATA.RF,1),size(DATA.RF,2));
seg.mask_lesion = double(seg.mask_lesion);
seg.mask_lesion(seg.mask_lesion == 0) = NaN;
MASK = seg.mask_lesion;
data_seg = DATA.RF.*MASK;
env_seg = env.*MASK;

[z_ROI, x_ROI] = find( MASK == 1 ); % These are the indexes
x_ROI_min = min(x_ROI); x_ROI_max = max(x_ROI); % For some reason the min it is one above the original
z_ROI_min = min(z_ROI); z_ROI_max = max(z_ROI); % For some reason the min it is one above the original
Zaxis = DATA.z(z_ROI_min:z_ROI_max);
Xaxis = DATA.x(x_ROI_min:x_ROI_max);

DATA_ROI_seg = data_seg(z_ROI_min:z_ROI_max, x_ROI_min:x_ROI_max,:); % Why not just use DATA.RF instead of making a previous mask
env_DATA_ROI_seg = env_seg(z_ROI_min:z_ROI_max, x_ROI_min:x_ROI_max,:);


[I2,J2,~] = size(DATA_ROI_seg);
axial_gate_length = nb_lambda_axial*lambda;             % [m]
delta_z = (Zaxis(end)-Zaxis(1))/length(Zaxis);          % [m]
nb_sample_ROI_axial = floor(axial_gate_length/delta_z); % [sample]

% Axial overlap
nb_sample_overlap_axial = fix(nb_sample_ROI_axial*overlap_axial); % [sample]
nb_ROI_axial = ceil((I2-nb_sample_overlap_axial)/(nb_sample_ROI_axial - nb_sample_overlap_axial)); % 26*(417-333)+333=>2493 == (26-1)*(417-333)+417=>2493
fprintf('\nAxial gate length = %d wl = %g mm = %g samp.\nNumber of axial region = %g\n',nb_lambda_axial,axial_gate_length*1e+3,nb_sample_ROI_axial,nb_ROI_axial);

% Lateral windows length
nb_lambda_lateral = nb_line_ROI_lateral*delta_lat/(lambda*10^(3)); % divide by 10^3 to pass to m as delta_lat is in mm
lateral_gate_length   = nb_lambda_lateral*lambda; % [m]

% Recouvrement inter-fenetre lateral
nb_sample_overlap_lateral = fix(nb_line_ROI_lateral*overlap_lateral); % [sample]
nb_ROI_lateral = fix((J2-nb_sample_overlap_lateral)/(nb_line_ROI_lateral - nb_sample_overlap_lateral)); % 164*(10-8)+8=>336 == (164-1)*(10-8)+10=>336
fprintf('\nLateral gate length = %2.1f wl = %g mm = %g lines\nNumber of lateral region = %g\n\n',nb_lambda_lateral, lateral_gate_length*1e3, nb_line_ROI_lateral, nb_ROI_lateral);

% Matrix initialization
X_ROI    = zeros(1,nb_ROI_lateral);
Z_ROI    = zeros(1,nb_ROI_axial);
%% Get SR ratio
delta_SNR = zeros(nb_ROI_axial, nb_ROI_lateral);
SNR = zeros(nb_ROI_axial, nb_ROI_lateral);
est_CF = zeros(nb_ROI_axial, nb_ROI_lateral);
l_CF = zeros(nb_ROI_axial, nb_ROI_lateral);
r_CF = zeros(nb_ROI_axial, nb_ROI_lateral);
% z = zeros(nb_ROI_axial, nb_ROI_lateral);

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
                % Window for each ROI
                mat_window = (ones(length(pixj),1)*wind')';
                if j == nb_ROI_lateral
                    w = 0;
                    while pixj(end) > J2
                        w = w + 1;
                        clear pixj;
                        pixj = (j-1)*(nb_line_ROI_lateral - nb_sample_overlap_lateral) + (1:nb_line_ROI_lateral-w);
                        warning('\nFor j = %d, ROI final sample (%d) is larger than the ROW size (%d)\n',j,pixj(end),J2);
                    end
                    mat_window = (ones(length(pixj),1)*wind')';
                end
                Z_ROI(i) = Zaxis( round((pixi(end) + pixi(1))/2 ) ); % [m]
%                 z(i,j)=Z_ROI(i);
                X_ROI(j) = Xaxis( round((pixj(end) + pixj(1))/2 ) ); % [m]
                  
                % Selection of data in data-block
                data_block = DATA_ROI_seg(pixi, pixj,:);
                env_block = env_DATA_ROI_seg(pixi, pixj,:);
                
                there_is_nan = double(isnan(data_block));
                [z_find_nan,x_find_nan] = find(there_is_nan == 1);
                  
                if isempty(z_find_nan) == 0 && isempty(x_find_nan) == 0
                  
                elseif (isempty(z_find_nan) == 1 && isempty(x_find_nan) == 1)
                    % dSNR
                    SNR(i,j) = mean(env_block(:))./std(env_block(:)); %(mean(env_block(:)))./sqrt(mean((env_block(:)- mean(env_block(:))).^2));
                    delta_SNR(i,j) = abs(SNR(i,j) - 1.91).*100./1.91;
                    
                    % Remove the continuous component
                    data_block = remove_DC_in_pixel_ref(squeeze(mean(data_block)), data_block, mat_window);
                    % Averaged power spectra
                    tf_data = power_spectrum_averaged_in_pixel_ref(data_block,P,0,0);

                    % CFDS
                    tf_data_dB = 10*log10(tf_data)-max(10*log10(tf_data));

%                     [m_int,m_int_index] = max(tf_data);
                    [~,m_int_index] = max(tf_data_dB);
%                     left_freq_index = find(tf_data(1:m_int_index) <= m_int/2,1,'last');
%                     right_freq_index = find(tf_data(m_int_index:end) <= m_int/2,1,'first');


                    left_freq_index = find(tf_data_dB(1:m_int_index) <= dB,1,'last');
                    right_freq_index = find(tf_data_dB(m_int_index:end) <= dB,1,'first');
                    right_freq_index = right_freq_index + m_int_index;
                  
                    if isempty(left_freq_index) || isempty(right_freq_index)
                        % Actually right_freq would never be empty, but just in case we add it
                        % left_freq could be empty when the lobe crosses the x axis
                        est_CF(i,j) = Xfreq(m_int_index);
                    else
                        est_CF(i,j) = Xfreq(round(left_freq_index+(right_freq_index-left_freq_index)/2));
                    end

                    if isempty(left_freq_index)
                        l_CF(i,j) = nan;
                    else
                        l_CF(i,j) = Xfreq(left_freq_index);
                    end

                    if isempty(left_freq_index)
                        r_CF(i,j) = nan;
                    else
                        r_CF(i,j) = Xfreq(right_freq_index);
                    end
                end
           end
end