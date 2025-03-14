function [Est_n_sample, Est_b_sample, Est_alphaz_sample,delta_SNR,SNR,est_CF,DATA,Xaxis,Zaxis,SR,z,...
    nb_ROI_axial,nb_ROI_lateral,mu_rpl_tv,mu_rpl_robust_tv] = qUS_estimation(Fs, BmodeWidth, phan_alpha, bsc_band,...
    x1, x2, y1, y2, c0, BW, tissue, phantom, line_ROI_lateral, lambda_axial, estim_method, delta_flag, DATA, delta_n, struct)
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

% Reference homogeneous phantom in Hz
PHAN.fs      = Fs;

% lateral_plot = 0:size(phantom,2)-1;
% % simulation
% 
% lateral_plot = lateral_plot-mean(lateral_plot); % I want it to start from 0
% % lateral_plot = lateral_plot*(transducer_footprint*(10^(-3))/size(phantom,2)); % lateral plot in m
% depth_plot   = (0:size(phantom,1)-1)*(c0/2)/PHAN.fs;

% PHAN.x       = lateral_plot;
% PHAN.z       = depth_plot;
% 
PHAN.RFline  = phantom;
% 
% 
% % Heterogeneous tissue
% DATA.RF    = tissue;
% 
% DATA.x     = PHAN.x; % Always the same number of lines
% DATA.z     = (0:size(DATA.RF,1)-1)*(c0/2)/PHAN.fs; % Depth resolution can change
% 
% last       = min([length(DATA.z), length(PHAN.z), length(DATA.RF), length(PHAN.RFline)]);
% 
% DATA.z     = DATA.z(:,1:last);
% % DATA.RF    = DATA.RF(1:last,:,1); % Not neccesary as we are given frames
% DATA.RF    = DATA.RF(1:last,:,:);
% 
% % Envelope

% env = abs(hilbert(DATA.RF(1:last,:,:)));
env = abs(hilbert(DATA.RF(:,:,:)));
% 
% PHAN.z       = PHAN.z(:,1:last);
% PHAN.RFline  = PHAN.RFline(1:last,:,:);
%%
% Number of lines in each Bmode
latlines  = size(DATA.RF,2);
% Distance in between adjacent lines
delta_lat = transducer_footprint/latlines;

P     = 2^10;
fct   = 1e+6; % MHz Factor
Xfreq = (0:(P/2-1))*PHAN.fs/P/fct; % [MHz] Getting half of frequency

fprintf('\nCalculation bandwidth:\nBW = %g - %g MHz, Fc = %g MHz\nwavelength (wl) = %g mm\n', BW(1), BW(2), mean([BW(1) BW(2)]), lambda*1e3);

Fi = find(Xfreq >= BW(1));
Ff = find(Xfreq >= BW(2));
band_ind = Fi(1) : Ff(1);      % [sample]
band = band_ind*PHAN.fs/P/fct; % [MHz] Bandwidth
window_type = 2;

% Axis to mm, the data is plotted in mm
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
%seg.mask_m_lesion = repmat(seg.mask_lesion, [1 1]); % Not used
%seg.mask = seg.mask_lesion; % Not used
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
PHAN_ROI_seg = PHAN.RFline(z_ROI_min:z_ROI_max, x_ROI_min:x_ROI_max,:);
 
[I2,J2,~] = size(DATA_ROI_seg);
axial_gate_length = nb_lambda_axial*lambda;             % [m]
%Lo2 = axial_gate_length*1e+2/2;                         % [cm]
delta_z = (Zaxis(end)-Zaxis(1))/length(Zaxis);          % [m]
nb_sample_ROI_axial = floor(axial_gate_length/delta_z); % [sample]

% Axial overlap
nb_sample_overlap_axial = fix(nb_sample_ROI_axial*overlap_axial);% [sample]
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
SR = zeros(nb_ROI_axial, nb_ROI_lateral, length(band_ind));
delta_SNR = zeros(nb_ROI_axial, nb_ROI_lateral);
SNR = zeros(nb_ROI_axial, nb_ROI_lateral);
est_CF = zeros(nb_ROI_axial, nb_ROI_lateral);
z = zeros(nb_ROI_axial, nb_ROI_lateral);

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
                
                Z_ROI(i) = Zaxis( round((pixi(end) + pixi(1))/2 ) )*1e2; % [cm]
                z(i,j)=Z_ROI(i);
                X_ROI(j) = Xaxis( round((pixj(end) + pixj(1))/2 ) )*1e2; % [cm]
                  
                % Selection of data in data-block
                phan_block = PHAN_ROI_seg(pixi, pixj,:);
                data_block = DATA_ROI_seg(pixi, pixj,:);
                env_block = env_DATA_ROI_seg(pixi, pixj,:);
                
                there_is_nan = double(isnan(data_block));
                [z_find_nan,x_find_nan] = find(there_is_nan == 1);
                  
                if isempty(z_find_nan) == 0 && isempty(x_find_nan) == 0
                  
                elseif (isempty(z_find_nan) == 1 && isempty(x_find_nan) == 1)
                    % SNR
                    SNR(i,j) = mean(env_block(:))./std(env_block(:));
                    delta_SNR(i,j) = abs(SNR(i,j) - 1.91).*100./1.91;
                    
                    % Remove the continuous component
                    phan_block = remove_DC_in_pixel_ref(squeeze(mean(phan_block)), phan_block, mat_window);
                    data_block = remove_DC_in_pixel_ref(squeeze(mean(data_block)), data_block, mat_window);

                    % Averaged power spectra
                    tf_phan = power_spectrum_averaged_in_pixel_ref(phan_block,P,0,0);
                    tf_data = power_spectrum_averaged_in_pixel_ref(data_block,P,0,0);
                    tf_data_dB = 20*log10(tf_data)-max(20*log10(tf_data));

                    % CFDS
%                     [m_int,m_int_index] = max(tf_data);
                    [~,m_int_index] = max(tf_data_dB);
%                     left_freq_index = find(tf_data(1:m_int_index) <= m_int/2,1,'last');
%                     right_freq_index = find(tf_data(m_int_index:end) <= m_int/2,1,'first');

                    left_freq_index = find(tf_data_dB(1:m_int_index) <= -20,1,'last');
                    right_freq_index = find(tf_data_dB(m_int_index:end) <= -20,1,'first');
                    right_freq_index = right_freq_index + m_int_index;
                    
                    if isempty(left_freq_index) || isempty(right_freq_index)
                        % Actually right_freq would never be empty, but just in case we add it
                        % left_freq could be empty when the lobe crosses the x axis
                        est_CF(i,j) = Xfreq(m_int_index);
                    else
                        est_CF(i,j) = Xfreq(round(left_freq_index+(right_freq_index-left_freq_index)/2));
                    end
                    
                    % Selection of the useful bandwidth
                    tf_phan = tf_phan(band_ind)';
                    tf_data = tf_data(band_ind);
                    
                    % Spectral ratio for BSC estimation

                    if delta_flag == 0
                        SR(i,j,:) = (tf_data./tf_phan').*transmission_membrane(band').*bsc_band;
                    else
                        SR(i,j,:) = (tf_data./tf_phan'); %.*transmission_membrane(band');
                    end
                end
           end
end


mu_rpl_tv =  nan;
mu_rpl_robust_tv = nan;

mu_b  = struct.mu(1);
mu_n  = struct.mu(2);
mu_a  = struct.mu(3);

if estim_method == 1
    %% Original RPL-TV
    % Indices initialization
    f = band;
    p = i;
    q = j;
    r = length(band);
    
    % Derivate matrix
    dy = diag(ones(p-1,1),1) - diag([ones(p-1,1);0]);
    Dy = sparse(kron(speye(q),dy)); %diag(ones(p*q -1,1),1) - diag(ones(p*q,1));
    
    % Depth component derivative
    dz = reshape(Dy*z(:),p,q);
    dz(end,:) = dz(end-1,:);
    
    % Spectrum Ratio
    SR = log(permute(SR,[3,1,2]));

    %% RPL-based algorithm variants
    Y = SR;
    
    % Matrices for RPL-based algorithms
    X = kron(speye(p*q),repmat(ones(r,1),[1,1]));
    Z = kron(speye(p*q),repmat(log(f'),[1,1]));
    W = -4*kron(speye(p*q),repmat(f',[1,1]));
    
    % Implementation parameters
    par_rpl.tol   = 1e-16;
    par_rpl.kmax  = 100;
    par_rpl.eps_f = 1e-16;
    par_rpl.m_est = 0; %Robust
    
    % Parameters for RPL-TV
%     mu_b  = 1e+0;%1e+0; %1e+0; % ORIGINAL 
%     mu_n  = 1e+3;%1e+3; %1e+2; % ORIGINAL 
%     mu_a  = 1e+3;%1e+3; %1e+3; % ORIGINAL 
    mu_rpl_tv    = [mu_b;mu_n;mu_a];
    
    x = rpl_tv_ini(Y,X,Z,W,mu_rpl_tv,par_rpl); % rpl_robust_tv(Y,X,Z,W,mu,c_rtv,sigma,u_0,par_rpl);
    
    
    % % Parameters for RPL-Robust-TV
    % mu_b_rb  = 1e+0; %1e+0;
    % mu_n_rb  = 1e+2; %1e+2;
    % mu_a_rb  = 1e+1; %1e+1;
    % mu_rpl_robust_tv    = [mu_b_rb;mu_n_rb;mu_a_rb];
    mu_rpl_robust_tv = [nan;nan;nan];
    % c_rtv = [0.1; 1; 1e308; 1e308];
    % sigma = [1.4484; 1.4484; 1.4484; 1.4484];
    % 
    % % for k = 1:n_it
    % %     % Regularized Power Law by Total Variation
    % %     %tic;
    % %     x(:,k) = rpl_tv_ini(Y,X,Z,W,mu,par_rpl); %rpl_tv(Y,X,Z,W,mu,par_rpl); rpl_tv_ini(Y,X,Z,W,mu,par_rpl);
    % %     
    % %     %t(k,1) = toc;
    % % end
    % x = rpl_robust_tv(Y,X,Z,W,mu_rpl_robust_tv,c_rtv,sigma,x,par_rpl);
    
    b = x(1:p*q); % x(1:p*q,k);
    n = x(1+p*q:2*p*q); % x(1+p*q:2*p*q,k)
    a = x(1+2*p*q:end); % x(1+2*p*q:end,k)
    
    b_ratio = reshape(exp(b),p,q);
    n_ratio = reshape(n,p,q);
    
    alpha_est = reshape(Dy*a,p,q)./dz;
    alpha_ratio = 8.686*alpha_est;

elseif estim_method == 2
    %% n prior

    SR = permute(SR,[3,1,2]);
    
    [~,p,q] = size(SR);
    
%     delta_n = -1.5; %2.8; %-6; % 2.5 4 nsam-nref
    % bsc_band -> n = delta_n del tejido
    % bsc_band -> n = n + delta_n del tejido

    comp_ref = comp_ref_n_bsc(delta_n,band,p,q);
    
    SR_comp = SR.*comp_ref;
    
    % indices initialization
    f = band;
    [r,p,q] = size(SR_comp);
    
    % log-spectrum Ratio
    Y = log(SR_comp);
    
    % matrices for RPL-based algorithms
    X = kron(speye(p*q),ones(r,1));
    W = -4*kron(speye(p*q),f');
    
    % Implementation parameters
    par_rpl.tol   = 1e-16;
    par_rpl.kmax  = 100;
    par_rpl.eps_f = 1e-16;
    par_rpl.m_est = 0; %Robust
    
    % Parameters for RPL-TV
%     mu_b  = 1e+0;%1e+0; %1e+0; % ORIGINALS
%     mu_n  = 1e+3;%1e+3; %1e+2; % ORIGINALS
%     mu_a  = 1e+3;%1e+3; %1e+3; % ORIGINALS

    mu_rpl_tv    = [mu_b;mu_n;mu_a];
    par_rpl.ini_method = 1;
    par_rpl.ini_tol = 1e-16;
    par_rpl.df_op = 1;
    
    % initialization for RPL-based methods
    u_0 = initialize_rpl_n_prior(Y,X,W,mu_rpl_tv,par_rpl);
    
    % RPL estimation
    [u_opt,~] = rpl_tv_n_prior(Y,X,W,mu_rpl_tv,u_0,par_rpl);
    
    dy = 0.5*(diag(ones(p-1,1),1) - diag(ones(p-1,1),-1));
    dy(1,1) = -1; dy(1,2) = 1; dy(end,end) = 1; dy(end,end-1) = -1;
    Dy = sparse(kron(speye(q),dy));
    
    b = u_opt(1:p*q);
    a = u_opt(p*q+1:2*p*q);
    
    % z = 1e+2*repmat(depth,1,q);
    dz = reshape(Dy*z(:),p,q);
    dz(end,:) = dz(end-1,:);
    
    b_ratio = reshape(exp(b),p,q);
    alpha_ratio  = reshape(8.686*Dy*a./dz(:),p,q);
    n_ratio = delta_n*ones(size(b_ratio));
end

Est_b_sample = b_ratio;
Est_n_sample = n_ratio;

if delta_flag == 0
    % Reference phantom values
    % b_ref     = phan_b; % Backscatter power law fit
    % n_ref     = phan_n; % Backscatter power law fit

    %Est_bsc_sample          = Est_b_sample.*(power(15,Est_n_sample)); % For f = 15e6 [Hz]
    Est_alphaz_sample       = alpha_ratio + phan_alpha;
else
    Est_alphaz_sample       = alpha_ratio;
end


end