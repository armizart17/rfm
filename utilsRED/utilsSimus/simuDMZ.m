clearvars;
% Cluster
addpath(genpath('/opt/MATLAB Add-Ons/')); 
addpath(genpath('./utils/'))
addpath(genpath(pwd));

% simulation settings
DATA_CAST           = 'gpuArray-single';
RUN_SIMULATION      = true;         % set to false to reload previous results instead of running simulation

% =========================================================================
% DEFINE THE K-WAVE GRID
% =========================================================================

% Setting the size of the perfectly matched layer (PML)
PML_X_SIZE = 40;    % [grid points]
PML_Y_SIZE = 40;    % [grid points]

% Setting the total number of grid points, not including the PML
Nx = 1580 - 2*PML_X_SIZE;  % [grid points] 9 cm 
Ny = 1580 - 2*PML_Y_SIZE;  % [grid points] 9 cm

% Calculating the space between the grid points
dx = 0.06e-3;       % [m]   
dy = dx;

% Creating the k-space grid
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% Define the radius of the phantom and inclusions
radius = 666;                    % r = 4 cm (666.67)       
radius_inclusion_1 = 166;         % r = 1 cm (166.67)
radius_inclusion_2 = 83;          % r = 0.5 cm (83)

% Create a circular mask for the main phantom
[X, Y] = meshgrid(1:Nx, 1:Ny);
center_x = Nx / 2;
center_y = 5 * Ny / 9;
circular_mask = ((X - center_x).^2 + (Y - center_y).^2) <= radius^2;
circular_mask2 = ((X - center_x).^2 + (Y - center_y).^2) <= (radius-3)^2;

% Create a circular mask for the first inclusion
center_inclusion_x_1 = 5 * Nx / 9;
center_inclusion_y_1 = 5 * Ny / 9;
circular_mask_inclusion_1 = ((X - center_inclusion_x_1).^2 + (Y - center_inclusion_y_1).^2) <= radius_inclusion_1^2;

% Create a circular mask for the second inclusion
center_inclusion_x_2 = 4 * Nx / 11; % Different position for the second inclusion       8       4
center_inclusion_y_2 = 5 * Ny / 9;
circular_mask_inclusion_2 = ((X - center_inclusion_x_2).^2 + (Y - center_inclusion_y_2).^2) <= radius_inclusion_2^2;

% =========================================================================
% DEFINE THE MEDIUM PARAMETERS
% =========================================================================

% Define the properties of the propagation medium
c0 = 1540;              % SoS [m/s] 
c1 = 1500;              % SoS [m/s]
medium.sound_speed = c0 * ones(Nx, Ny);
medium.sound_speed(:,:) = medium.sound_speed(:,:) .* circular_mask + c1 * ~circular_mask;

rho = 1000;                                         % [kg/m^3]

alpha_coeff = 0.5;                                  % [dB/(MHz^y cm)]
alpha_coeff_inclusion_1 = 0.8;
alpha_coeff_inclusion_2 = 0.7;                      % New inclusion alpha value
alpha_coeff_water = 0.0022;

medium.alpha_power  = 1.05;

BonA = 6;
BonA_inclusion_1 = 11; 
BonA_inclusion_2 = 9;                               % New inclusion BonA value
BonA_water = 5.1;

background_map_mean = 1;
background_map_std = 0.02;

% Creating the time array
t_end = (Nx * dx) * 2.2 / c0;     % [s]
kgrid.makeTime(c0, [], t_end);

% =========================================================================
% DEFINE THE INPUT SIGNAL
% =========================================================================

% Define properties of the input signal
tone_burst_freq     = 4e6;        % [Hz] {5}
tone_burst_cycles   = 10;
source_strength_vector = [80 400]*1e3;

for ss = 1:length(source_strength_vector)
    
    source_strength = source_strength_vector(ss);  % [Pa]
    input_signal_norm = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);
    input_signal = (source_strength) * input_signal_norm;

    % =========================================================================
    % DEFINE THE ULTRASOUND TRANSDUCER
    % =========================================================================

    % Creating an instance of kWaveArray and set the dimension to 2D
    dim = '2D';                 
    number_elements = 128;
    element_width   = 8*dx;         
    element_length  = dx;  	        
    element_pitch   = 8*dx;  	   
    radius = inf;

    karray = kWaveArray;
    positions = 0 - (number_elements * element_pitch / 2 - element_pitch / 2) + (0:number_elements-1) * element_pitch;
    
    for ind = 1:number_elements
        angle = 90;
        karray.addRectElement([kgrid.x_vec(160), positions(ind)], element_width, element_length, angle);
    end
    
    source.p_mask = karray.getArrayBinaryMask(kgrid);
    input_signal = repmat(input_signal,number_elements,1);
    source.p = karray.getDistributedSourceSignal(kgrid, input_signal);
    sensor.mask = karray.getArrayBinaryMask(kgrid);

    % =========================================================================
    % DEFINE THE MEDIUM PARAMETERS
    % =========================================================================

    % Initialize density map with water density
    BonA_map = BonA_water * ones(Nx, Ny);
    % Apply the phantom density
    BonA_map(circular_mask) = BonA;
    % Apply the first inclusion density
    BonA_map(circular_mask_inclusion_1) = BonA_inclusion_1;
    % Apply the second inclusion density
    BonA_map(circular_mask_inclusion_2) = BonA_inclusion_2;

    % Initialize alpha coefficient map with water value
    alpha_coeff_map = alpha_coeff_water * ones(Nx, Ny);
    % Apply the phantom alpha coefficient
    alpha_coeff_map(circular_mask) = alpha_coeff;
    % Apply the first inclusion alpha coefficient
    alpha_coeff_map(circular_mask_inclusion_1) = alpha_coeff_inclusion_1;
    % Apply the second inclusion alpha coefficient
    alpha_coeff_map(circular_mask_inclusion_2) = alpha_coeff_inclusion_2;

    % Define the region of interest (ROI) that includes the phantom and its inclusions
    roi_mask = circular_mask | circular_mask_inclusion_1 | circular_mask_inclusion_2;

    % Extract the bounding box around the ROI
    [rows, cols] = find(roi_mask);
    min_row = min(rows);
    max_row = max(rows);
    min_col = min(cols);
    max_col = max(cols);

    % Extract the sub-matrix containing the phantom and its inclusions
    roi_BonA_map = BonA_map(min_row:max_row, min_col:max_col);
    roi_alpha_coeff_map = alpha_coeff_map(min_row:max_row, min_col:max_col);
    
    % =========================================================================
    % RUN THE SIMULATION FOR EACH ROTATION
    % =========================================================================
    
    for num = 1:1

        filename = ['./data/Views_density_map_sam_2D_circ_pha_mi_inc_mi_V10/views_circ_water_dens_u_', num2str(background_map_mean), 'std', strrep(num2str(background_map_std), '.', 'p'), '_', num2str(num), '.mat'];             
        load(filename);

        roi_density_map = density_map(min_row:max_row, min_col:max_col);

        for rot = 0:10:350
            disp(['Computing PLANE WAVE slice num=', num2str(num), ' with rotation=', num2str(rot)]);
            
            % Rotate and update BonA map
            rotated_roi_BonA_map = imrotate(roi_BonA_map, rot, 'bilinear', 'crop');
            BonA_map(min_row:max_row, min_col:max_col) = rotated_roi_BonA_map;
            BonA_map(~circular_mask2) = BonA_water;

            % Rotate and update alpha coefficient map
            rotated_roi_alpha_coeff_map = imrotate(roi_alpha_coeff_map, rot, 'bilinear', 'crop');
            alpha_coeff_map(min_row:max_row, min_col:max_col) = rotated_roi_alpha_coeff_map;
            alpha_coeff_map(~circular_mask2) = alpha_coeff_water;

            % Rotate and update density map
            rotated_roi_density_map = imrotate(roi_density_map, rot, 'bilinear', 'crop');
            density_map(min_row:max_row, min_col:max_col) = rotated_roi_density_map;
            density_map(~circular_mask2) = rho;
            
            % Set the input settings
            input_args = {...
                'PMLInside', false, 'PMLSize', [PML_X_SIZE, PML_Y_SIZE], ...
                'DataCast', DATA_CAST, 'DataRecast', true, 'PlotSim', false};
            
            % Update the medium properties for the simulation
            medium.BonA = BonA_map;
            medium.density = density_map;
            medium.alpha_coeff = alpha_coeff_map;

            % Run the simulation
            sensor_data = kspaceFirstOrder2DC(kgrid, medium, source, sensor, input_args{:});
            combined_sensor_data_num = karray.combineSensorData(kgrid, sensor_data);
            rf(:,:,num) = bf_planewave(combined_sensor_data_num, 1/kgrid.dt, dx, karray.number_elements);

            % Save the results
            file_out = ['V_', num2str(rot), '_circ_water_karray_pw_rf', num2str(source_strength / 1000), 'kPa_sam_map_BonA_inc_hom', num2str(BonA), ...
            'map_alphacoeff_hom', strrep(num2str(alpha_coeff), '.', 'p'), '_f', strrep(num2str(medium.alpha_power), '.', 'p'), ...
            'map_dens_hom_u', num2str(background_map_mean), 'std', strrep(num2str(background_map_std), '.', 'p')];

            pathOut = './data/out/BonA_Att_Dens_karray_water_circ_views10_mi_V11/';
            if ~exist("pathOut", "dir")
                mkdir(pathOut);
            end    
    
            save([pathOut, file_out], 'rf', 'BonA_map', 'alpha_coeff_map', 'density_map');
            clear rf;
        end
    end
end
