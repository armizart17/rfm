
alpha_value     = 0.5;
numSim          = 1;
nameFolderOut   = 'test';
typeSimu        = 'homo';



%% simuMultiPurposePWI_v1(alpha_value, numSim, nameFolderOut, typeSimu)
% function simuMultiPurposePWI_v1(alpha_value, numSim, nameFolderOut, typeSimu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Classical homogeneous reference approx size widht (3.8cm), height (5.47cm)
% PLANE WAVE IMAGING v1.0 (use for 1 inclusion and 1 reference)
% INPUT
%       alpha_value     : is for alpha_attenuation, when inclusion [alpha_bg alpha_inc]
%       numSim          : number of simulations
%       nameFolderOut   : name of folder to save (It will save data in ./out/nameFolderOut/
%       typeSimu        : type of simu (i.e. 'homo' and 'circle' for now, 
%                           TBD:'layers', check patterns function)
% OUTPUT              
% SAVED DATA in ./out/nameFolderOut/rf(iSim)_a_(alpha_bg)_(alpha_inc).mat'
% EXAMPLE: simuMultiPurposePWI_v1([1.0 0.5], 1, 'PWIv1', 'circle')
% By EMZ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FOR CLUSTER ("simuMultiPurposePWI_v1.sh" file)
% #!/bin/bash
% #SBATCH --output="slurm-%j.out"
% ##SBATCH --time=22:00:00  # 22 hours
% #SBATCH --partition=debug
% #SBATCH --nodes=1  # Ensure only one node is used
% #SBATCH --ntasks=1
% #SBATCH --cpus-per-task=10  # w1,2,7 (20), w8,9,10(16)
% #SBATCH --gres=gpu:1  # Request 1 GPU (modify if needed)
% #SBATCH --mem=MaxMemPerNode  # Use all available memory
% #SBATCH --exclusive  # Reserve the entire node
% #SBATCH --exclude=worker2  # Exclude specific workers
% 
% # Get the number of CPUs dynamically
% TOTAL_CPUS=$SLURM_CPUS_PER_TASK  # Get the number of allocated CPUs
% 
% # Capture all command-line arguments
% args=("$@")
% 
% # Convert to MATLAB-compatible syntax (quote string arguments)
% matlab_args=""
% for arg in "${args[@]}"; do
%     if [[ "$arg" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
%         # If the argument is numeric, use it as is
%         matlab_args+="$arg,"
%     else
%         # If the argument is a string, wrap it in single quotes
%         matlab_args+="'$arg',"
%     fi
% done
% 
% # Remove the last comma
% matlab_args=${matlab_args%,}
% 
% # Run MATLAB with the processed arguments
% srun --ntasks=1 --nodes=1 matlab -nosplash -nodesktop -nodisplay -r "c=parcluster('local'); numWorkers=min(c.NumWorkers, feature('numcores')); maxNumCompThreads(numWorkers); delete(gcp('nocreate')); parpool('local', numWorkers); simuMultiPurposePWI_v1($matlab_args); exit;"

%%
[~, pcname] = system('hostname');

if strcmp(pcname(1:end-1), 'C084285') % PC LIM
    pathNow = 'D:\emirandaz\qus\qus-lim';
else % Cluster
    cd('../'); % when you run it usually goes inside ./codes/script.m
    addpath(genpath(pwd))
    addpath(genpath('/opt/MATLAB Add-Ons'))
    fprintf('PC CLUSTER\n');
end


% list_rho_mean = [1000 1000 1000 1000 1000];
% list_rho_std = [0.01 0.02 0.03 0.04 0.05];

% PRE-SET THE DENSITY MAP
rho_mean = 1000;
rho_std  = 0.02;
Nx = 944; Ny = 655; Nz = 140;
density_map = rho_mean .* (1 + rho_std * randn(Nx, Ny, Nz));

OutputDir = fullfile(pwd, 'out', nameFolderOut);
if ~exist(OutputDir,"dir"); mkdir(OutputDir); disp('OutputDir created'); end

% delete(gcp)
% parpool

% c = parcluster('local');  
% numWorkers = min(c.NumWorkers, feature('numcores'));  % Use the max allowed
% fprintf('Num of workers : %d\n', numWorkers);
% delete(gcp('nocreate'));
% parpool('local', numWorkers);


for iSim = 1:numSim

clearvars -except iSim alpha_value numSim nameFolderOut list_rho_mean list_rho_std OutputDir typeSimu density_map
% clc; 
rng shuffle;

% DATA_CAST = 'single';     % set to 'single' or 'gpuArray-single' to speed up computations
% DATA_CAST = 'gpuArray-single';     % set to 'single' or 'gpuArray-single' to speed up computations
if gpuDeviceCount > 0
    DATA_CAST = 'gpuArray-single';
else
    DATA_CAST = 'single';
end

RUN_SIMULATION  = true;         % set to false to reload previous results instead of running simulation

%% Create the computational grid
elem_pitch = 0.30e-3;

scale_fs = 2;
% pml_x_size = scale_fs*40;                % [grid points]
% pml_y_size = scale_fs*14;                % [grid points]
% pml_z_size = scale_fs*14;                % [grid points]

pml_x_size = 40;                % [grid points]
pml_y_size = 10;                % [grid points]
pml_z_size = 10;                % [grid points]

% set total number of grid points not including the PML
% Nx = scale_fs*810 - 2 * pml_x_size;      % [grid points]
% Ny = scale_fs*540 - 2 * pml_y_size;      % [grid points] ;% /4 for going from 128 to 41 aperture
% Nz = scale_fs*100 - 2 * pml_z_size;      % [grid points]

Nx = 1024 - 2 * pml_x_size;      % [grid points]
Ny = 675 - 2 * pml_y_size;      % [grid points]
Nz = 160 - 2 * pml_z_size;      % [grid points]

% element_width = 4*scale_fs;             % 4 for at least 0.075mm
% ratio = element_width;
% dy = elem_pitch/element_width*2;          % grid point spacing in the y direction [m]
% dx = dy;
% dz = dx;

dx = 0.06e-3;
dy = dx;                        % [m]
dz = dx;                        % [m]

% kgrid = kWaveGrid(Nx, dx, Ny, dy);
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

%% Signals settings

c0 = 1540;                      % [m/s]
rho0 = 1000;                    % [kg/m^3]
t_end   = (Nx*dx)*2/c0;     % [s]
kgrid.makeTime(c0, [], t_end);
fs      = 1/kgrid.dt;

source_strength = 1e6;
tone_burst_freq = 6.66e6;        % [Hz]
tone_burst_cycles = 3.5;

angleMax = 15;       % maximum steering angle in degrees
nAngles  = 3;       % total number of steering angles

steering_angle_vector = linspace(-angleMax, angleMax, nAngles)';

for gg = 1:length(steering_angle_vector)
    clear transducer
    steering_angle = steering_angle_vector(gg);

    input_signal_norm = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);
    input_signal = (source_strength ./ (c0 * rho0)) .* input_signal_norm;
    
    %% Transducer settings
    % physical properties of the transducer
    %transducer.number_elements = 32;  	% total number of transducer elements
    transducer.number_elements  = 128;  	% total number of transducer elements
    transducer.element_width    = 5;       % width of each element [grid points]
    transducer.element_length   = 100;  	% length of each element [grid points]
    transducer.element_spacing  = 0;  	% spacing (kerf  width) between the elements [grid points]
    transducer.radius           = inf;            % radius of curvature of the transducer [m]
    
    % calculate the width of the transducer in grid points
    transducer_width = transducer.number_elements * transducer.element_width ...
        + (transducer.number_elements - 1) * transducer.element_spacing;
    
    % use this to position the transducer in the middle of the computational grid
    % transducer.position = round([1, Ny/2 - transducer_width/2]);
    transducer.position = round([1, Ny/2 - transducer_width/2, Nz/2 - transducer.element_length/2]);
        
    % properties used to derive the beamforming delays
    transducer.sound_speed      = c0;                    % sound speed [m/s]
    transducer.focus_distance   = Inf;              % focus distance [m]
    %transducer.elevation_focus_distance = 19e-3;    % focus distance in the elevation plane [m]
    transducer.elevation_focus_distance = Inf;    % focus distance in the elevation plane [m]
    transducer.steering_angle   = steering_angle;                  % steering angle [degrees]
    
    % apodization
    transducer.transmit_apodization = 'Rectangular';
    transducer.receive_apodization  = 'Rectangular';
    
    % define the transducer elements that are currently active
    %number_active_elements = 32;
    transducer.active_elements = ones(transducer.number_elements, 1);

    transducer.input_signal = input_signal;
    % create the transducer using the defined settings
    transducer = kWaveTransducer(kgrid, transducer);
        
    transducer.properties;

    % =======================================
    % MEDIUM PROPERTIES
    % =======================================

    % SoS
    % medium.sound_speed = patterns(sos_mean, sos_std, 'homo', [], Nx, Ny);
    medium.sound_speed      = 1540;
    medium.sound_speed_ref  = 1540;
    
    % Density
    % medium.density = patterns(rho_mean, rho_std, 'homo', [], Nx, Ny);
    medium.density          = density_map;
    % medium.density = rho_mean + rho_mean.*rho_std.*randn(Nx, Ny);
    
    % Attenuation
    offset = 5;
    medium.alpha_mode       = 'no_dispersion';
    medium.alpha_power      = 1.05;
        
    if strcmp(typeSimu, 'homo')
    % Homogeneous
        % medium.alpha_coeff = patterns(alpha_mean, alpha_std, 'homo', [], Nx, Ny);
        medium.alpha_coeff  = alpha_value; 
    end
    
    if strcmp(typeSimu, 'circle')
        % Inclusion 
        % alpha_mean is [alpha_mean_bg alpha_mean_inc]
        % alpha_std is [alpha_std_bg alpha_std_inc], default [0 0]
        
        alpha_coeff_map = ones(Nx,Ny,Nz);
        
        alpha_mean = alpha_value;
        alpha_std  = [0 0];
        
        alpha_params.radius_disk = 10e-3;  % [m]
        alpha_params.center_depth = 20e-3; % [m]
        alpha_params.offset = offset;
        alpha_params.dx = dx;
        
        [alpha_map2D, ~, ~] = patterns(alpha_mean,...
                alpha_std, 'circle', alpha_params, Nx, Ny);

        for mm = 1:Nz
            alpha_coeff_map(:,:,mm) = alpha_map2D; 
        end
        medium.alpha_coeff  = alpha_coeff_map;
    end

    % set the input settings
    input_args = {...
        'PMLInside', false, 'PMLSize', [pml_x_size, pml_y_size, pml_z_size], ...
        'DataCast', DATA_CAST, 'DataRecast', true, 'PlotSim', false};
    
    % set medium position
    medium_position = 1;

    if gpuDeviceCount > 0
        sensor_data = kspaceFirstOrder3DG(kgrid, medium, transducer, transducer, input_args{:});
    else
        sensor_data = kspaceFirstOrder3D(kgrid, medium, transducer, transducer, input_args{:});      
    end

    rf_prebf = sensor_data';

    % SAVE DATA FILE
    % Create the filename
    if isscalar(alpha_value) % homo case
        % If alpha_value is a scalar
        filename = sprintf('rf_prebf%d_ang%.1f_a_%.1f', iSim, steering_angle, alpha_value);
    else
        % If alpha_value is a vector
        alpha_str = strjoin(arrayfun(@(x) sprintf('%.1f', x), alpha_value, 'UniformOutput', false), '_');
        filename = sprintf('rf_prebf%d_ang%.1f_a_%s', iSim, steering_angle, alpha_str);
    end

    % Replace '.' with 'p' in the filename
    filename = strrep(filename, '.', 'p');

    % Save file
    save(fullfile(OutputDir, filename+".mat"), "rf_prebf", "fs", "medium", "alpha_value", "c0", "steering_angle");

end
end

%% PATTERNS
function [med, med_mean, med_std] = patterns(var_mean, var_std, patt_type,var_params, Nx, Ny)
    switch patt_type
        case 'homo'
            med_mean = var_mean;
            med_std = var_std;
        case 'circle'
            circle_ind = logical(makeDisc(Nx, Ny,...
                round(var_params.center_depth/var_params.dx)+var_params.offset, Ny/2, ...
                round(var_params.radius_disk/var_params.dx)));
            med_mean = var_mean(1)*(~circle_ind) + var_mean(2)*(circle_ind);
            med_std = var_std(1)*(~circle_ind) + var_std(2)*(circle_ind);
        case 'layers_vert'
            layer_pos = var_params.layer_pos/var_params.dx;
            [X,~] = meshgrid(1:Ny, 1:Nx);
            med_mean = var_mean(1)*(X<layer_pos) ...
                + var_mean(2)*(X>=layer_pos);
            med_std = var_std(1)*(X<layer_pos) ...
                + var_std(2)*(X>=layer_pos);
        case 'layers_horz'
            layer_pos = var_params.layer_pos/var_params.dy;
            [~,Y] = meshgrid(1:Ny, 1:Nx);
            med_mean = var_mean(1)*(Y<layer_pos) ...
                + var_mean(2)*(Y>=layer_pos);

            med_std = var_std(1)*(Y<layer_pos) ...
                + var_std(2)*(Y>=layer_pos);
        case 'layers_vert3'
            layer_pos = var_params.layer_pos/var_params.dx;
            [X,~] = meshgrid(1:Ny, 1:Nx);
            med_mean = var_mean(1)*(X<layer_pos(1)) ...
                + var_mean(2)*((X<layer_pos(2))&(X>=layer_pos(1)))...
                + var_mean(3)*(X>=layer_pos(2));
            med_std = var_std(1)*(X<layer_pos(1)) ...
                + var_std(2)*((X<layer_pos(2))&(X>=layer_pos(1)))...
                + var_std(3)*(X>=layer_pos(2));
        case 'layers_horz3'
            layer_pos = var_params.layer_pos/var_params.dy;
            [~,Y] = meshgrid(1:Ny, 1:Nx);
            med_mean = var_mean(1)*(Y<layer_pos(1)) ...
                + var_mean(2)*((Y<layer_pos(2))&(Y>=layer_pos(1)))...
                + var_mean(3)*(Y>=layer_pos(2));
            med_std = var_std(1)*(Y<layer_pos(1)) ...
                + var_std(2)*((Y<layer_pos(2))&(Y>=layer_pos(1)))...
                + var_std(3)*(Y>=layer_pos(2));
        case 'circleN'
            circle_ind = nan(Nx, Ny, length(var_params.center_depth));
            bg_mask = ones(Nx, Ny);

            med_mean = zeros(Nx, Ny);
            med_std = zeros(Nx, Ny);
            for cc = 1:length(var_params.center_depth)
                circle_ind(:,:,cc) = logical(makeDisc(Nx, Ny,...
                    round(var_params.center_depth(cc)/var_params.dx)+var_params.offset, var_params.center_lat(cc), ...
                    round(var_params.radius_disk(cc)/var_params.dx)));
                med_mean = med_mean + var_mean(cc)*(circle_ind(:,:,cc));
                med_std = med_std + var_std(cc)*(circle_ind(:,:,cc));
                bg_mask = bg_mask & not(circle_ind(:,:,cc));
            end
            med_mean = med_mean + var_mean(end)*(bg_mask);
            med_std = med_std + var_std(end)*(bg_mask);

        otherwise
            med_mean = 0;
            med_std = 0;
            warning('Unexpected type')
    end
    med = med_mean + med_mean.*med_std.*randn(Nx, Ny);
end

