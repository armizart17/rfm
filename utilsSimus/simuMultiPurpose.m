function simuMultiPurpose(alpha_value, numSim, nameFolderOut, typeSimu)
% function simuMultiPurpose(alpha_value, numSim, nameFolderOut, typeSimu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Classical homogeneous reference approx size widht (3.8cm), height (5.47cm)
% INPUT
%       alpha_value     : is for alpha_attenuation, when inclusion [alpha_bg alpha_inc]
%       numSim          : number of simulations
%       nameFolderOut   : name of folder to save (It will save data in ./out/nameFolderOut/
%       typeSimu        : type of simu (i.e. 'homo' and 'circle' for now, 
%                           TBD:'layers', check patterns function)
% OUTPUT              
% SAVED DATA in ./out/nameFolderOut/rf(iSim)_a_(alpha_bg)_(alpha_inc).mat'
% By EMZ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FOR CLUSTER ("simuMultiPurpose.sh" file)
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
% srun --ntasks=1 --nodes=1 matlab -nosplash -nodesktop -nodisplay -r "c=parcluster('local'); numWorkers=min(c.NumWorkers, feature('numcores')); maxNumCompThreads(numWorkers); delete(gcp('nocreate')); parpool('local', numWorkers); simuMultiPurpose($matlab_args); exit;"


% cd('../'); % when you run it usually goes inside ./codes/script.m
addpath(genpath(pwd))
addpath(genpath('/opt/MATLAB Add-Ons'))

list_rho_mean = [1000 1000 1000 1000 1000];
list_rho_std = [0.01 0.02 0.03 0.04 0.05];
% list_rho_std = [0.01 0.02 0.03 0.04 0.05];

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

clearvars -except iSim alpha_value numSim nameFolderOut list_rho_mean list_rho_std OutputDir typeSimu
clc; rng shuffle;

% DATA_CAST = 'single';     % set to 'single' or 'gpuArray-single' to speed up computations
DATA_CAST = 'gpuArray-single';     % set to 'single' or 'gpuArray-single' to speed up computations


%% Create the computational grid
elem_pitch = 0.30e-3;

scale_fs = 2;
pml_x_size = scale_fs*40;                % [grid points]
pml_y_size = scale_fs*14;                % [grid points]

% set total number of grid points not including the PML
Nx = scale_fs*810 - 2 * pml_x_size;      % [grid points]
Ny = scale_fs*540 - 2 * pml_y_size;    % [grid points] ;% /4 for going from 128 to 41 aperture

PML_size = [pml_x_size pml_y_size];       % size of the PML in grid points TONY

element_width = 4*scale_fs;             % 4 for at least 0.075mm
ratio = element_width;
dy = elem_pitch/element_width;          % grid point spacing in the y direction [m]
dx = dy;

kgrid = kWaveGrid(Nx, dx, Ny, dy);

Nx_tot = Nx;
Ny_tot = Ny;
rx = ones(Nx_tot,1)*linspace(-Ny_tot*dy/2,Ny_tot*dy/2,Ny_tot);
rz = linspace(0,Nx_tot*dx,Nx)'*ones(1,Ny_tot);

%% Signals and Transducer settings
c0 = 1540;
t_end = (Nx*dx)*2/c0;     % [s]
kgrid.makeTime(c0, [], t_end);
fs = 1/kgrid.dt;

%[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed,0.40875,1.2*Nx*dx*2/1540);

focal_distance = 20e-3;   % center of circle
focal_number = 2;
nAperture = (focal_distance/focal_number)/dy;
nAperture = ratio*floor(nAperture/ratio);
nApertureEle =  nAperture/ratio;

nAperture = nApertureEle*ratio;

nLines = floor(Ny/ratio); % vary slightly plm_y to get nLines=128
% nLines = 128;

bf_data_final = nan(kgrid.Nt ,nLines);

%% Parameters
offset = 5;
rho0 = 1000;

% SoS
% medium.sound_speed = patterns(sos_mean, sos_std, 'homo', [], Nx, Ny);
medium.sound_speed = 1540;
medium.sound_speed_ref = 1540;

% Density
rho_mean = 1000;
rho_std  = 0.02;
% medium.density = patterns(rho_mean, rho_std, 'homo', [], Nx, Ny);
medium.density = rho_mean .* (1 + rho_std * randn(Nx, Ny));
% medium.density = rho_mean + rho_mean.*rho_std.*randn(Nx, Ny);

% Attenuation
offset = 5;

if strcmp(typeSimu, 'homo')
% Homogeneous
% medium.alpha_coeff = patterns(alpha_mean, alpha_std, 'homo', [], Nx, Ny);
medium.alpha_coeff = alpha_value; 
end

if strcmp(typeSimu, 'circle')
% Inclusion 
% alpha_mean is [alpha_mean_bg alpha_mean_inc]
% alpha_std is [alpha_std_bg alpha_std_inc], default [0 0]

alpha_mean = alpha_value;
alpha_std  = [0 0];

alpha_params.radius_disk = 10e-3;  % [m]
alpha_params.center_depth = 20e-3; % [m]
alpha_params.offset = offset;
alpha_params.dx = dx;

[medium.alpha_coeff, alpha_coeff_mean, alpha_coeff_std] = patterns(alpha_mean,...
        alpha_std, 'circle', alpha_params, Nx, Ny);
end

medium.alpha_power = 1.1;
medium.alpha_mode = 'no_dispersion';

source_strength = 1e6;
tone_burst_freq = 6.66e6;        % [Hz]
tone_burst_cycles = 3.5;

%%
for ii = 1:nLines

    jj = ratio*ii;
    axis_x = rz(:,1);
    axis_y = rx(1,:);

    disp(['Lines: ',num2str(ii),' de ',num2str(nLines)]);
    src_ini = max(1,jj - nAperture/2) ;
    src_fin = min(Ny,jj + nAperture/2-1) ;
    % src_ini = round(Ny/2 - nAperture/2) ;
    % src_fin = round(Ny/2 + nAperture/2-1) ;

    [temp,pos] = min(abs(axis_x-focal_distance));
    focus_point = [axis_y(jj) axis_x(pos)];

    aperture_point_src = [axis_y(src_ini:src_fin)' axis_x(offset)*ones(src_fin-src_ini+1,1)];
    
%         figure (6); plot(aperture_point_src(:,1),aperture_point_src(:,2),'sb');
%         hold on; plot(focus_point(:,1),focus_point(:,1),'sr'); hold off;

%% Define Sensor properties

    sensor.mask = zeros(Nx, Ny);
    % Need a slight offset here to prevent backwards propagating pulse
    sensor.mask(offset, src_ini:ratio:src_fin) = 1;
    %sensor.directivity_size = 0.2698e-3;
    %sensor.directivity_angle = zeros(size(sensor.mask));
    %sensor.directivity_size = 0.4e-3;
    sensor.directivity_size = 10*kgrid.dx;
    sensor.directivity_angle = zeros(size(sensor.mask));

%% Define Source parameters

    %amp = 100000; % [au]
    source.p_mask = zeros(Nx, Ny);
    source.p_mask(offset,src_ini:src_fin) = 1;
    apod = boxcar(nnz(source.p_mask));

    input_signal_norm = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);
    input_signal = (source_strength ./ (c0 * rho0)) .* input_signal_norm;

    %excitation = amp*toneBurst(1/kgrid.dt,6.66*1e6,5);
    %src_exc = apod(:)*excitation(:)';
    src_exc = apod(:)*input_signal(:)';

    angle = 0;
    %source.p1 = steer_delay(src_exc,angle,dy,dt,1540);
    %figure; imagesc(src_exc);
    %figure; imagesc(source.p1);
    %source = steer_delay_focus(src_exc,angle,dy,kgrid.dt,1540,aperture_point_src,focus_point);
    %source3 = steer_delay_focus(src_exc,angle,dy,kgrid.dt,1540,aperture_point_src,[0 Inf]);
    source.p = src_exc;
    %figure; imagesc(source.p2);
    %source.p = source.p2;

    
%% Set the input arguments:
% Force the PML to be outside the computational grid; switch off p0 smoothing 
% within kspaceFirstOrder2D

    %PML_alpha = 2;   % Default is 2
    % input_args = {'PMLInside', false, 'PMLAlpha', PML_alpha, 'PMLSize', PML_size, 'PlotPML', false,...
       % 'Smooth', false,'PlotSim',false, 'DataCast',DATA_CAST, 'DataRecast', true};
    input_args = {...
        'PMLInside', false, ...
        'PMLSize', PML_size, ... %'PMLAlpha', PML_alpha, ...
        'DataCast', DATA_CAST, ...
        'DataRecast', true, ...
        'PlotSim',false};

    % run the simulation
    % colormap gray
    sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

    sensor_data = sensor_data';

%% Beamform
    max_apert = 64; % elements
    f_number = 2;
    c_bf = 1540;
    bf_data = BFangle(sensor_data,max_apert,fs,c_bf,elem_pitch,'rect',f_number,0);

    if src_ini <= 1
        index = size(bf_data,2) - floor(nApertureEle/2);
    elseif src_fin>= Ny
        index = floor(nApertureEle/2)+1;
    else
        index = floor(nApertureEle/2)+1;
    end

    bf_data_final(:,ii) = bf_data(:,index);

end

%% Plot BF data

normZero = @(x) x-max(x(:));
rf2Bmode = @(x) 20*log10(abs(hilbert(x)));

offset_plot = 100;

axAxis = 0:size(bf_data_final,1)-1; axAxis = axAxis*1/fs*c0/2;
latAxis = 0:size(bf_data_final,2)-1;  latAxis = latAxis-mean(latAxis); latAxis = latAxis *elem_pitch;

rf  = bf_data_final(offset_plot:end,:);
z   = axAxis(offset_plot:end);
x   = latAxis;

% figure; imagesc(x*1e3, z*1e3, normZero(rf2Bmode(rf)),[-60 0])
% title('B-mode'); colormap gray
% xlabel('Lateral distance (mm)')
% ylabel('Depth (mm)')
% axis image

% Create the filename
if isscalar(alpha_value) % homo case
    % If alpha_value is a scalar
    filename = sprintf('rf%d_a_%.1f', iSim, alpha_value);
else
    % If alpha_value is a vector
    alpha_str = strjoin(arrayfun(@(x) sprintf('%.1f', x), alpha_value, 'UniformOutput', false), '_');
    filename = sprintf('rf%d_a_%s', iSim, alpha_str);
end

% Replace '.' with 'p' in the filename
filename = strrep(filename, '.', 'p');

% Save file
save(fullfile(OutputDir, filename+".mat"), "rf", "x", "z", "fs", "medium", "alpha_value");
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

