function simuGauss_TUFFC25(alpha_power, alpha_coeff_sam, alpha_coeff_ref, nameFolderOut)
% (1) Homogeneous simulation Gauss sample and 
% (2) classical homogeneous reference (PRELOADED DENSITY MAPS)
% #!/bin/bash
% #SBATCH --output="./pipelinelim/out/slurm-%j.out"
% #SBATCH --time=22:00:00  # 22 hours
% #SBATCH --partition=debug
% #SBATCH --gres=gpu:1
% #SBATCH --mem=MaxMemPerNode # Request maximum available memory
% #SBATCH --ntasks=1
% #SBATCH --cpus-per-task=12
% #SBATCH --exclude=worker1,worker2 # Exclude specific workers
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
% srun matlab -nosplash -nodesktop -nodisplay -r "cd('./codes'); simuGauss_TUFFC25($matlab_args); exit;"



rng("default")

cd('../'); % when you run it usually goes inside ./codes/script.m

rng shuffle;
addpath(genpath(pwd))
addpath(genpath('/opt/MATLAB Add-Ons'))
% DATA_CAST = 'single';     % set to 'single' or 'gpuArray-single' to speed up computations
DATA_CAST = 'gpuArray-single';     % set to 'single' or 'gpuArray-single' to speed up computations

% delete(gcp)
% parpool

c = parcluster('local');  
numWorkers = min(c.NumWorkers, feature('numcores'));  % Use the max allowed
fprintf('Num of workers : %d\n', numWorkers);
delete(gcp('nocreate'));
parpool('local', numWorkers);


gen_sam = true;
gen_ref = true;

file_dens_ref = './densitymaps/TUFFC25/map_sd0p005_i1.mat';
file_dens_sam = './densitymaps/TUFFC25/map_sd0p005_mask9_s1_i1.mat';

pathRF = fullfile(pwd, 'out', nameFolderOut);

if ~exist(pathRF,"dir"); mkdir(pathRF); disp('FolderOut created'); end
%% create the computational grid
if gen_sam

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
rho_std  = 0.005;

% medium.density = patterns(rho_mean, rho_std, 'homo', [], Nx, Ny);
% medium.density = rho_mean .* (1 + rho_std * randn(Nx, Ny));
% medium.density = rho_mean + rho_mean.*rho_std.*randn(Nx, Ny);

% CONVOLUTION **
% kernel_size     = 9;
% kernel_sigma    = 1;
% mask = fspecial('gaussian', [kernel_size kernel_size], kernel_sigma);
% 
% density_pad = rho_mean .* (1 + rho_std * randn(Nx+2*pad, Ny+2*pad));
% density_pad_mask = conv2(density_pad, mask, 'same');
% density_map_sam = density_pad_mask(1+pad:end-pad, 1+pad:end-pad);


load(file_dens_sam); % MAP reference
density_map_sam = density; clear density;

% SAMPLE
medium.density = density_map_sam;

% Attenuation
medium.alpha_power = alpha_power;
medium.alpha_mode = 'no_dispersion';
% medium.alpha_coeff = patterns(alpha_mean, alpha_std, 'homo', [], Nx, Ny);
medium.alpha_coeff = alpha_coeff_sam;

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

% Create the filename
filename = sprintf('rfsam_%.3f', alpha_power);

% Replace '.' with 'p' in the filename
filename = strrep(filename, '.', 'p');

% Save file
save(fullfile(pathRF, filename), "rf", "x", "z", "fs", "medium");
end

if gen_ref 

%% Create the computational grid
kgrid = kWaveGrid(Nx, dx, Ny, dy);


%% Signals and Transducer settings

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
rho_std  = 0.005;

% REFERENCE
% density_map_ref = density_pad(1+pad:end-pad, 1+pad:end-pad,:);
% medium.density = density_map_ref;

load(file_dens_ref); % MAP reference
density_map_ref = density; clear density;

medium.density = density_map_ref;

% Attenuation
medium.alpha_power = alpha_power;
medium.alpha_mode = 'no_dispersion';
% medium.alpha_coeff = patterns(alpha_mean, alpha_std, 'homo', [], Nx, Ny);
medium.alpha_coeff = alpha_coeff_ref;

source_strength = 1e6;
tone_burst_freq = 6.66e6;        % [Hz]
tone_burst_cycles = 3.5;

%% Bucle
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

% Create the filename
filename = sprintf('rfref_%.3f', alpha_power);

% Replace '.' with 'p' in the filename
filename = strrep(filename, '.', 'p');

% Save file
save(fullfile(pathRF, filename), "rf", "x", "z", "fs", "medium");

end

end