% check emz


% Script to simulate a 3D linear array. Adapted from k-wave example on
% how to setup a linear array transducer using the kWaveArray class.


clearvars; clc
addpath(genpath(pwd));

% =========================================================================
% DEFINE LITERALS
% =========================================================================

% select which k-Wave code to run
%   1: MATLAB CPU code
%   2: MATLAB GPU code
%   3: C++ code
%   4: CUDA code
model           = 4;

% medium parameters
c0              = 1540;     % sound speed [m/s]
rho0            = 1000;     % density [kg/m^3]

% source parameters
source_f0       = 6.66e6;      % source frequency [Hz]
source_amp      = 1e6;      % source pressure [Pa]
source_cycles   = 3.5;      % number of toneburst cycles
element_num_x   = 32;       % number of elements
element_num_y   = 8;        % number of elements
element_width   = 0.25e-3;  % width [m]
element_length  = 0.25e-3;  % elevation height [m]
element_pitch   = 0.3e-3;   % pitch [m]

% grid parameters
grid_size_x     = 3.2e-2;    % [m]
grid_size_y     = 1.2e-2;    % [m]
grid_size_z     = 1.2e-2;    % [m]

% transducer position
translation     = [-1.5e-2, 0, 0];
rotation        = [0, 90, 0];

% computational parameters
ppw             = 6;        % number of points per wavelength
cfl             = 0.3;      % CFL number

yPositionArray = (-1.5:1:1.5)*element_num_y*element_pitch;

%% GRID
% --------------------

% calculate the grid spacing based on the PPW and F0
dx = c0 / (ppw * source_f0);   % [m]

% compute the size of the grid
Nx = roundEven(grid_size_x / dx);
Ny = roundEven(grid_size_y / dx);
Nz = roundEven(grid_size_z / dx);

% create the computational grid
kgrid = kWaveGrid(Nx, dx, Ny, dx, Nz, dx);

% create the time array
depth           = grid_size_x/2 - translation(1);
t_end           = depth*2/c0;     % [s];    % total compute time [s]
kgrid.makeTime(c0, cfl, t_end);

%% MEDIUM
% --------------------

cx = 0; cz = 1.5e-2; rInc = 0.45e-2;
acsBack = 0.5; acsInc = 1;
stdBack = 1/100; stdInc = 4/100;

rz = kgrid.x - translation(1);
rx = kgrid.z;
incMask = (rx - cx).^2 + (rz - cz).^2 < rInc^2;
scatterMap = stdBack*randn(size(rx));
scatterMap(incMask) = scatterMap(incMask)/stdBack*stdInc;
attenuationMap = acsBack*ones(size(rx));
attenuationMap(incMask) = acsInc;

% assign medium properties
medium.sound_speed = c0;
medium.density = rho0*(1 + scatterMap);
medium.alpha_coeff = attenuationMap;
medium.alpha_power = 1.1;
% medium.alpha_power = 1.0;
% medium.alpha_mode = 'no_dispersion'; 
medium.sound_speed = c0*ones(size(rx));
clear attenuationMap scatterMap rx rz incMask

%% Looping transducer positions
for ii = 1:4
    %% SOURCE
    % create time varying source signals (one for each physical element)
    source_sig = source_amp .* toneBurst(1/kgrid.dt, source_f0, source_cycles, ...
        'SignalOffset', zeros(element_num_y*element_num_x,1));

    % create empty kWaveArray
    karray = kWaveArray('BLITolerance', 0.05, 'UpsamplingRate', 10);

    % Axis for B-mode coordinates
    xElem = zeros(element_num_x,1); yElem = zeros(element_num_y,1);

    % add rectangular elements in B-mode cords (xy)
    for iy = 1:element_num_y
        y_pos = 0 - (element_num_y * element_pitch / 2 - element_pitch / 2) + (iy - 1) * element_pitch;
        for ix = 1:element_num_x
            x_pos = 0 - (element_num_x * element_pitch / 2 - element_pitch / 2) + (ix - 1) * element_pitch;
            % add element (set rotation angle to match the global rotation angle)
            karray.addRectElement([x_pos, y_pos, 0], element_width, element_length, rotation);

            xElem(ix) = x_pos;
        end
        yElem(iy) = y_pos;
    end

    % move the array
    translation(2) = yPositionArray(ii);
    yElem = yElem + yPositionArray(ii);
    karray.setArrayPosition(translation, rotation)

    % assign binary mask
    source.p_mask = karray.getArrayBinaryMask(kgrid);

    % plot the off-grid source mask
    voxelPlot(single(source.p_mask));
    title('Off-grid source mask');


    % assign source signals
    source.p = karray.getDistributedSourceSignal(kgrid, source_sig);

    %% SENSOR
    % set sensor mask to record central plane
    sensor.mask = karray.getArrayBinaryMask(kgrid);

    %% SIMULATION
    % --------------------

    % set input options
    input_args = {...
        'PMLSize', 'auto', ...
        'PMLInside', false, ...
        'PlotPML', false, ...
        'DisplayMask', 'off'};

    % run code
    switch model
        case 1
            % MATLAB CPU
            sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, ...
                input_args{:}, ...
                'DataCast', 'single', ...
                'PlotScale', [-1, 1] * source_amp/10);

        case 2
            % MATLAB GPU
            sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, ...
                input_args{:}, ...
                'DataCast', 'gpuArray-single', ...
                'PlotScale', [-1, 1] * source_amp/10);

        case 3
            % C++
            sensor_data = kspaceFirstOrder3DC(kgrid, medium, source, sensor, input_args{:});

        case 4
            % C++/CUDA GPU
            sensor_data = kspaceFirstOrder3DG(kgrid, medium, source, sensor, input_args{:});
    end

    % combine sensor data
    rf_prebf = karray.combineSensorData(kgrid, sensor_data);
    rf_prebf = reshape(rf_prebf,element_num_x,element_num_y,size(rf_prebf,2));

    %% SAVE DATA
    fs = 1/kgrid.dt;
    t0 = source_cycles/source_f0/2;
    save("rf_prebf_pos"+ii+".mat",...
        'rf_prebf','xElem','yElem','fs','t0');
end