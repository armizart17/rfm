
function [] = simulateLiverRef(baseDir)

addpath(genpath(pwd))
simuNames = {'ref1','ref2','ref3','ref4','ref5','ref6'};
mkdir(baseDir)

% medium parameters
c0              = 1540;     % sound speed [m/s]
rho0            = 1000;     % density [kg/m^3]

% source parameters
source_f0       = 4e6;   % source frequency [Hz]
source_amp      = 1e6;      % source pressure [Pa]
source_cycles   = 3.5;      % number of toneburst cycles
source_focus    = 7.7e-2;   % focal length [m]
element_pitch   = 0.3e-3;   % pitch [m]
element_width   = 0.25e-3;  % width [m]
focal_number    = 4;
nLines          = 64;       % Number of beams;

% grid parameters
grid_size_x     = 11e-2;    % [m]
grid_size_y     = 4e-2;     % [m]

% transducer position
translation     = [-5e-2, 0];
rotation        = 0;

% computational parameters
DATA_CAST       = 'gpuArray-single'; % set to 'single' or 'gpuArray-single'
ppw             = 6;        % number of points per wavelength, 4 to 8
depth           = 10e-2;     % imaging depth [m]
cfl             = 0.3;      % CFL number, could be 0.3 or 0.5
PMLSize         = [43,40];

%% For looping simulations

for iSim = 1:length(simuNames)
    %% GRID

    % calculate the grid spacing based on the PPW and F0
    dx = c0 / (ppw * source_f0);   % [m]
    
    % compute the size of the grid
    Nx = roundEven(grid_size_x / dx);
    Ny = roundEven(grid_size_y / dx);
    
    % create the computational grid
    kgrid = kWaveGrid(Nx, dx, Ny, dx);
    
    % create the time array
    t_end           = depth*2/c0;     % [s];    % total compute time [s]
    kgrid.makeTime(c0, cfl, t_end);

    %% MEDIUM
    rz = kgrid.x - translation(1);
    rx = kgrid.y;

    medium.sound_speed = c0*ones(size(rx));
    medium.density = rho0*ones(size(rx));
    medium.alpha_coeff = 0.5;
    medium.alpha_power = 1;
    medium.alpha_mode = 'no_dispersion';
    medium.sound_speed_ref = c0;
    nScatterers = getNumScatterers(c0,source_f0,focal_number, ...
        source_cycles,grid_size_x,grid_size_y,24);
    medium.sound_speed = generateScattering(medium.sound_speed, ...
        0.05*c0*ones(size(rx)), nScatterers);


    figure('Units','centimeters', 'Position',[5 5 25 10]),
    tiledlayout(1,3),
    nexttile,
    imagesc(100*rx(1,:),100*rz(:,1),medium.sound_speed, c0*[0.9 1.1])
    xlabel('x [cm]'), ylabel('z [cm]')
    title('Sound speed')
    c = colorbar; ylabel(c,'m/s')
    axis image
    
    nexttile,
    imagesc(100*rx(1,:),100*rz(:,1),medium.density)
    xlabel('x [cm]'), ylabel('z [cm]')
    title('Density')
    c = colorbar; ylabel(c,'kg/m^3')
    colorbar,
    axis image
    
    t3 = nexttile;
    imagesc(100*rx(1,:),100*rz(:,1),medium.alpha_coeff, [0.4 1.7])
    xlabel('x [cm]'), ylabel('z [cm]')
    title('Absorption')
    c = colorbar; ylabel(c,'dB/cm/MHz')
    axis image
    colormap(t3,"turbo")

    save(fullfile(baseDir,"refCords.mat"),"rx","rz")

    %% SOURCE
    aperture = source_focus/focal_number;
    element_num = floor(aperture/element_pitch);
    
    % set indices for each element
    ids = (0:element_num-1) - (element_num-1)/2;
    
    % set time delays for each element to focus at source_focus
    time_delays = -(sqrt((ids .* element_pitch).^2 + source_focus.^2) - source_focus) ./ c0;
    time_delays = time_delays - min(time_delays);
    
    % create time varying source signals (one for each physical element)
    source_sig = source_amp .* toneBurst(1/kgrid.dt, source_f0, ...
        source_cycles, 'SignalOffset', round(time_delays / kgrid.dt));
    
    % create empty kWaveArray
    karray = kWaveArray('BLITolerance', 0.05, 'UpsamplingRate', 10);
    
    % add rectangular elements
    for ind = 1:element_num
        
        % set element y position
        y_pos = 0 - (element_num * element_pitch/2 - element_pitch/2) + (ind-1) * element_pitch;
        
        % add element (set rotation angle to match the global rotation angle)
        karray.addRectElement([0, y_pos], element_width/4, element_width, rotation);
    end

    %% For Looping
    yCords = ( (0:nLines-1) - (nLines-1)/2 )* element_pitch; % Lateral cord of each element
    % bf_data_final = zeros(kgrid.Nt ,nLines);
    clear rf_prebf
    for iLine = 1:nLines
        %% SOURCE  
        % move the array
        translation(2) = yCords(iLine);
        karray.setArrayPosition(translation, rotation)
    
        % assign binary mask
        source.p_mask = karray.getArrayBinaryMask(kgrid);
            
        % assign source signals
        source.p = karray.getDistributedSourceSignal(kgrid, source_sig);
        clc
        disp(['Line: ',num2str(iLine),' of ',num2str(nLines)]);
    
        %% SENSOR
        
        % set sensor mask to record central plane
        sensor.mask = karray.getArrayBinaryMask(kgrid);
        sensor.directivity_size = 10*kgrid.dx;
        sensor.directivity_angle = zeros(size(sensor.mask));
        
        %% SIMULATION
        
        % set input options
        input_args = {...
            'PMLInside', false, ...
            'PMLSize', PMLSize, ... 
            'DataCast', DATA_CAST, ...
            'DataRecast', true, ...
            'PlotSim', true};
        
        % MATLAB CPU
        sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, ...
            input_args{:});

                % combine sensor data
        combined_sensor_data = karray.combineSensorData(kgrid, sensor_data);

        % beamforming
        fs = 1/kgrid.dt;        
        rf_prebf(:,:,iLine) = combined_sensor_data';
    end

    %% VISUALISATION
    offset = 1;
    
    axAxis = 0:kgrid.Nt-1; axAxis = axAxis*kgrid.dt*c0/2;
    latAxis = yCords;
    
    % rf = bf_data_final(offset:end,:);
    z = axAxis(offset:end);
    x = latAxis;
    
    save(fullfile(baseDir,['rf_prebf_',simuNames{iSim},'.mat']),...
        'rf_prebf','x','z','fs','time_delays');

end

end




function map = generateScattering(propertyMap, stdMap, nScatterers)
[M,N] = size(propertyMap);
scatterProb = nScatterers/M/N;
scatterMap = rand(M,N)<scatterProb;
randMap = propertyMap + stdMap.*randn(M,N);
map = propertyMap;
map(scatterMap) = randMap(scatterMap);
end

function nScatterers = getNumScatterers(c0,f0,fn,sc,x,y,sprc)
%   c0:     speed of sound      [m/s]
%   f0:     source frequency    [Hz]
%   fn:     focal number
%   sc:     source cycles
%   x,y:    grid size           [m]
%   sprc:   scatterers/res cell [1/m^2]
    wl = c0/f0; % Wavelength
    scatterDensity = sprc/(fn*wl)/(wl*sc/2);
    nScatterers = scatterDensity*x*y;
end