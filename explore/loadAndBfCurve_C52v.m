function [rf,bMode,xp,zp,z0p,xr,zr,r0,sampleFreq] = loadAndBfCurve_C52v(fileName, presetName, ax)
    [~, ~, ext] = fileparts(fileName);
    % Whether is a '.vrs' or '.mat' recover the RF Channel Data
    if (ext == ".vrs")
        channelData = vsv.file.readVRSFile(fileName);
        
    elseif (ext == ".mat")
        channelData = load(fileName);
        channelData = cell2mat(channelData.RcvData);
    end

    presetData = load(presetName);
    presetData = presetData.preSet;
    
    Trans = presetData.Trans;
    Receive = presetData.Receive;
    P = presetData.P;
    txNum = presetData.Refocus.TXNum;
    RcvZone = presetData.Refocus.RcvZone;

    % Check if the field exists before assignment
    if isfield(presetData, 'Refocus') && isfield(presetData.Refocus, 'TXNum')
        txNum = presetData.Refocus.TXNum; % Safe assignment
    elseif isfield(presetData, 'TXNum')
        txNum = presetData.TXNum; % Safe assignment
    end

    centralFreq = Receive(1).demodFrequency*1e6; % Central freq. [Hz]
    sampleFreq = Receive(1).decimSampleRate*1e6; % Sample freq. [Hz]

    nPulses = P.numRays;
    endDepth = P.endDepth;
    nElements = Trans.numelements;
    nSamples = Receive(1).endSample - Receive(1).startSample + 1;
    soundSpeed = presetData.speedOfSound; % [m/s]
    wvl = soundSpeed/centralFreq;


    % Aux variables
    t = (0:(nSamples-1))/sampleFreq; % [sec.] time domain 0:T:(N_sample-1)*T
    z = soundSpeed*t/2*1e3;
    elementPosX = Trans.ElementPos(:, 1)*wvl; % x [m] X center of Trans.element
    elementPosZ = Trans.ElementPos(:, 3)*wvl; % z [m] Z center of Trans.element
    phi = Trans.ElementPos(:, 4); % phi [rad] angle for every Trans.element
    focus = soundSpeed*t(:)/2;
    tgcDbCm = 0.1 * (centralFreq * 1e-6) ^ 1;
    tgcNpM = tgcDbCm / 8.686 * 100;
    r = exp(tgcNpM * z*1e-3); % [m]

    % Get reception delays
    rxDelays = getRXDelays(elementPosX, elementPosZ, phi, focus, nElements, nPulses, soundSpeed);
    rxSamples = round(rxDelays*sampleFreq);
    rxSamples(rxSamples<1) = 1;
    rxSamples(rxSamples>nSamples) = nSamples;

    % Apodization
    maxAperSize = 65;
    apodAperture = getApodCoeff(maxAperSize, nElements, nPulses);

    % Organize and BMode
    rfChannel = zeros(nSamples, nElements, nPulses);
    plotDepth = floor(endDepth*2*4);
    bMode = zeros(plotDepth, nPulses); 
%     rfChannel(:, :, :) = getOrderedData(channelData, Receive, nSamples, nElements, nPulses, 1);

    try
        % Intentar con la asignación de tres dimensiones
        for i = 1:txNum
            rfChannel(:, :, :, i) = getOrderedData(channelData, Receive, nSamples, nElements, nPulses, i);
        end
        
        samples2Recon = [nSamples, RcvZone*2*4, 1];
        samples2Recon(2) = [];
        for i = 1:txNum
            rfChannel(samples2Recon(i+1):samples2Recon(i), :, :, end) = rfChannel(samples2Recon(i+1):samples2Recon(i), :, :, i);
        end
        rfChannel = rfChannel(:, :, :, end);
%         rfChannel(:, :, :) = getOrderedData(channelData, Receive, nSamples, nElements, nPulses, 1);
    catch
        % Si hay un error, probar con la otra asignación
        rfChannel = getOrderedDataAnt(channelData(:, :), Receive, nSamples, nElements, nPulses, txNum);
    end
    
    [bMode, rf] = getBModeC52v(rfChannel(:, :, :), rxSamples, apodAperture, nSamples, nPulses, nElements, plotDepth, r);
   
    % For polar grid
    param = getparamC52v();
    z2 = (0:size(bMode, 1) - 1)/sampleFreq*soundSpeed/2; % [m]
    
    % CS
    % [xp, zp] = toPolarGrid(size(bMode,[1,2]), z2(end), param);  
    % EMZ
    % [xp,zp,z0p] = impolgrid(size(bMode,[1,2]), z2(end), param);
    
    % vLast
    [xp,zp,z0p,xr,zr,r0] = toPolarGridRec(size(bMode,[1,2]), z2(end), param);  

    zp = zp - Trans.radiusMm*1e-3*(1-cos(phi(1)));
    
    [plane, patientCode, frameN] = getAcqData(fileName);
    % plotBMode(bMode(:, :, end), ax, xp, zp, plane, patientCode, frameN);

    rectangleROI = gobjects(0);
    rectangleCount = 0;
    assignin('base', 'rectangleROI', rectangleROI);
    assignin('base', 'rectangleCount', rectangleCount);
end

%% Aux Functions
function plotBMode(Image, axes, xp, zp, plane, patientCode, frameN)
        pcolor(axes, xp*1000, zp*1000, Image);
        axes.Color = [0 0 0];
        axes.FontSize = 13;
        cbar = colorbar;
        ylabel(cbar, '[dB]');
        %cbar.Ticks = [];
        clim([-55, 0]);
        colormap gray
        fontSize = 18;
        caption = append('Paciente: ', patientCode, ' - Plano: ', plane, ...
            ' - Fotograma: ', frameN);
        title(caption, 'FontSize', fontSize);
        ylabel('[mm]', 'FontSize', 12);
        xlabel('[mm]', 'FontSize', 12);
        shading interp
        axis equal ij tight
end

% Dynamic Focusing Delays
function rx_delays = getRXDelays(element_pos_x, element_pos_z, phi, focus, n_elements, n_pulses, sound_speed)
    rx_delays = zeros(length(focus), n_elements, n_pulses); % length(focus) is same as length(t)
    for n = 1:n_pulses
        xfocus = element_pos_x(n) + sin(phi(n)) * focus;
        zfocus = element_pos_z(n) + cos(phi(n)) * focus;
        for e = 1:n_elements
            rx_delays(:,e,n) = (focus + ...
                sqrt((zfocus- element_pos_z(e)).^2 + (xfocus - element_pos_x(e)).^2))/sound_speed;
        end
    end
end

% Dynamic aperture matrix
function  apodAperture = getApodCoeff(maxAperSize, nElements, nPulses)
    apodAperture = zeros(nElements, nPulses);
    for i = 1:nPulses
        aperCenter = i;
        halfAperSize = floor(maxAperSize/2);
        aIndex = -halfAperSize:halfAperSize;
        aperture = aperCenter + aIndex;
        aperture = aperture(aperture>=1);
        aperture = aperture(aperture<=nElements);
        apodAperture(aperture, i) = 1;
    end
end

% Organize data
function rf_channel = getOrderedData(BuffData, Receive, n_samples, n_elements, n_pulses, k)
    rf_channel = zeros(n_samples, n_elements, n_pulses);
    for n = 1:n_pulses  
        % Select RF Channel Data from every Buffer
        rf_channel(:, :, n) = BuffData(Receive(n + n_pulses*(k-1)).startSample:Receive(n + n_pulses*(k-1)).endSample, :); % RF Data from Buffer
    end
end

% Organize data Ant
function rf_channel = getOrderedDataAnt(BuffData, Receive, n_samples, n_elements, n_pulses, k)
    rf_channel = zeros(n_samples, n_elements, n_pulses);
    for n = 1:n_pulses  
        % Select RF Channel Data from every Buffer
        rf_channel(:, :, n) = BuffData(Receive(k*(n-1)+1).startSample:Receive(k*(n-1)+1).endSample, :); % RF Data from Buffer
    end
end

% Get BMode & RF @EMZ
function [bMode, rf_data] = getBModeC52v(rf_channel, rx_samples, apodAperture, n_samples, n_pulses, n_elements, plotDepth, tgcr)
    rf_data = zeros(n_samples, n_pulses);
    for n = 1:n_pulses
        % Delaying
        for e = 1:n_elements
            rf_channel(:, e, n) = rf_channel(rx_samples(:, e, n), e, n);
        end
        % Apodization
        rf_channel(:, :, n) = rf_channel(:, :, n) .* apodAperture(n, :);
        % Summing
        rf_data(:, n) = sum(rf_channel(:, :, n), 2);
    end

    % Depth segmentation
    rf_data = rf_data(1:plotDepth, :);

    % TGC
    tgcr = tgcr(1:plotDepth);
    rf_data = bsxfun(@times, tgcr', rf_data);

    rf_data_padded = zeros(floor(size(rf_data, 1)*1.5), size(rf_data, 2));
    rf_data_padded(1:size(rf_data, 1), :) = rf_data;

    % BMode
    bMode = abs(hilbert(rf_data_padded));
    bMode = 20*log10(bMode);
    bMode = bMode(1:size(rf_data ,1), :);
    bMode = bMode-max(bMode(:));
end

% Get Acquisition Plane
function [plane, patientCode, frameN] = getAcqData(fileName)
    [~, fName, ~] = fileparts(fileName);
    fParts = strsplit(fName, '_');
    frameN = '1';
    patientCode = fParts{1};
    plane = getPlane(fParts{2});
    if (length(fParts) > 3)
        frameN = fParts{end};
    end       
end


function acqPlane = getPlane(ac_plane)
    switch ac_plane
        case 'LCM'
            acqPlane = 'Línea Clavicular Media';
        case 'IOLAI'
            acqPlane = 'Intercostal Oblicuo Línea Axilar Inferior';
        case 'IHR'
            acqPlane = 'Interfaz Hepatorrenal';
        case 'LHI'
            acqPlane = 'Lóbulo Hepático Izquierdo';
        case 'PL'
            acqPlane = 'Plano Libre';
        otherwise
           acqPlane = 'Plano Legado (antiguo)';
    end
end

% C5-2v Parameters (Verasonics)
function param = getparamC52v()
        param.fc = 3.57e6; % Transducer center frequency [Hz]
        param.kerf = 48e-6; % Kerf [m]
        param.width = 460e-6; % Width [m]
        param.pitch = 508e-6; % Pitch [m]
        param.Nelements = 128;
        param.bandwidth = 79; % Fractional bandwidth [%]
        param.radius = 49.57e-3; % Array radius [m]
        param.height = 13.5e-3; % Elevation height [m]
        param.focus = 60e-3; % Elevation focus [m]
end

function [xPolar, zPolar] = toPolarGrid(siz,zmax,param)                     
    N = param.Nelements;
    R = param.radius;
    p = param.pitch;
    L = 2*R*sin(asin(p/2/R)*(N-1)); % chord length
    d = sqrt(R^2-L^2/4); % apothem
    z0 = -d;
    [th,r] = meshgrid(...
        linspace(atan2(L/2,d),atan2(-L/2,d),siz(2))+pi/2,...
        linspace(R+p,-z0+zmax,siz(1)));
    [xPolar,zPolar] = pol2cart(th,r);

    zPolar = zPolar+z0;
end

function [xPolar,zPolar,z0,xFull,zFull,r0] = toPolarGridRec(siz,zmax,param)                     
    N = param.Nelements;
    R = param.radius;
    p = param.pitch;
    L = 2*R*sin(asin(p/2/R)*(N-1)); % chord length
    d = sqrt(R^2-L^2/4); % apothem
    z0 = -d;
    [th,r] = meshgrid(...
        linspace(atan2(L/2,d),atan2(-L/2,d),siz(2))+pi/2,...
        linspace(R+p,-z0+zmax,siz(1)));
    [xPolar,zPolar] = pol2cart(th,r);

    zPolar = zPolar+z0;

    % NEW EMZ
    th_r = -(linspace(atan2(L/2,d),atan2(-L/2,d),siz(2)))*180/pi;
    r_r = linspace(R+p,-z0+zmax,siz(1));
    
    xFull = th_r; % [deg]
    r0 = r_r(1);
    zFull = (r_r-r0); % [m]

end
