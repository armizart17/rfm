function [a_gaus, a_expo, a_fluid, ...
    EAC_gauss, EAC_expo, EAC_fluid, ...
    reg_gaus_k, reg_expo_k, reg_fluid_k, gof] = QUS_estim(BSCcurve_thy, band_k, params)

    if nargin<3
        params.gauss = 1;
        params.expo = 1;
        params.fluid = 1;
        params.plot = 1;
    end
    
    a_gaus  = NaN;
    a_expo  = NaN;
    a_fluid = NaN;
    
    %% ESD estimation
    N     = length(band_k);
    % The maximum radius for the fit, in um:
    max_rad = 200;
    II_max = max_rad - 1;
    a = linspace(1, max_rad, II_max)*1E-6;
    
    FF_Gauss        = @(ka)  exp(-0.827*power(ka, 2));
    sbesselj        = @(z,x) sqrt(pi./(2*x)).*besselj(z+0.5, x);
    FF_fluid        = @(ka)  power(3./(2*ka).*sbesselj(1, 2*ka), 2);
    FF_expo         = @(ka)  power( 1 + 4/(3.63/2)^2*power(ka, 2), -2);
    
    goodnessoffit   = @(x,xref) 1 - (norm(x - xref).^2)/(norm(x - mean(x)).^2);
    
    amatrix         = a'*ones(1, length(band_k));
    kband_matrix    = ones(length(a), 1)*band_k;
    BSC_matrix      = repmat(BSCcurve_thy, length(a), 1);
    
    if params.gauss
        F_gaus      = kband_matrix.^4 .*FF_Gauss(kband_matrix.*amatrix);
        X_gaus      = log10(BSC_matrix./F_gaus);
        Temp_gaus   = var(X_gaus');
    end
    
    if params.expo
        F_expo      = kband_matrix.^4 .*FF_expo(kband_matrix.*amatrix);
        X_expo      = log10(BSC_matrix./F_expo);
        Temp_expo   = var(X_expo');
    end
    
    if params.fluid
        F_flui      = kband_matrix.^4 .*FF_fluid(kband_matrix.*amatrix);
        X_flui      = log10(BSC_matrix./F_flui);
        Temp_flui   = var(X_flui');
    end
    
    if params.gauss
        [~, I_gaus] = min(Temp_gaus);
        if I_gaus == 1 || I_gaus == II_max
            a_gaus = NaN;
        else
            a_gaus = a(I_gaus);
        end
    end
    
    if params.expo
        [~, I_expo] = min(Temp_expo);
        if I_expo == 1 || I_expo == II_max
            a_expo = NaN;
        else
            a_expo = a(I_expo);
        end
    end
    
    if params.fluid
        [~, I_flui] = min(Temp_flui);
        if I_flui == 1 || I_flui == II_max
            a_fluid = NaN;
        else
            a_fluid = a(I_flui);
        end
    end
    
    if params.plot
        figure(123),
        semilogy(a*1e6, Temp_gaus, '-'); hold on;
        semilogy(a*1e6, Temp_expo, '-'); hold on;
        semilogy(a*1e6, Temp_flui, '-'); hold off;
        legend( 'Gaus', 'Expo', 'Flui')
        title('Gaussian form factor');
        xlabel('a_{eff} [micron]')
        ylabel('MSE')
    end
    
    % Now estimating the EAC:
    if isnan(a_gaus) == 0
        Vs          = 4*pi/4*power(a_gaus,3);
        factor      = power(band_k,4).*FF_Gauss(band_k*a_gaus)*Vs^2/(16*pi^2);
        EAC_gauss   = mean(BSCcurve_thy./factor);
        reg_gaus_k  = EAC_gauss*factor;
        gof.gauss   = goodnessoffit(BSCcurve_thy, reg_gaus_k);
    else
        EAC_gauss   = NaN;
        reg_gaus_k  = nan(1,length(band_k));
        gof.gauss   = NaN;
    end
    
    if isnan(a_expo) == 0
        Vs          = 4*pi/4*power(a_expo,3);
        factor      = power(band_k,4).*FF_expo(band_k*a_expo)*Vs^2/(16*pi^2);
        EAC_expo    = mean(BSCcurve_thy./factor);
        reg_expo_k  = EAC_expo*factor;
        gof.expo    = goodnessoffit(BSCcurve_thy, reg_expo_k);
    else
        EAC_expo    = NaN;
        reg_expo_k  = nan(1,length(band_k));
        gof.expo    = NaN;
    end
    
    if isnan(a_fluid) == 0
        Vs          = 4*pi/4*power(a_fluid,3);
        factor      = power(band_k,4).*FF_fluid(band_k*a_fluid)*Vs^2/(16*pi^2);
        EAC_fluid   = mean(BSCcurve_thy./factor);
        reg_fluid_k = EAC_fluid*factor;
        gof.fluid   = goodnessoffit(BSCcurve_thy, reg_fluid_k);
    else
        EAC_fluid   = NaN;
        reg_fluid_k = nan(1,length(band_k));
        gof.fluid   = NaN;
    end
    
    if params.plot
        fprintf('\nFit GAUSS parameter:\n Diameter = %g microns (GoF = %1.3f)\n',2*a_gaus*1e6, gof.gauss);
        fprintf('Fit EXPO parameter:\n Diameter = %g microns (GoF = %1.3f)\n',  2*a_expo*1e6, gof.expo);
        fprintf('Fit FLUID parameter:\n Diameter = %g microns (GoF = %1.3f)\n',  2*a_fluid*1e6, gof.fluid);
    end

end