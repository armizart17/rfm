function [Slope_BSC, Inter_BSC, Midband, lin_curve, Integ_BSC] = BSC_parameters(band, BSC)

    % Calcul of BSC parameters: Slope, Intercept, Intercept and Slope evaluation
    [Sparam, Cparam, ~, ~] = fit_linear(band', log10(BSC') ,1);
    lin_curve = Sparam*band + Cparam;
    
    Slope_BSC = Sparam;
    Inter_BSC = power(10,Cparam); 
    Midband = power(10,lin_curve(round(length(band)/2)));
    Integ_BSC = 1/(length(band)*abs(band(end) - band(1)))* sum(BSC);

end