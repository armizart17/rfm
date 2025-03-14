function [ATT] = attenuation_phantoms_Np(FREQ, CHOICE)
% function [ATT] = attenuation_phantoms_Np(FREQ, CHOICE)
% FREQ : integer or a vector

    switch CHOICE
        case 1
            % phantom "low attenuation"
            ATT  = (0.0141*FREQ.^2 + 0.0838*FREQ - 0.0662)/8.6858;    % [Np/cm]
        case 2
            % phantom "high attenuation" 
            ATT = (0.0085*FREQ.^2 +0.5866*FREQ -0.3725)/8.6858;      % [Np/cm]
        case 3
            % phantom 3 - "Agar"
            ATT =  (0.0076*FREQ.^2 + 0.1189*FREQ -0.0319)/8.6858;        % [Np/cm]
        case 4
            % phantom 4 - "Agar & Milk"
            ATT =  (0.0057*FREQ.^2 + 0.4432*FREQ -0.1000)/8.6858;        % [Np/cm]
        case 5 % high attenuation phantom using Madsen values 
            madsenHighAtt = [1.2000,2.7000,4.4600,6.3600,7.9800,10.2800];
            madsenFreqs = [2.5,5,7.5,10,12,15];
            % spline interpolation
            ATT = interp1(madsenFreqs,madsenHighAtt./8.6858,FREQ,'spline');
        case 6 % new low att phantom using Madsen values
            madsenLowAtt = [0.076,0.327,0.729,1.296,2,3.345,5.411];
            madsenFreqs = [1,2.5,5,7.5,10,12,15];
            % spline interpolation
            ATT = interp1(madsenFreqs,madsenLowAtt./8.6858,FREQ,'spline');
    end

end
