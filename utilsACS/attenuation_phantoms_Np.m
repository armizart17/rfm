% FREQ : integer or a vector
function [ATT] = attenuation_phantoms_Np(FREQ, CHOICE, SLOPE)
% BASED ON:
% function [ATT] = attenuation_phantoms_Np(FREQ, CHOICE)
% MODIFICATION FOR DATA UIUC (HAN, AIGUO) returns in Np/cm
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
    case 5 %high attenuation phantom using Madsen values - spline interpolation
        madsenHighAtt=[1.2000,2.7000,4.4600,6.3600,7.9800,10.2800];
        madsenFreqs=[2.5,5,7.5,10,12,15];
        ATT=interp1(madsenFreqs,madsenHighAtt./8.6858,FREQ,'spline');
    case 6 %new low att phantom using Madsen values
        madsenLowAtt=[0.076,0.327,0.729,1.296,2,3.345,5.411];
        madsenFreqs=[1,2.5,5,7.5,10,12,15];
        ATT=interp1(madsenFreqs,madsenLowAtt./8.6858,FREQ,'spline');

    % PHANTOM UIUC V1 (Han, Aiguo)
    case 20  % P20 phantom
       
        alpha = 0.551509191332541;
        n = 1.107857152971458;
        beta = 0.208884829068225; 

        ACC  = alpha*FREQ.^n + beta;    % [dB/cm]   
        ATT = ACC/8.6858; % [Np/cm] 

    case 22 % P22 phantom

        alpha = 0.222127445307033;
        n = 1.151897288058939;
        beta = 0.064230439662590;  

        ACC  = alpha*FREQ.^n + beta;    % [dB/cm]   
        ATT = ACC/8.6858; % [Np/cm] 
    
    case 111
        % phantom "linear dependency"
        ATT  = SLOPE*FREQ /8.6858;             % [Np/cm]
    case 222
        % phantom "power dependency" (COILA EXVIVO)
        a = 0.28; m = 1.34; % 6 BACKGROUND
        ATT  = (a*FREQ.^m)/8.6858;             % [Np/cm]

end