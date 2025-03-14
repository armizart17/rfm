function [ACS, SLD_lin] = fitPix_ACS(SLD_og, SLD_den)
% function [ACS] = fitPix_ACS(SLD_og, SLD_den)
%   Detailed explanation goes here     
    [M, N, P] = size(SLD_og.SLD_term);
    ACS = zeros(M,N);
    SLD_lin = zeros(M,N,P);
    for u = 1 : M
        for v = 1 : N
            % Denoised
            SLD_vect_den = squeeze(squeeze(SLD_den(u,v,:)));
            [beta_reg , tau_reg, y_lin_reg, R_reg] = fit_linear(SLD_og.band(:),SLD_vect_den, 2); 
            clear SLD_vect_den;
            ACS(u,v) = - beta_reg/(4*SLD_og.zp_zd)*8.6858; % [dB/MHz/cm]
            SLD_lin(u,v,:) = y_lin_reg;
        end
    end

end