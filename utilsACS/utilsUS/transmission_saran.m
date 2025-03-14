function t = transmission_saran(band)

    c_s = 2400; %[m/s]
    c_p = 1540; %[m/s]
    rho_s = 1690; %[kg/m^3]
    rho_p = 1000; %[kg/m^3]
    Z_p = c_p*rho_p;
    Z_s = c_s*rho_s;
    
    f0 = band*1e6;  %[Hz]
    alpha = 5*(band).^1.5; %[Np/MHz^1.5/m]
    L = 25.1e-6; %[m]
    
    
    t = abs(     2*Z_p ./...
        (    (2*Z_p)*cos((2*pi*f0/c_s - 1i*alpha)*L) ...
        + 1i *(Z_s+Z_p^2/Z_s)*sin((2*pi*f0/c_s - 1i*alpha)*L) )) .^4;

end

% function t = transmission_saran(band)
% a = -0.087146134643529;
% b = 0.642790369805658;
% c = 1.052562238394482;
% t = (a*band.^(b)+c).^4;
% end