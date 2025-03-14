function [freqs, bsc] = duke_bsc_band(c0, band, beads_per_mm3, anderson_coeffs)
    freqs = band;
    k0s = 2*pi*freqs/c0;
    
    % The scatterer radii in mm:
    a_values = (38:1:45)/2000; %um a mm, r -> d
    % The radii pdf:
    n_values = ones(size(a_values));
    n_values = n_values/sum(n_values);
    N_vals = length(a_values);
    
    bsc = zeros(length(band), 1);
    
    for index = 1:N_vals
        bsc = bsc + n_values(index)*pi*(a_values(index))^2*ppval(anderson_coeffs, k0s*a_values(index));
    end
    
    scats_per_mm3 = beads_per_mm3/4/pi;
    bsc = scats_per_mm3*bsc;
    
    % Finally, converting to cm^-1:
    bsc = 10*bsc;
end