function comp_factor = comp_ref_a(delta_a, j_sam, band, depth, q)
% function comp_factor = comp_ref_a(delta_a, j_sam, band, depth, q)
    Np2dB = 20*log10(exp(1)); % 8.6859
    comp_factor = exp(-4*delta_a*100*repmat((band.^j_sam)*depth',[1,1,q])/Np2dB);
end