function comp_factor = comp_mod_freq_a(a_ref, j_sam, j_ref, band, depth, q)
% function comp_factor = comp_mod_freq_a(a_ref, j_sam, j_ref, band, depth, q)
% useful gamma corrector when m_s != m_r,
% from alpha_s*f^m_s, alpha_r*f^m_r. Otherwise is 1
% see Paper DoF 
    Np2dB = 20*log10(exp(1)); % 8.6859
    comp_factor = exp(4*a_ref*100*repmat((band.^j_sam-band.^j_ref)*depth',[1,1,q])/Np2dB);
end