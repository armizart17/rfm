function comp_factor = comp_ref_s_bsc(delta_s,band,p,q)
    comp_factor = exp(+delta_s*repmat(band.^2,[1,p,q]));
end