function comp_factor = comp_ref_n_bsc(delta_n, band, p, q)
% function comp_factor = comp_ref_n_bsc(delta_n, band, p, q)
    % comp_factor = permute(repmat(band,p,1,q),[2,1,3]).^(-delta_n);
    band = band(:); % make sur it's col vector
    comp_factor = repmat(band,[1,p,q]).^(-delta_n);
end