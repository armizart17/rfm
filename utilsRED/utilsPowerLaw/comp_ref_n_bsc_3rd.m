function comp_factor = comp_ref_n_bsc_3rd(delta_n, band, p, q)
% function comp_factor = comp_ref_n_bsc(delta_n, band, p, q)
    band = band(:); % Ensure it's a column vector [r,1]
    comp_factor = repmat(band, [1, p, q]).^(-delta_n); % Size: [r, p, q]
    comp_factor = permute(comp_factor, [2, 3, 1]);    % Size: [p, q, r]

    % equivalent
    % comp_factor = reshape(band.^(-delta_n), [1, 1, numel(band)]) .* ones(p, q, 1);


end