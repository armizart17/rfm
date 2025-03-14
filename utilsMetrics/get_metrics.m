function metrics = get_metrics(img, mask_inc, mask_back, method)
% function metrics = get_metrics(img, mask_inc, mask_back, method)

    % Mean and std
    metrics.mean_inc    = mean(img(mask_inc));
    metrics.mean_back   = mean(img(mask_back));
    metrics.std_inc     = std(img(mask_inc));
    metrics.std_back    = std(img(mask_back));
    
    % CNR
    cnr_eps = 1E-5; % to avoid dividing by 0
    metrics.cnr         = abs(metrics.mean_inc - metrics.mean_back) /...
                         ( cnr_eps + sqrt(metrics.std_inc.^2 + metrics.std_back.^2) );
    % Method name
    metrics.method      = method;

end