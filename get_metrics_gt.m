function metrics = get_metrics_gt(img, mask_inc, mask_back, method, gt_inc, gt_back)
% function metrics = get_metrics_gt(img, mask_inc, mask_back, method, gt_inc, gt_back)
    
    % Regions vectorized
    reg_inc = img(mask_inc); 
    reg_back = img(mask_back);

    % Mean and std
    metrics.mean_inc    = mean(reg_inc);
    metrics.mean_back   = mean(reg_back);
    metrics.std_inc     = std(reg_inc);
    metrics.std_back    = std(reg_back);
    
    % CNR
    cnr_eps = 1E-5; % to avoid dividing by 0
    metrics.cnr         = abs(metrics.mean_inc - metrics.mean_back) /...
                         ( cnr_eps + sqrt(metrics.std_inc.^2 + metrics.std_back.^2) );

    % MPE (bias) 
    err_inc = (reg_inc - gt_inc)/gt_inc;
    err_back = (reg_back - gt_back)/gt_back;

    metrics.mpe_inc      = mean(err_inc);
    metrics.sdpe_inc     = std(err_inc);

    metrics.mpe_back      = mean(err_back);
    metrics.sdpe_back     = std(err_back);

    % MAE (Mean Absolute Error)
    abs_err_inc = abs(reg_inc - gt_inc)/gt_inc;
    abs_err_back = abs(reg_inc - gt_back)/gt_back;

    metrics.mae_inc      = mean(abs_err_inc);

    metrics.mae_back     = mean(abs_err_back);

    % RMSE
    err_inc_sq = (reg_inc - gt_inc).^2;
    err_back_sq = (reg_back - gt_back ).^2;

    metrics.rmse_inc     = sqrt(mean(err_inc_sq, 'omitnan'));  % rmse(reg_inc, gt_inc)
    metrics.rmse_back    = sqrt(mean(err_back_sq, 'omitnan')); % rmse(reg_back, gt_back)

    % NRMSE (Normalized RMSE)
    metrics.nrmse_inc     = metrics.rmse_inc / mean(gt_inc);
    metrics.nrmse_back    = metrics.rmse_back / mean(gt_back);

    % Method name
    metrics.method        = method;

end