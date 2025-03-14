function metrics = get_metrics_homo_gt(img, mask_homo, gt_homo, method)
% function metrics = get_metrics_homo_gt(img, mask_homo, gt_homo, method)

    reg_homo = img(mask_homo);

    metrics.method        = method;

    % Mean and std
    metrics.mean_homo    = mean(reg_homo);
    metrics.std_homo     = std(reg_homo);
    metrics.cv_homo      = 100*abs(metrics.std_homo / metrics.mean_homo);
    
    % MPE (bias) 
    err_homo = (reg_homo - gt_homo)./gt_homo;
    
    metrics.mpe_homo      = mean(err_homo);
    metrics.sdpe_homo     = std(err_homo);
    
    % MAE (Mean Absolute Error)
    abs_err_homo = abs(reg_homo - gt_homo)./abs(gt_homo);
    
    metrics.mae_homo      = mean(abs_err_homo);
    
    % RMSE
    err_homo_sq = (reg_homo - gt_homo).^2;
    metrics.rmse_homo     = sqrt(mean(err_homo_sq, 'omitnan'));  % rmse(reg_homo, gt_homo)
    % metrics.rmse_homo     = rmse(reg_homo, gt_homo);
    
    % NRMSE (Normalized RMSE) 
    % metrics.nrmse_homo    = metrics.rmse_homo / mean(gt_homo); % if GT is the same for all
    % RRMSE = sqrt(mean((y_predicted - y_actual).^2) )/ sqrt(mean(gt_homo.^2)) * 100; % In percentage
    metrics.nrmse_homo    = sqrt(mean( ( (reg_homo - gt_homo)./ gt_homo  ).^2) );


end