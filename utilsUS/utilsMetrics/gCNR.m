function [gCNR_value] = gCNR(env1, env2, lim, samples)
% function [gCNR_value] = gCNR(env1, env2, lim, samples)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PYTHON CODE I. Salazar
% def gcnr(env1, env2):
%     _, bins = np.histogram(np.concatenate((env1, env2)), bins=128)
%     f, _ = np.histogram(env1, bins=bins, density=True)
%     g, _ = np.histogram(env1, bins=bins, density=True)
%     f /= f.sum()
%     g /= g.sum()
%     return 1 - np.sum(np.minimum(f, g))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin == 2 
        samples = 32;
        lim = [0 1];
        disp('Default bins = 256');
    end
    combined_env = [env1, env2];
    lim_v = linspace(lim(1), lim(2), samples);
    [~, bins] = histcounts(combined_env, lim_v);
%     [~, bins] = histcounts(combined_env, numbins);
    
    f = histcounts(env1, bins, 'Normalization', 'probability');
    g = histcounts(env2, bins, 'Normalization', 'probability');
    
    f = f / sum(f);
    g = g / sum(g);
    
    gCNR_value = 1 - sum(min(f, g));
end

