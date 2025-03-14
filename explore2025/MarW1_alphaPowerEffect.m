%% SPECTRAL PARAMETERS
pars.P = 4096; % NFFT only for calculate BSC_RPM_ok 10wl
% pars.P = 8192; % 15wl
pars.bw          = [3 8.5]; % [MHz] % old
pars.bw          = [3 9]; % [MHz] % I think better for BSC performance
pars.bw          = [3 8.75]; % [MHz] % new
pars.overlap     = 0.8;
pars.blocksize   = 12; % wavelengths

% new
pars.z_roi       = [5 45]*1E-3; % all
pars.x_roi       = [-18 18]*1E-3;
pars.window_type = 3; %  (1) Hanning, (2) Tuckey, (3) Hamming, (4) Tchebychev
pars.saran_layer = false;

%% GENERAL REGULARIZTATION SETTINGS
% Implementation parameters
par_rpl.tol        = 1e-16;
par_rpl.kmax       = 100;
par_rpl.eps_f      = 1e-16;
par_rpl.m_est      = 0; %Robust
par_rpl.ini_tol    = 1e-16;
par_rpl.df_op      = 1;
par_rpl.ini_method = 1; % METHOD LEAST SQUARES INITIALIZATION 

% Parameters for RPL-TV
mu_rpl_tv    = [1E3; 1E3; 1E4]; % [mu_g, mu_s, mu_a]

%% UTILS
Np2dB       = 20*log10(exp(1));
dB2Np       = 1/Np2dB;
%%
calc2dStats = {@(x) mean(x(:)), @(x) std(x(:)), @(x) 100 * std(x(:)) / mean(x(:))};
methods = {'3-DoF', '2-DoF-a', '2-DoF-g', '2-DoF-s'};

pathData = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\NonUniformBSC2D\';

% j_sam_values = [1.1, 1.2, 1.3, 1.4, 1.5];
% j_ref_values = [1.1, 1.2, 1.3, 1.4, 1.5];

j_sam_values = [1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9];
j_ref_values = [1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9];

alpha_sam = 0.5;
alpha_ref = 0.5;

% Preallocate results
num_sam = length(j_sam_values);
num_ref = length(j_ref_values);

maps_results_all = cell(num_sam, num_ref, length(methods));
bsc_results_all  = cell(num_sam, num_ref, length(methods));

tic
% Loop over all j_sam and j_ref combinations
for iSam = 1:num_sam
    for iRef = 1:num_ref
        j_sam = j_sam_values(iSam);
        j_ref = j_ref_values(iRef);

        % Update folder and file names
        folderDataSam = sprintf('');
        % folderDataSam = strrep(folderDataSam, '.', 'p');
        rf_sam_name = strcat('rf9_', sprintf('%.2g', j_sam));
        rf_sam_name = strrep(rf_sam_name, '.', 'p');
        rf_sam_name = strcat(rf_sam_name, '.mat');
        SAM = load(fullfile(pathData, folderDataSam, rf_sam_name));

        folderDataRef = sprintf('');
        % folderDataRef = strrep(folderDataRef, '.', 'p');
        rf_ref_name = strcat('rf7_', sprintf('%.2g', j_ref));
        rf_ref_name = strrep(rf_ref_name, '.', 'p');
        rf_ref_name = strcat(rf_ref_name, '.mat');
        REF = load(fullfile(pathData, folderDataRef, rf_ref_name));

        % Set values
        SAM.alpha_power = j_sam;
        SAM.acs = alpha_sam; % [dB/cm/MHz] 
        REF.alpha_power = j_ref; 
        REF.acs = alpha_ref; % [dB/cm/MHz]

        % Compute delta priors
        delta_alpha_prior   = alpha_sam - alpha_ref; 
        delta_g_prior       = log(1);
        delta_s_prior       = 0.0245;

        % Power spectra estimation
        spectralData_sam = calc_powerSpectra_vSimple(SAM, pars);
        S_sam = spectralData_sam.powerSpectra;
        spectralData_ref = calc_powerSpectra_vSimple(REF, pars);
        S_ref = spectralData_ref.powerSpectra;

        % Compute Spectral Ratio
        SR_emz = S_sam ./ S_ref;
        SR = permute(SR_emz, [3,1,2]); clear SR_emz

        band    = spectralData_sam.band;
        depth   = spectralData_sam.depth;
        [r,p,q] = size(SR);
        f = band(:); 

        comp_freq_a = comp_mod_freq_a(alpha_ref,j_sam,j_ref,band,depth,q);

        % Matrix
        X = kron( speye(p*q), ones(size(f)) );
        Z = kron( speye(p*q), -f.^2 );
        W = kron( speye(p*q), -4*f.^j_sam );
        
        dy = 0.5*(diag(ones(p-1,1),1) - diag(ones(p-1,1),-1));
        dy(1,1) = -1; dy(1,2) = 1; dy(end,end) = 1; dy(end,end-1) = -1;
        Dy = sparse(kron(speye(q),dy));
        z = 1E2*repmat(depth,1,q); % 1E2*spectralData_sam.depth * ones(1, q); % 2d array
        dz = reshape(Dy*z(:),p,q);
        dz(end,:) = dz(end-1,:);  

        % Loop over methods
        for iMet = 1:length(methods)
            estim_method = methods{iMet};

            if strcmp(estim_method, '2-DoF-a')
   
                comp_ref    = comp_ref_a(-delta_alpha_prior,j_ref,band,depth,q);                
                SR_comp     = SR .* comp_ref .* comp_freq_a;
                Y = log(SR_comp);

                u_0 = initialize_rpl_a_prior(Y, X, Z, mu_rpl_tv, par_rpl);
                [u_opt,~] = rpl_tv_a_prior(Y, X, Z, mu_rpl_tv, u_0, par_rpl);
                
                g = u_opt(1:p*q);
                s = u_opt(p*q+1:2*p*q);
                a_Np2dB = delta_alpha_prior*ones(p*q, 1);

            elseif strcmp( estim_method, '2-DoF-s')
                comp_ref    = comp_ref_s_bsc(delta_s_prior, band, p, q);
                SR_comp     = SR .* comp_ref .* comp_freq_a;
                Y = log(SR_comp);

                u_0 = initialize_rpl_n_prior(Y, X, W, mu_rpl_tv, par_rpl);
                [u_opt,~] = rpl_tv_n_prior(Y, X, W, mu_rpl_tv, u_0, par_rpl);

                g = u_opt(1:p*q);
                a = u_opt(p*q+1:2*p*q);
                s = delta_s_prior*ones(p*q, 1);
                a_Np2dB = Np2dB*Dy*a./dz(:);

            elseif strcmp( estim_method, '2-DoF-g')
                comp_ref    = comp_ref_g_bsc(delta_g_prior);
                SR_comp     = SR .* comp_ref .*comp_freq_a;
                Y = log(SR_comp);

                u_0 = initialize_rpl_b_prior(Y, Z, W, mu_rpl_tv, par_rpl);
                [u_opt,~] = rpl_tv_b_prior(Y, Z, W, mu_rpl_tv, u_0, par_rpl);
                
                s = u_opt(1:p*q);
                a = u_opt(p*q+1:2*p*q);
                g = delta_g_prior*ones(p*q, 1);
                a_Np2dB = Np2dB*Dy*a./dz(:);

            elseif strcmp(estim_method, '3-DoF')
                
                SR_comp = SR .* comp_freq_a;
                Y = log(SR_comp);

                u_0 = initialize_rpl(Y, X, Z, W, mu_rpl_tv, par_rpl);
                [u_opt,~] = rpl_tv(Y, X, Z, W, mu_rpl_tv, u_0, par_rpl);

                g = u_opt(1:p*q);
                s = u_opt(p*q+1:2*p*q);
                a = u_opt(2*p*q+1:3*p*q);
                a_Np2dB = Np2dB*Dy*a./dz(:);
            end

            % Compute final parameters
            g_ratio     = reshape(exp(g), p, q);
            g_ratio_dB  = 10*log10(g_ratio);
            alpha_ratio = reshape(a_Np2dB, p, q);
            s_ratio     = reshape(s, p, q); 
            acs_sam     = alpha_ratio + alpha_ref;

            % Compute statistics
            [m_a, s_a, cv_a] = deal(calc2dStats{1}(acs_sam), calc2dStats{2}(acs_sam), calc2dStats{3}(acs_sam));
            [m_g, s_g, cv_g] = deal(calc2dStats{1}(g_ratio_dB), calc2dStats{2}(g_ratio_dB), calc2dStats{3}(g_ratio_dB));
            [m_s, s_s, cv_s] = deal(calc2dStats{1}(s_ratio), calc2dStats{2}(s_ratio), calc2dStats{3}(s_ratio));

            % Save maps results
            maps_results_all{iSam, iRef, iMet} = {acs_sam, g_ratio_dB, s_ratio};

            % Compute and save BSC results
            g_est = median(g_ratio(:));
            s_est = median(s_ratio(:));
            bsc_est_gauss = g_est .* exp(-s_est .* band.^2);

            coeffs_pl = polyfit(log(spectralData_sam.band), log(bsc_est_gauss), 1);
            d_n_pl = coeffs_pl(1);
            d_b_pl = exp(coeffs_pl(2));
            bsc_fit_powlaw = d_b_pl * band.^d_n_pl;

            bsc_results_all{iSam, iRef, iMet} = {bsc_est_gauss, bsc_fit_powlaw};

            fprintf('=============== Method: %s ===============\n', estim_method);
            fprintf('j_sam=%.1f, j_ref=%.1f, α_s    : %.3f ± %.4f\n', j_sam, j_ref, m_a, s_a);
            fprintf('j_sam=%.1f, j_ref=%.1f, Δb [dB]: %.3f ± %.4f\n', j_sam, j_ref, m_g, s_g);
            fprintf('j_sam=%.1f, j_ref=%.1f, Δn     : %.3f ± %.4f\n', j_sam, j_ref, m_s, s_s);
   
        end
    end
end
tt = toc;
fprintf('Elapsed time %.4f\n', tt)
%% Save results

dirResults = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\qus-lim\out\alpha_power';
if (~exist(dirResults)) mkdir(dirResults); end
fileName = 'results_1p1_1p9';

save(fullfile(dirResults, fileName+".mat"), ...
    'bsc_results_all', 'maps_results_all', ...
    'par_rpl', 'j_sam_values', 'j_ref_values', 'alpha_sam', 'alpha_ref');
% save('bsc_maps_results.mat', 'bsc_results_all', 'maps_results_all');


%% CALCULATE GROUNDTRUTH RPM METHOD

Np2dB       = 20*log10(exp(1));
dB2Np       = 1/Np2dB;

calc2dStats = {@(x) mean(x(:)), @(x) std(x(:)), @(x) 100 * std(x(:)) / mean(x(:))};
methods = {'3-DoF', '2-DoF-a', '2-DoF-g', '2-DoF-s'};

pathData = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\NonUniformBSC2D\';

% j_sam_values = [1.1, 1.2, 1.3, 1.4, 1.5];
% j_ref_values = [1.1, 1.2, 1.3, 1.4, 1.5];

j_sam_values = [1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9];
j_ref_values = [1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9];

alpha_sam = 0.5;
alpha_ref = 0.5;

% Preallocate results
num_sam = length(j_sam_values);
num_ref = length(j_ref_values);

clear maps_results_all bsc_results_all

maps_results_all = cell(num_sam, num_ref);
bsc_results_all  = cell(num_sam, num_ref);

tic
% Loop over all j_sam and j_ref combinations
for iSam = 1:num_sam
    for iRef = 1:num_ref
        j_sam = j_sam_values(iSam);
        j_ref = j_ref_values(iRef);

        % Update folder and file names
        folderDataSam = sprintf('');
        % folderDataSam = strrep(folderDataSam, '.', 'p');
        rf_sam_name = strcat('rf9_', sprintf('%.2g', j_sam));
        rf_sam_name = strrep(rf_sam_name, '.', 'p');
        rf_sam_name = strcat(rf_sam_name, '.mat');
        SAM = load(fullfile(pathData, folderDataSam, rf_sam_name));

        folderDataRef = sprintf('');
        % folderDataRef = strrep(folderDataRef, '.', 'p');
        rf_ref_name = strcat('rf7_', sprintf('%.2g', j_ref));
        rf_ref_name = strrep(rf_ref_name, '.', 'p');
        rf_ref_name = strcat(rf_ref_name, '.mat');
        REF = load(fullfile(pathData, folderDataRef, rf_ref_name));

        % Set values
        SAM.alpha_power = j_sam;
        SAM.acs = alpha_sam; % [dB/cm/MHz] 
        REF.alpha_power = j_ref; 
        REF.acs = alpha_ref; % [dB/cm/MHz]

        BSC     = calculateBSC_RPM_fast(SAM, REF, pars); % fast TBD**
        bsc_rpm = BSC.BSCcurve_Uni(:,2); % median

        % Perform linear regression  bsc = d_s . f^2 + ln(d_g) 
        freq     = BSC.band;
        coeffs   = polyfit(-freq.^2, log(bsc_rpm), 1); % Fit y = mx + c
        d_s      = coeffs(1); % Slope = d_s -0.0319 (mean), -0.0317 (median)
        d_g      = exp(coeffs(2)); % Intercept = ln(d_g) 

        % Option 2: BSC = g.exp(-s.f^2) or b.f^n
        % Matrix [log(BSC)] = [ log(f) | 1] * [ n ; log(b) ] rr = MM * qq % POW LAW
        % Matrix [log(BSC)] = [ -f.^2  | 1] * [ s ; log(b) ] rr = MM * qq % GAUSS
        % MM           = [ log(freq), ones( size(freq) ) ]; % POW LAW
        % MM           = [ freq.^2,   ones( size(freq) ) ]; % GAUSS 
        % rr           = log(bsc);
        % qq           = cgs(MM' * MM, MM' * rr, 1e-16, 100); % qq = MM \ rr;
        % d_n          = qq(1);
        % d_b           = exp(qq(2));
       
        % Display results
        
        fprintf('======RPM Gauss (g.exp(-s.f^2))=====\n')
        fprintf('j_sam=%.1f, j_ref=%.1f\n', j_sam, j_ref);
        fprintf('d_g          = %f\n', d_g);
        fprintf('g_s/g_r [dB] = %f\n', 10*log10(d_g));
        fprintf('Δs           = %f\n', d_s);
        fprintf('===================================\n')
        
        bsc_rpm_gauss = d_g*exp(-d_s* freq.^2);

        bsc_results_all{iSam, iRef} = {bsc_rpm, bsc_rpm_gauss, 10*log10(d_g), d_s};
                      
    end
end
tt = toc;
fprintf('Elapsed time %.4f\n', tt)

%%
dirResults = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\qus-lim\out\alpha_power';
if (~exist(dirResults)) mkdir(dirResults); end
fileName = 'resultsRPM_1p1_1p9';

save(fullfile(dirResults, fileName+".mat"), ...
    'bsc_results_all', ...
    'par_rpl', 'j_sam_values', 'j_ref_values', 'alpha_sam', 'alpha_ref');
% save('bsc_maps_results.mat', 'bsc_results_all', 'maps_results_all');