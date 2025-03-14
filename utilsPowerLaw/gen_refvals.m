function ref_vals = gen_refvals(band,data,test_params)

    src         = test_params.src;
    aref_num    = 4; %test_params.aref_num; %
    bscref_file = test_params.bscref_file;

    if src < 3
        ref_vals = data.ref_maps;
    elseif src < 6
        ref_vals = {1;0;0};
    elseif src == 3
        %ref_vals = {1;0;0};
        %ref_vals = {0.0042;1.692;0.6022};
        %ref_vals = {0.0023;2.0563;0.5131};
        ref_vals = {0.0023;2.0563;0.5169};
    else
        alpha_freqs       = attenuation_phantoms_Np(band,aref_num);
        [alpha_ref,~,~,~] = fit_linear(band,alpha_freqs',1);
        a_ref             = 8.686*alpha_ref;
        bsc_refdata       = load(bscref_file);
        bsc_faran         = bsc_refdata.bsc_faran;
        freqs             = bsc_refdata.freqs;
        bsc_ref           = interp1(freqs,bsc_faran,band)';
        p2                = polyfit(log(band),log(bsc_ref),1); 
        b_ref             = exp(p2(2));
        n_ref             = p2(1);
        ref_vals          = {b_ref;n_ref;a_ref};
    end

end