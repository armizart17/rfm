function cf = dfc_factor(agl,band,depth,lateral,data)

    ref_t    = data.ref_t;
    f_length = data.f_length;
    rad_t    = data.rad_t;
    c0       = data.c0;

    p = length(depth);
    q = length(lateral);
    r = length(band);

    z = repmat(depth,[1,r])';

    g_p  = (pi*rad_t^2)./(f_length*c0./(band*1e+6));
    G_p  = repmat(g_p,[1,p]);
    cond = (1./(1+pi./G_p) <= z./f_length).*(z./f_length <= 1./(1-pi./G_p));
    u_f  = G_p.*(f_length./z-1);

    D_z   = cond.*(0.46*pi*rad_t^2).*exp((-0.46/pi)*(u_f.^2))./(z.^2) + (1-cond).*(1.07*pi*rad_t^2).*(u_f.^-2)./(z.^2);
    D_ref = abs(1-exp(-1j*G_p).*(besselj(0,G_p) + 1j*besselj(1,G_p))).^2;

    cf = repmat((D_ref./D_z)/(agl/ref_t^2),[1,1,q]);

end