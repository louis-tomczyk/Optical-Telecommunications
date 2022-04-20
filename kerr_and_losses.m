function S_t = kerr_and_losses(E_t)

    global gamma
    global dl
    global alpha_lin
    global Nsamples

    z_eff       = (1-exp(-alpha_lin*dl))/alpha_lin;
    Phase_NL    = gamma*(abs(E_t)).^2.*z_eff;

    S_t         = E_t.*exp(-1i*Phase_NL).*exp(-alpha_lin*dl/2);;

    p = Phase_NL(floor(Nsamples)/2);
endfunction
