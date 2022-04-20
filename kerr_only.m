function S_t = kerr_only(E_t,dl)

    global Nsamples
    global gamma
    
    Phase_NL = gamma*(abs(E_t)).^2*dl;
    S_t = E_t.*exp(1i*Phase_NL);

endfunction