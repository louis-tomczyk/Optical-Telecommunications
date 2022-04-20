function S_t = dispersion_and_kerr(E_t,dl)
  
  S_t = dispersion(E_t,dl/4);
  S_t = kerr_only(S_t,dl/2);
  S_t = dispersion(E_t,dl/4);
  
endfunction