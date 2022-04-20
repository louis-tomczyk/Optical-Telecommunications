function S_t = NLSE(E_t,what,n_steps)
  
  global N_steps
  
  switch what
    case "dispersion"
      S_t = dispersion(E_t);
      for k=2:n_steps
        S_t = dispersion(S_t);
      endfor
      
    case "kerr only"
      S_t = kerr_only(E_t);
      for k=2:n_steps
        S_t = kerr_only(S_t);
      endfor
      
    case "kerr and losses"
      S_t = kerr_only(E_t);
      for k=2:n_steps
        S_t = kerr_and_losses(S_t);
      endfor
      
    case "dispersion and kerr"
      S_t = dispersion_and_kerr(E_t);
      for k=2:n_steps
        S_t = dispersion_and_kerr(S_t);
      endfor
  endswitch
endfunction