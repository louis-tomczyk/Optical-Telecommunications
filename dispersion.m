function S_t = dispersion(E_t,dl)
  
    global beta2
    global FrequencyArray
    global TimeArray
    global SampleRate
    
    H = exp(1i*beta2*(2*pi*FrequencyArray).^2/2*dl);
    
    E_f = FT(E_t,TimeArray);
    S_f = H.*E_f;
    S_t = iFT(S_f,SampleRate);
    
endfunction
