function beta2 = D2beta2(D,lambda)
  
  c       = 299792458;
  omega0  = 2*pi*c/lambda;  
  beta2   = -2*pi*c*D/omega0^2;
  
endfunction