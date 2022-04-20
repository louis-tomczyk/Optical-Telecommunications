% ================================
%   IDENTIFICATION
% ================================
% author --- louis tomczyk
% institution --- Telecom Paris
% date --- 2022/04/05
% contact --- louis.tomczyk.work@gmail.com
% ================================
%   CODE
% ================================
% Name of file --- FT.m
% Inputs ---  x : the function we want the Fourier Transform
%             t : the time array
% Outputs --- xfft : the Fourier Transform with the right normalization coefficient
% ================================
%   BIBLIOGRAPHY
% ================================

% ================================
%   REMARKS
% ================================
% This FFT is in agreement with Parseval-Plancherel Theorem stating that
% the power contained in the signal, in time or frequency domain, is the same.

% --- let introduce the function t->x(t) we note it "x_t".
% --- the time is noted "t"
% --- its Fourier transform is given by x(f) = F[t->x(t)](f) we note it "x_f"
% x_t  = ...  
% x_f  = FT(x_t,t)

% --- the power in both domains is given by Riemann Integrals (rectangles) using
% --- the L^2(R) norm, where "R" stands for all the real numbers from -infinity
% --- +infinity
% x_t = (abs(x_t)).^2;
% x_f = (abs(x_f)).^2;

% --- finally we sum all the retangles 
% power_t = sum(x_t)*dT
% power_f = sum(x_f)*dF

function xfft = FT(x,t)
  tmax      = t(end);
  Nsamples  = length(x);
  index     = 0:Nsamples-1;
  norm      = (-1).^index/Nsamples*tmax;
  xfft      = fftshift(norm.*fft(x));
endfunction