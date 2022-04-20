% ================================
%   IDENTIFICATION
% ================================
% author --- louis tomczyk
% institution --- Telecom Paris
% date --- 2022/04/07
% contact --- louis.tomczyk.work@gmail.com
% ================================
%   CODE
% ================================
% Name of file --- iFT.m
% Inputs ---  xfft : the function we want the Inverse Fourier Transform
%             SampleRate : also called Sampling Frequency
% Outputs --- x : the inverse Fourier Transform with the right normalization coefficient
% ================================
%   BIBLIOGRAPHY
% ================================

% ================================
%   REMARKS
% ================================
% This iFFT is in agreement with Parseval-Plancherel Theorem stating that
% the power contained in the signal, in time or frequency domain, is the same.

% --- let introduce the function t->x(t) we note it "x_t".
% --- the frequency is noted "t"
% --- its Fourier transform is given by x(f) = F[t->x(t)](f) we note it "x_f"
% --- its inverse Fourier transform is given by x(t) = F⁻¹[f->x(f)](t) we note it "x_t"

% x_t = ...
% x_f = FT(x_t)
% x_tt = iFT(x_f)

% if the iFT function works, then the plots should exactly overlap:
% figure
% hold on
% plot(TimeArray-TimeArray_Mid,tmp_t,'linewidth',10,'k')
% plot(TimeArray-TimeArray_Mid,tmp_tt)

function x_t = iFT(x_f,SampleRate)
  norm       = SampleRate;
  x_t        = norm.*fftshift(ifft(fftshift(x_f)));
endfunction