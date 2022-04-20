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
% Name of file ---
% Inputs ---
% Outputs ---
% ================================
%   BIBLIOGRAPHY
% ================================

% --------------------------------
%   Maintenance
% --------------------------------

%addpath /home/louis/Documents/6___Telecom_Paris/Codes/optilux_copy

clear all
close all
clc

% --------------------------------
%   Simulation parameters
% --------------------------------
Simulation_Parameters

%%% pulse in time domain
%%figure
%%plot(TimeArray-TimeArray_Mid,sig_opt_in_t)
%%axis([-2*t_FWHM,2*t_FWHM])
%%
%%% pulse in frequency domain
%%figure
%%plot(FrequencyArray,sig_opt_in_f)
%%axis([-2*nu_FWHM,2*nu_FWHM])

% --------------------------------
%   Signals
% --------------------------------

%%sig_in_t = u;
%%sig_in_f = FT(sig_in_t,TimeArray);
%%
%%sig_in_f_power = (abs(sig_in_f)).^2;
%%sig_in_t_power = (abs(sig_in_t)).^2;

% --- verifying FT
%sig_t = u;
%sig_t = (abs(u)).^2;
%sig_f = (abs(FT(u,TimeArray))).^2;
%
%power_t = sum(sig_t)*dT;
%power_f = sum(sig_f)*dF;
%
%assert(abs(power_t-power_f)<1e-2)

% --- verifying IFT
%sig_t = u;
%sig_f = FT(sig_t,TimeArray);
%sig_tt = iFT(sig_f,SampleRate);
%
%
%sig_t = (abs(sig_t)).^2;
%sig_tt = (abs(sig_tt)).^2;
%
%assert(abs(sig_t-sig_tt)<1e-2)


% --------------------------------
%   Dispersion only
% --------------------------------
%%% --- signals in
%%sig_in_t = u;
%%sig_in_f = FT(sig_in_t,TimeArray);
%%
%%sig_in_t_power = (abs(sig_in_t)).^2;
%%sig_in_f_power = (abs(sig_in_f)).^2;
%%
%%% --- dispersion
%%for k=1:N_steps
%%  sig_out_t = dispersion(sig_in_t,dl);
%%endfor
%%
%%% --- signals out
%%sig_out_f = FT(sig_out_t,TimeArray);
%%sig_out_f_power = empower(sig_out_f);
%%sig_out_t_power = empower(sig_out_t);
%%
%%% --- pulse in time domain
%%figure
%%hold on
%%scatter(TimeArray-TimeArray_Mid,sig_in_t_power,'s')
%%scatter(TimeArray-TimeArray_Mid,sig_out_t_power)
%%axis([-10*t_FWHM,10*t_FWHM])
%%
%%% --- pulse in frequency domain
%%figure
%%hold on
%%scatter(FrequencyArray,sig_in_f_power,100,'s')
%%scatter(FrequencyArray,sig_out_f_power)
%%axis([-2*nu_FWHM,2*nu_FWHM])


% --------------------------------
%   Kerr only
% --------------------------------

%%SOTP = zeros(N_steps,Nsamples);
%%SOFP = zeros(N_steps,Nsamples);
%%Phase_NL = zeros(1,Nsamples);
%%
%%% --- kerr
%%sig_out_t = kerr_and_losses(sig_in_t);
%%sig_out_f = FT(sig_out_t,TimeArray);
%%
%%SOTP(1,:) = empower(sig_out_t);
%%SOFP(1,:) = empower(sig_out_f);
%%
%%for k=2:N_steps
%%%  sig_out_t = kerr_and_losses(sig_out_t);
%%  sig_out_t = kerr_only(sig_out_t,dl);
%%  sig_out_f = FT(sig_out_t,TimeArray);
%%  SOFP(k,:) = empower(sig_out_f);
%%  SOTP(k,:) = empower(sig_out_t);
%%endfor
%%
%%N_steps_star = round(3.5*pi*LNL/dl)
%%
%%l = 0:dl:L;
%%z_eff       = (1-exp(-alpha_lin*l))/alpha_lin;
%%P0          = sig_in_t_power(floor(Nsamples/2));
%%PL          = sig_in_t_power(end);
%%%Phase_NL    = gamma*P0.*z_eff;
%%Phase_NL    = gamma*P0*l;
%%
%%% --- pulse in time domain
%%figure
%%pause(2)
%%for k=1:N_steps_star
%%  subplot(1,2,1)
%%    plot((TimeArray-TimeArray_Mid)*1e12,SOTP(k,:),'color','black','linewidth',5)
%%    axis([-10*t_FWHM*1e12,10*t_FWHM*1e12,0,1])
%%    xlabel("time [ps]",'fontsize',20)
%%    ylabel("power [a.u.]",'fontsize',20)
%%    title("Time Domain",'fontsize',20)
%%    text(2*t_FWHM*1e12,0.9,sprintf("distance = %i [km]",k*dl/1000),'fontsize',20)
%%    
%%  subplot(1,2,2)
%%    plot(FrequencyArray*1e9,SOFP(k,:),'color','black','linewidth',5)
%%    axis([-20*nu_FWHM*1e9,20*nu_FWHM*1e9,0,1.9*10^(-22)])
%%    xlabel("frequency [GHz]",'fontsize',20)
%%    ylabel("power [a.u.]",'fontsize',20)
%%    title("Frequency Domain",'fontsize',20)
%%    text(2*nu_FWHM*1e9,1.5e-22,sprintf("accumulated phase = %.3f pi",Phase_NL(k+1)/pi),'fontsize',20)
%%
%%  pause(0.5)
%%endfor

% --------------------------------
%   Kerr and Dispersion
% --------------------------------
SOT = zeros(N_steps,Nsamples);
SOF = zeros(N_steps,Nsamples);
SOTP = zeros(N_steps,Nsamples);
SOFP = zeros(N_steps,Nsamples);

sig_in_t = u;
sig_in_f = FT(sig_in_t,TimeArray);

SOT(1,:) = sig_in_t;
SOF(1,:) = sig_in_f;
SOTP(1,:) = (abs(sig_in_t)).^2;
SOFP(1,:) = (abs(sig_in_f)).^2;

sig_out_t = dispersion_and_kerr(sig_in_t,dl);
sig_out_f = FT(sig_out_t,TimeArray);

SOT(2,:) = sig_in_t;
SOF(2,:) = sig_in_f;
SOTP(2,:) = (abs(sig_in_t)).^2;
SOFP(2,:) = (abs(sig_in_f)).^2;

figure
for k=3:N_steps
%  sig_out_t = kerr_only(sig_in_t,k*dl);
  sig_out_t = dispersion_and_kerr(sig_out_t,dl);
  sig_out_f = FT(sig_out_t,TimeArray);
  
  SOF(k,:) = sig_out_f;
  SOT(k,:) = sig_out_t;
  
  SOFP(k,:) = empower(sig_out_f);
  SOTP(k,:) = empower(sig_out_t);
endfor


% --- pulse in time domain


for k=1:N_steps
  subplot(1,2,1)
    plot((TimeArray-TimeArray_Mid)*1e12,SOTP(k,:),'color','black','linewidth',5)
    axis([-10*t_FWHM*1e12,10*t_FWHM*1e12,0,1])  %--- power
%      axis([-10*t_FWHM*1e12,10*t_FWHM*1e12,0,1])% --- field
    xlabel("time [ps]",'fontsize',20)
    ylabel("power [a.u.]",'fontsize',20)
    title("Time Domain",'fontsize',20)
    text(2*t_FWHM*1e12,0.9,sprintf("distance = %i [km]",k*dl/1000),'fontsize',20)
    
  subplot(1,2,2)
    plot(FrequencyArray*1e9,SOFP(k,:),'color','black','linewidth',5)
      axis([-10*nu_FWHM*1e9,10*nu_FWHM*1e9,0,3e-22])% --- power
%      axis([-20*nu_FWHM*1e9,20*nu_FWHM*1e9,0,1])% --- field
    xlabel("frequency [GHz]",'fontsize',20)
    ylabel("power [a.u.]",'fontsize',20)
    title("Frequency Domain",'fontsize',20)
%    text(2*nu_FWHM*1e9,1.5e-22,sprintf("accumulated phase = %.3f pi",Phase_NL(k+1)/pi),'fontsize',20)

  pause(1e-3)
endfor
