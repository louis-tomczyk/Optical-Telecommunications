clear all
close all
clc

alpha_dB  = 0.2/1000;% [dB/m]
alpha_lin = log(10)/10*alpha_dB;% [m⁻¹]
P0        = 1e0;% [W]
gamma     = 1.3/1000;% [W⁻¹.m⁻¹]
delta_phi = linspace(pi/6,pi/2,3);% [rad]
delta_z   = 0;linspace(0,1e4,10);% [m]
%%  
%%for k=1:length(delta_phi)
%%  dz2(k,:)  = -1/alpha_lin*log(1-alpha_lin*delta_phi(k)/gamma/P0*exp(alpha_lin*delta_z));
%%endfor
%% 
%%markers = cellstr(["o";"s";"h";"^";"d"]); 
%%figure
%%hold on
%%for k=1:length(delta_phi)
%%  plot(delta_z/1000,dz2(k,:)/1000,'color',[1,1,1]*exp(-0.4*k),"linewidth",20*exp(-0.5*k),"marker",markers{k})
%%endfor
%%xlabel("delta z [km]","fontsize",20,"fontweight","bold")
%%ylabel("dz2 [km]","fontsize",20,"fontweight","bold")
%%set(gca,"fontsize",15)
%%

dphi  = linspace(pi/12,pi/2,1e2);
L     = (floor(logspace(3,7,4)/1e3))*1e3;
leg_items = {};
for k=1:length(L)
  leg_items{k} = strcat(num2str(L(k)/1e3),' [km]');
endfor
markers = cellstr(["o";"s";"h";"d"]); 

for j=1:length(L)
  for k=1:length(dphi)
    N(j,k)    = round(gamma*P0/(alpha_lin*dphi(k))*(1-exp(-alpha_lin*L(j))));
  endfor
endfor

figure
plot(dphi,N,markers,"markersize",10)
leg = legend(leg_items)
xlabel("Δϕ [rad]","fontsize",20,"fontweight","bold")
ylabel("Number of points","fontsize",20,"fontweight","bold")
set(gca,"fontsize",15)
set(leg,"fontsize",15,'fontweight','bold')
title({"Number of points required","to compute the NLSE with constant phase shift"},"fontsize",20)
%%
%%n1 = 1:N(1);
%%n2 = 1:N(2);
%%n3 = 1:N(3);
%%
%%for k=1:length(dphi)
%%  zn1(k,:) = -1/alpha_lin*log(1-n1*alpha_lin*dphi(k)/(gamma*P0));
%%  zn2(k,:) = -1/alpha_lin*log(1-n2*alpha_lin*dphi(k)/(gamma*P0));
%%  zn3(k,:) = -1/alpha_lin*log(1-n3*alpha_lin*dphi(k)/(gamma*P0));
%%endfor
%%
%%figure
%%hold on
%%plot(n1,zn1(1,:)/1000)
%%plot(n1,zn1(2,:)/1000)
%%plot(n1,zn1(3,:)/1000)


