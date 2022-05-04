%   ex17: Optical propagation within an optical fiber impaired by only
%   group-velocity dispersion (GVD).

clear all
close all
clc

%% Global parameters
ex01_inigstate_bis
Nsymb = 16;             % number of symbols
Nt    = 32;             % number of discrete points per symbol
Np    = 50

%% Tx parameters
symbrate    = 10;           % symbol rate [Gbaud].
tx.rolloff  = 0.4;          % pulse roll-off
tx.emph     = 'asin';       % digital-premphasis type
modfor      = 'ook';        % modulation format
PdBm        = 0;            % power [dBm]
Plin        = 10.^(PdBm/10);% [mW]
lam         = 1550;         % carrier wavelength [nm]

%% Channel parameters

Dtot = 0;                   % residual dispersion per span [ps/nm]

%% Transmission fiber
Lengths       = linspace(1e3,1e5,Np);  % length [m]
ft.lambda     = 1550;                   % wavelength [nm] of fiber parameters
ft.alphadB    = 0;                      % attenuation [dB/km]

ft.disp        = 17;                     % dispersion [ps/nm/km] @ ft.lambda
ft.slope      = 0;                      % slope [ps/nm^2/km] @ ft.lambda
ft.n2         = 2.5656e-20;             % nonlinear index [m^2/W]
ft.aeff       = 80;                     % effective area [um^2]

ft.gamma      = n22gamma(ft.n2,lam*1e-9,ft.aeff*1e-12)

%% Rx parameters
rx.modformat = modfor;      % modulation format
rx.sync.type = 'da';        % time-recovery method
rx.oftype = 'gauss';        % optical filter type
rx.obw = Inf;               % optical filter bandwidth normalized to symbrate
rx.eftype = 'rootrc';       % optical filter type
rx.ebw = 0.5;               % electrical filter bandwidth normalized to symbrate
rx.epar = tx.rolloff;       % electrical filter extra parameters
rx.type = 'bin';            % binary pattern

%% Init
Nsamp   = Nsymb*Nt;         % overall number of samples
fs      = symbrate*Nt;      % sampling rate [GHz]
inigstate(Nsamp,fs);        % initialize global variables: Nsamp and fs.




Nt      = 5;
m       = 0;
Max_u   = 100;
eps     = 15;      % the higher it gets, the lower the surface of the shaded area gets
f       = csvread("student_law_coefficients.txt");

while Max_u>eps

  RMSEs           = zeros(Nt,Np);
  
  for j=1:Nt
  
    E               = lasersource(Plin,lam,struct('pol','single'));  % y-pol does not exist
    [patx, patbinx] = pattern(Nsymb,'rand',struct('format',modfor));
    [sigx, normx]   = digitalmod(patbinx,modfor,symbrate,'costails',tx);
    E               = iqmodulator(E, sigx,struct('norm',normx));
    E_in            = E;

    for k=1:Np
      ft.length   = Lengths(k);
      E_out       = fiber(E_in,ft);
      RMSEs(j,k)  = RMSE(E_in.field,E_out.field)*100;
    endfor
  endfor

  %% uncertainties
  RMSEs   = RMSEs/Plin;
  RMSEs_m = mean(RMSEs);
  t       = student_coefficient(Nt,95,f);
  RMSEs_u = t*sqrt(1/(Nt^2-Nt)*sum((RMSEs-RMSEs_m).^2));
  Max_u   = 2*mean(RMSEs_u)
  Nt      = Nt+5;
  m       = m+1;
endwhile

Lengthss = [Lengths,fliplr(Lengths)];
Y2l     = RMSEs_m-RMSEs_u;
Y2h     = RMSEs_m+RMSEs_u;
Max_Y2d = max(Y2h-Y2l);
B       = [Y2l,fliplr(Y2h)];

% fit
k         = 0;
poly_rmse = 100;
eps       = 5; % the lower it gets, the more precise is the fitting

%%while poly_rmse>eps
%%  
%%  p         = polyfit(log(Lengths),RMSEs_m,1);
%%  xp        = linspace(Lengths(1),Lengths(end),Np);
%%  yp        = polyval(p,xp);
%%  poly_rmse = RMSE(RMSEs_m/100,yp);
%%  k         = k+1;
%%  
%%endwhile

A = RMSEs_m(1);
B = RMSEs_m(end);
L0 = Lengths(1);
Lend = Lengths(end);

alpha = (A-B)/log(L0/Lend);
beta = A-alpha*log(L0);

xp = log(Lengths)
yp = alpha*log(Lengths)+beta;

ts  = 1/GSTATE.FSAMPLING; % sampling time [ns]
tim = linspace(0,(GSTATE.NSAMP-1)*ts,length(E_in.field));

%% plot
figure
  subplot(1,2,1)
    hold on
    fill(Lengthss/1e3,10*log10(B),[0.2,0.2,0.2],"facealpha",0.1)
    plot(Lengths/1e3,10*log10(RMSEs_m),'k','linewidth',5)
    plot(exp(xp)/1e3,yp,'b--','linewidth',3)
    xlabel("fiber length [km]","fontsize",20,"fontweight","bold")
    ylabel("10.log10(RMSE) [dB]","fontsize",20,"fontweight","bold")
    axis([Lengths(1)/1e3,Lengths(end)/1e3])
    title("Evolution of the RMSE/P_{in} ratio with the dispersion")
    text(2,30,sprintf("fiber length = %.1f [km]",ft.length/1e3),"fontsize",15)
    text(2,28,sprintf("baud rate = %.1f [Gbaud]",symbrate),"fontsize",15)
    text(2,26,sprintf("attenuation = %.1f [dB/km]",ft.alphadB),"fontsize",15)
    text(2,24,sprintf("non linearity = %.1f [1/W/km]",ft.gamma*1e3),"fontsize",15)
    set(gca,"fontsize",20)
  subplot(1,2,2)
    hold on
    plot(tim,empower(E_in.field),'k','linewidth',3)
    plot(tim,empower(E_out.field),'color',[0.8,0.8,0.8],'linewidth',3)
    xlabel("time [ns]","fontsize",20,"fontweight","bold")
    ylabel("power [mW]","fontsize",20,"fontweight","bold")
    title("Time profile","fontsize",20)
    legend("TX","RX")
    set(gca,"fontsize",20)
