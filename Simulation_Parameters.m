%--------------------------------------------------------------------------
% Rates parameters
%--------------------------------------------------------------------------
global ReferenceFrequency
global SampleRate
global SymbolRate
global BitsPerSymbol
global Nsamples
global Nsymbols

% in [Hz]
ReferenceFrequency  = 193.41e12; 

% in [baud]
SymbolRate          = 25e9; 
BitsPerSymbol    = 512;
Nsymbols            = 128;

SampleRate          = BitsPerSymbol*SymbolRate;
Nsamples            = Nsymbols*BitsPerSymbol;

%--------------------------------------------------------------------------
% Time and Frequency Axis
%--------------------------------------------------------------------------
global TimeArray
global dT
global FrequencyArray
global dF

% Time interval between two consecutive samples in [s]
dT = 1/SampleRate;

% Frequency interval between two consecutive samples in [Hz]
dF = SymbolRate/Nsymbols;

TimeArray       = (0:(Nsamples-1))*dT;
FrequencyArray  = (-Nsamples/2:Nsamples/2-1)*dF;
TimeArray_Mid   = floor(Nsamples)/2*dT;

%--------------------------------------------------------------------------
% Physical parameters
%--------------------------------------------------------------------------
%--------------------------
% --- fundamental parameters
%--------------------------

global c
global omega0

c       = 299792458;              % [m/s]
lambda  = c/ReferenceFrequency;   % [m]
omega0  = ReferenceFrequency*2*pi;% [rad/s]

%--------------------------
% --- pulse parameters
%--------------------------
global T0

T0      = 1e-11;                % [s]
P0      = 1/(T0*sqrt(2*pi));    % [W]
t_FWHM  = 2*sqrt(log(2))*T0;    % [s]
nu_FWHM = sqrt(log(2))/(pi*T0); % [Hz]

%--------------------------
% --- fibre parameters
% CF ( ENSSAT, Peucheret ) Linear propagation in optical fibres
%--------------------------
global L
global Aeff
global alpha_lin
global beta2
global gamma

L = 2.5e4;

what_fibre = "HNLF"%input("\t type of fibre?\t",'s');

switch what_fibre
  case "SMF"
    alpha_dB  = 2e-4;                 % [m⁻¹]
    alpha_lin = log(10)/10*alpha_dB;  % [dB/m]
    
    D         = 17e-6;                % [s/m²]
    S         = 0.058e3;              % [s/m³]
    beta2     = D2beta2(D,lambda);    % [s²/m/rad²]
    LD        = T0^2/abs(beta2);      % [m]
    
    Aeff      = 80e-12;                     % [m²]
    gamma     = 1.3e-3;                     % [1/W/m]
    n2        = gamma2n2(gamma,lambda,Aeff);% [m²/W]
    
  case "DSF"
    alpha_dB  = 2e-4;                 % [m⁻¹]
    alpha_lin = log(10)/10*alpha_dB;  % [dB/m]
    
    D         = 0;                    % [s/m²]
    S         = 0.08e3;               % [s/m³]
    beta2     = D2beta2(D,lambda);    % [s²/m/rad²]
    
    Aeff      = 50e-12;                     % [m²]
    gamma     = 2.1e-3;                     % [1/W/m]
    n2        = gamma2n2(gamma,lambda,Aeff);% [m²/W]
    
  case "NZDSF"
    alpha_dB  = 2e-4;                 % [m⁻¹]
    alpha_lin = log(10)/10*alpha_dB;  % [dB/m]
    
    D         = 4.5e-6;               % [s/m²]
    S         = 0.045e3;              % [s/m³]
    beta2     = D2beta2(D,lambda);    % [s²/m/rad²]
    LD        = T0^2/abs(beta2);      % [m]
    
    Aeff      = 50e-12;                     % [m²]
    gamma     = 2.3e-3;                     % [1/W/m]
    n2        = gamma2n2(gamma,lambda,Aeff);% [m²/W]

  case "LEAF"
    alpha_dB  = 2e-4;                 % [m⁻¹]
    alpha_lin = log(10)/10*alpha_dB;  % [dB/m]
    
    D         = 4.0e-6;               % [s/m²]
    S         = 0.09e3;               % [s/m³]
    beta2     = D2beta2(D,lambda);    % [s²/m/rad²]
    LD        = T0^2/abs(beta2);      % [m]
    
    Aeff      = 72e-12;                     % [m²]
    gamma     = 1.5e-3;                     % [1/W/m]
    n2        = gamma2n2(gamma,lambda,Aeff);% [m²/W]

  case "HNLF"
    alpha_dB  = 15e-3;                 % [m⁻¹]
    alpha_lin = log(10)/10*alpha_dB;  % [dB/m]
    
    D         = -0.06e-6;               % [s/m²]
    S         = 0.011e3;               % [s/m³]
    beta2     = D2beta2(D,lambda);    % [s²/m/rad²]
    LD        = T0^2/abs(beta2);      % [m]
    
    Aeff      = 13.5e-12;                     % [m²]
    gamma     = 7.4e-3;                     % [1/W/m]
    n2        = gamma2n2(gamma,lambda,Aeff);% [m²/W]
endswitch




% typical non linearities distance in [m]
% the total power is calculated using Riemann Integration point of view
% i.e. small rectangles : width (dT)*height (sig_opt).
% the factor *2 comes from the fact that we centered the pulse
% TimeArray-TimeArray_Mid around 0...
%%sig_opt = (abs(sqrt(P0)*u)).^2;
%%total_power_in = sum(sig_opt)*dT*2
%%LNL = 1/(gamma*total_power_in);

%--------------------------
% --- pulse shape
%--------------------------

Type = "gaussian";%input("pulse shape type? gaussian or soliton? ","s");

if strcmp(Type,"gaussian")==1
  u   = exp(-(TimeArray-TimeArray_Mid).^2./T0^2);
  
elseif strcmp(Type,"soliton")==1
  P0  = abs(beta2)/(gamma*T0^2);
  u   = sqrt(P0)*1./cosh((TimeArray-TimeArray_Mid)./T0);
  
endif

%--------------------------
% --- Computation parameters 
%--------------------------

global dl
global N_steps

% length step in [m]
dl = 1000;

% total number of steps
N_steps = L/dl;