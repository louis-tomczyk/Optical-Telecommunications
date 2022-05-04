%   ex11: evaluate the eye diagram of a phase-modulated signal and the
%   corresponding best eye-opening.

clear all
close all
clc

%% Global parameters

Nsymb = 1024;               % number of symbols
Nt    = 32;                 % number of discrete points per symbol

%% Tx parameters
symbrate  = 10;             % symbol rate [Gbaud]
R         = linspace(0,1,1)
tx.emph   = 'asin';         % digital-premphasis type
modfor    = 'qpsk';         % modulation format
PdBm      = 0;              % power [dBm]       
lam       = 1550;           % carrier wavelength [nm]  

%% Rx parameters
rx.modformat  = modfor;   % modulation format
rx.sync.type  = 'da';     % time-recovery method
rx.oftype     = 'gauss';  % optical filter type
rx.obw        = Inf;      % optical filter bandwidth normalized to symbrate
rx.eftype     = 'rootrc'; % optical filter type
rx.ebw        = 0.5;      % electrical filter bandwidth normalized to symbrate

%% Init global variables
Nsamp = Nsymb*Nt;        % overall number of samples
fs    = symbrate*Nt;     % sampling rate [GHz]
inigstate(Nsamp,fs);     % initialize global variables: Nsamp and fs.

%% Tx side

Plin    = 10.^(PdBm/10);          % [mW]
E       = lasersource(Plin,lam);  % electric field
[Ex,Ey] = pbs(E);                 % split in two orthogonal polarizations

patx    = pattern(Nsymb,'rand',struct('format',modfor));
paty    = pattern(Nsymb,'rand',struct('format',modfor));

tx.rolloff      = R;      % pulse roll-off
[elecx, normx]  = digitalmod(patx,modfor,symbrate,'rootrc',tx);
[elecy, normy]  = digitalmod(paty,modfor,symbrate,'rootrc',tx);

Ex  = iqmodulator(Ex, elecx,struct('norm',normx));
Ey  = iqmodulator(Ey, elecy,struct('norm',normy));
E   = pbc(Ex,Ey); % combine creating a PDM signal

%% Rx side

rx.epar = R;    % electrical filter extra parameters
rsig    = rxfrontend(E,lam,symbrate,rx);    % front-end
eyeop   = evaleye([patx paty],angle(rsig),symbrate,modfor,struct('plot',true));
fprintf('Eye opening: %.2f [dB]\n',eyeop)
