%% Informations
% % ================================
% %   IDENTIFICATION
% % ================================
% % author --- louis tomczyk
% % institution --- Telecom Paris
% % date --- 2022-05-10
% % contact --- louis.tomczyk@telecom-paris.com
% % ================================
% %   CODE
% % ================================
% % Name of file ---
% % Inputs ---
% % Outputs ---
% % ================================
% %   BIBLIOGRAPHY
% % ================================
% 
% %%   maintenance
clear
close all
clc
clc
clc

folder = fileparts(which(mfilename));
folder = folder(1:end-11);
addpath(genpath(folder));

%% parameters

% Parameters - loop
Nt  = 1e3;
Ntt = 1e0;

% Parameters - fiber
L   = 1e3; % [m]
Np  = 20;
PMD = 10*10^(-13.5); % [s/sqrt(m)]
DGD = zeros(Ntt,Nt); % [s]

% Parameters - random shifts
alphaRand   = 1;
phiRand     = 1;
dlRand      = 1;

if alphaRand == 0
    a_min = 0;
    a_max = 0;
else
    a_min = 0;
    a_max = 2*pi;
end


if phiRand == 0
    p_min = 0;
    p_max = 0;
else
    p_min = 0;
    p_max = 2*pi;
end

dL_mean = L/Np;


a_var   = a_max-a_min;
p_var   = p_max-p_min;
dL_var  = dL_mean/5*dlRand;

% Parameters - spectral properties
lambda  = 1550e-9;
nu      = lambda2nu(lambda);
omega   = 2*pi*nu;
dnu     = 12.5e9;
domega  = 2*pi*dnu;

%% PMD calculations

parfor i = 1:Ntt
    for j = 1:Nt
        for k=1:Np
          if k ==1
            dL      = dL_mean+dL_var*randn; % [m]
            a       = a_min+a_var*rand; % [rad]
            p       = p_min+p_var*rand; % [rad]

            B       = diag([exp(1i*PMD*omega*sqrt(dL)/2+1i*p),...
                            exp(-1i*PMD*omega*sqrt(dL)/2-1i*p)],0);
            Bm      = diag([exp(1i*PMD*(omega-domega)*sqrt(dL)/2+1i*p),...
                            exp(-1i*PMD*(omega-domega)*sqrt(dL)/2-1i*p)],0);
            BM      = diag([exp(1i*PMD*(omega+domega)*sqrt(dL)/2+1i*p),...
                            exp(-1i*PMD*(omega+domega)*sqrt(dL)/2+1i*p)],0);    
            R       = [cos(a),sin(a);...
                       -sin(a),cos(a)];

            T       = R*B;
            Tm      = R*Bm;
            TM      = R*BM;
              
            tmp     = T;
            tmp_m   = Tm;
            tmp_M   = TM;
    
          else
            dL      = dL_mean+dL_var*randn; % [m]
            a       = a_min+a_var*rand; % [rad]
            p       = p_min+p_var*rand; % [rad]

            B       = diag([exp(1i*PMD*omega*sqrt(dL)/2+1i*p),...
                            exp(-1i*PMD*omega*sqrt(dL)/2-1i*p)],0);
            Bm      = diag([exp(1i*PMD*(omega-domega)*sqrt(dL)/2+1i*p),...
                            exp(-1i*PMD*(omega-domega)*sqrt(dL)/2-1i*p)],0);
            BM      = diag([exp(1i*PMD*(omega+domega)*sqrt(dL)/2+1i*p),...
                            exp(-1i*PMD*(omega+domega)*sqrt(dL)/2+1i*p)],0);    
            R       = [cos(a),sin(a);...
                       -sin(a),cos(a)];

            T       = R*B;
            Tm      = R*Bm;
            TM      = R*BM;
              
            tmp     = tmp*T;
            tmp_m   = tmp_m*Tm;
            tmp_M   = tmp_M*TM;
          end
        end
        
        dTdomega    = abs((tmp_M-tmp_m)/(2*domega));
        Tinv        = T^(-1);
        U           = dTdomega*Tinv;
        e           = eig(U);
        DGD(i,j)    = abs(abs(e(1)-e(2)));
    end
end

%% Histogram display
tau_m   = PMD*sqrt(L);
DGD_m   = mean(DGD);
DGD_s   = 1/Ntt*sum((DGD-DGD_m).^2);

tau_rel_err = (mean(DGD_m)-tau_m)/tau_m;
Nbins   = 50;
fig = figure;
    set(fig,"visible","on")
    H  = histogram(DGD*1e12,Nbins,'Normalization','pdf');
    xlabel("$\Delta\tau$ [ps]","interpreter","latex")
    ylabel("$p(\Delta\tau)$","interpreter","latex")
    title("Probability Distribution Function $(p)$ of the Differential Group Delays $(\Delta\tau)$","interpreter","latex")
    Values   = H.Values;
    BinEdges = H.BinEdges;
    Bins     = (BinEdges(1:end-1)+BinEdges(2:end))/2;
    
    DGD_mdn  = median(DGD(1,:)*1e12);
    DGD_rms  = rms(DGD(1,:)*1e12);

    axis([Bins(1),1.75*Bins(end),0,1.1*max(Values)])
    text(1.1*Bins(end),0.9*max(Values),strcat("DGD median = ",num2str(DGD_mdn,"%0.1f")," [ps]"),"interpreter","latex")
    text(1.1*Bins(end),0.85*max(Values),strcat("DGD rms = ",num2str(DGD_rms,"%0.1f")," [ps]"),"interpreter","latex")

    text(1.1*Bins(end),0.75*max(Values),strcat("L = ",num2str(L/1e3,"%i")," [$km$]"),"interpreter","latex")
    text(1.1*Bins(end),0.7*max(Values),strcat("PMD = ",num2str(PMD*10^(13.5),"%0.1f")," [$ps/\sqrt{km}$]"),"interpreter","latex")
    text(1.1*Bins(end),0.65*max(Values),strcat("$\lambda$ = ",num2str(lambda*1e9),' [nm]'),"interpreter","latex")
    
    text(1.1*Bins(end),0.55*max(Values),strcat("N seg = ",num2str(Np)),"interpreter","latex")
    text(1.1*Bins(end),0.5*max(Values),strcat("N rea = ",num2str(Nt*Ntt/1000),'k'),"interpreter","latex")

    if alphaRand == 0
        text(1.1*Bins(end),0.4*max(Values),strcat("$\alpha = cte$"),"interpreter","latex")
        name = "alphaCte";
    else
        text(1.1*Bins(end),0.4*max(Values),strcat("$\alpha \hookrightarrow \mathcal{U}[0,2\pi]$"),"interpreter","latex")
        name = "alphaRand";
    end
    if phiRand == 0
        text(1.1*Bins(end),0.35*max(Values),strcat("$\phi = cte$"),"interpreter","latex")
        name = strcat(name," phiCte");
    else
        text(1.1*Bins(end),0.35*max(Values),strcat("$\phi \hookrightarrow \mathcal{U}[0,2\pi]$"),"interpreter","latex")
        name = strcat(name," phiRand");
    end
    if dlRand == 0
        text(1.1*Bins(end),0.3*max(Values),strcat("$dL = cte$"),"interpreter","latex")
        name = strcat(name," dlCte");
    else
        text(1.1*Bins(end),0.3*max(Values),strcat("$dL \hookrightarrow \mathcal{N}[L/N_p,L/(5N_p)]$"),"interpreter","latex")
        name = strcat(name," dlRand");
    end

    saveas(H,strcat("DGD distribution - ",name, sprintf(" - L_%i PMD_%.1f mdn_%.1f.png",L/1e3,PMD*10^(13.5),DGD_mdn)))

%% True Maxwelian Fit --- naive version
% alpha0  = 1e-13;
% alpha   = alpha0;
% dalpha  = 1*1e-13;
% 
% x       = linspace(0,25*1e-12,Nbins);
% y       = sqrt(2/pi)/alpha^3*(x.^2).*exp(-x.^2/(2*alpha^2));
% max_tmp = max(Values)/max(y);
% y       = y*max_tmp;
% 
% rmse_tmp=10;
% rmse_eps=0.1;
% 
% count       = 1;
% count_max   = round(round(Bins(end)*1e-12/dalpha)/3);
% RMSE_val    = rmse_tmp;
% figure
%     while rmse_tmp>rmse_eps
% 
%         clf
% 
%         alpha   = alpha+dalpha;
%         y       = sqrt(2/pi)/alpha^3*(x.^2).*exp(-x.^2/(2*alpha^2));
%         max_tmp = max(Values)/max(y);
%         y       = y*max_tmp;
% 
%         hold on
%         scatter(Bins, Values)
%         plot(x*1e12,y)
%         axis([0,25,0,.21])
%         title(sprintf("%i --- %.1f",count,rmse_tmp))
% 
%         rmse_tmp= RMSE(Values,y)*100;
%         RMSE_val=[RMSE_val,rmse_tmp];
%  
%         count     = count+1;
% 
%         if count >count_max
%             break;
%         end
%         pause(1e-2)
%     end
% 
% [val,ind]   = min(RMSE_val);
% XX          = linspace(0,count,length(RMSE_val));
% alpha_opt   = alpha0+ind*dalpha;
% y_opt       = sqrt(2/pi)/alpha_opt^3*(x.^2).*exp(-x.^2/(2*alpha_opt^2));
% max_opt     = max(Values)/max(y_opt);
% y_opt       = y_opt*max_opt;
% 
% figure
% subplot(1,2,1)
%     hold on
%     plot(XX,RMSE_val)
%     scatter(XX(ind),RMSE_val(ind),100,"markerfacecolor",'k')
%     text(XX(ind),2*RMSE_val(ind),'h')
% subplot(1,2,2)
%     hold on
%     scatter(Bins, Values)
%     plot(x*1e12,y_opt)


%% True Maxwellian fit --- adaptative version 1

alpha0  = 1e-13;
alpha   = alpha0;
dalpha  = 1*1e-13;
R       = [0,alpha0];

x       = linspace(0,25*1e-12,Nbins);
y       = sqrt(2/pi)/R(end)^3*(x.^2).*exp(-x.^2/(2*R(end)^2));
max_tmp = max(Values)/max(y);
y       = y*max_tmp;

rmse_tmp=10;
rmse_eps=1.0;

count       = 1;
count_max   = round(round(Bins(end)*1e-12/dalpha)/3);
RMSE_val    = rmse_tmp;
figure
    while rmse_tmp>rmse_eps

        clf

        R       = [R,AGD(R(end-1:end),0.1,0)];
        y       = sqrt(2/pi)/(alpha0+R(end))^3*(x.^2).*exp(-x.^2/(2*(alpha0+R(end))^2));
        max_tmp = max(Values)/max(y);
        y       = y*max_tmp;

        hold on
        scatter(Bins, Values)
        plot(x*1e12,y)
        axis([0,25,0,.21])
        title(sprintf("%i --- %.1f",count,rmse_tmp))

        rmse_tmp= RMSE(Values,y)*100;
        RMSE_val=[RMSE_val,rmse_tmp];
 
        count     = count+1;

        if count >count_max
            break;
        end
        pause(1e-2)
    end

close all
[val,ind]   = min(RMSE_val);
XX          = linspace(1,count,count);
alpha_opt   = alpha0+R(end);
y_opt       = sqrt(2/pi)/alpha_opt^3*(x.^2).*exp(-x.^2/(2*alpha_opt^2));
max_opt     = max(Values)/max(y_opt);
y_opt       = y_opt*max_opt;

figure
subplot(1,2,1)
    hold on
    plot(XX,RMSE_val)
    scatter(XX(ind),val,100,"markerfacecolor",'k')
    text(0.25*XX(ind),val,sprintf("RMSE opt = $%.1f$ \\%%",val),"interpreter","latex")
    xlabel("number of epochs","interpreter","latex")
    ylabel("Root Mean Squared Error $\%$","interpreter","latex")
subplot(1,2,2)
    hold on
    scatter(Bins, Values)
    plot(x*1e12,y_opt)
    xlabel("$\Delta\tau$ [ps]","interpreter","latex")
    ylabel("$p(\Delta\tau)$","interpreter","latex")    
