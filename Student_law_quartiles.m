clear all;
close all;
clc;

k   = [5:5:20,30:10:80,100,120,121];
i   = k+1;

t95 = [2.015,1.812,1.753,...
      1.725,1.708,1.697,...
      1.684,1.676,1.671,...
      1.664,1.660,1.658,1.645];
      
t99 = [3.365,2.764,2.602,...
      2.528,2.485,2.457,...
      2.423,2.403,2.390,...
      2.374,2.364,2.358,2.326];

xp    = linspace(0,121,122);

% power      
p95l  = polyfit(k,log(t95),7);
p99l  = polyfit(k,log(t99),7);

yp95l = polyval(p95l,xp);
yp99l = polyval(p99l,xp);

yp95lb= exp(yp95l)(i);
yp99lb= exp(yp99l)(i);

rmse95l= RMSE(yp95lb,t95)
rmse99l= RMSE(yp99lb,t99)

yp95lf= exp(yp95l)(xp(4):xp(end));
yp99lf= exp(yp99l)(xp(4):xp(end));

% polynomial
p95   = polyfit(k,t95,7);
p99   = polyfit(k,t99,7);

yp95  = polyval(p95,xp);
yp99  = polyval(p99,xp);

yp95b = yp95(i);
yp99b = yp99(i);

rmse95= RMSE(yp95b,t95)
rmse99= RMSE(yp99b,t99)

%figure
%hold on
%plot(k,t95)
%plot(k,yp95b)
%
%plot(k,t99)
%plot(k,yp99b)
%
%
M = [xp(4:end)',yp95lf',yp99lf'];
csvwrite("student_law_coefficients.txt",M)