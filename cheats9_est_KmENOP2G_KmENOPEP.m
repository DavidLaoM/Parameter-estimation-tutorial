% Km and vmax fixed. Estimate KmP2G abd KmPEP.
clc, clear, close all

    % Setup
rng default
load('testdata.mat');
xdata = data.x_ENO;
ydata = data.v_ENO;
tdata = data.t_ENO;

    % Generate the noised data
    ydata_1 = ydata + ydata .* 0.01 .* randn(size(ydata));
    ydata_2 = ydata + ydata .* 0.02 .* randn(size(ydata));
    ydata_5 = ydata + ydata .* 0.05 .* randn(size(ydata));
    
    % Delete part of dataset (minimum region)
    tdata2      = [tdata(1:6); tdata(12:end)];
    ydata2      = [ydata(1:6); ydata(12:end)];
    ydata2_1    = [ydata_1(1:6); ydata_1(12:end)];
    ydata2_2    = [ydata_2(1:6); ydata_2(12:end)];
    ydata2_5    = [ydata_5(1:6); ydata_5(12:end)];
    xdata2 = [xdata(1:6,:); xdata(12:end,:)];
   
% Setup
options = optimoptions('lsqnonlin','Display','iter','MaxFunEvals',1000);
pENO    = [365.806 6.7 0.04 0.5];
p0      = [365.806 6.7 500 500];
lb      = [365.806 6.7 0 0];
ub      = [365.806 6.7 1000 1000];

[p21,resnorm21,residual21,exitflag21,output21,lambda21,jacobian21] = lsqnonlin(@ENOFitCost,p0,lb,ub,options,xdata2,ydata2);
disp('p21');
[p22,resnorm22,residual22,exitflag22,output22,lambda22,jacobian22] = lsqnonlin(@ENOFitCost,p0,lb,ub,options,xdata2,ydata2_1);
disp('p22');
[p23,resnorm23,residual23,exitflag23,output23,lambda23,jacobian23] = lsqnonlin(@ENOFitCost,p0,lb,ub,options,xdata2,ydata2_2);
disp('p23');
[p24,resnorm24,residual24,exitflag24,output24,lambda24,jacobian24] = lsqnonlin(@ENOFitCost,p0,lb,ub,options,xdata2,ydata2_5);
disp('p24');
disp(p21);
disp(p22);
disp(p23);
disp(p24);

% Estimations take 33, 36, 40 51 iterations. Parameter 4 estimation
% continues off.

%% Kinetics
% Enolase
function v = ENO(p,x)
    v = (p(1).*(x(:,1) - x(:,2)./p(2)))./(p(3).*(1 + x(:,1)./p(3) + x(:,2)./p(4)));
end
% p(1) = VmENO;     365.806 
% p(2) = KeqENO;    6.7
% p(3) = KmENOP2G;  0.04
% p(4) = KmENOPEP;  0.5
% x(1) = P2G;
% x(2) = PEP;

% Enolase cost funciton
function e = ENOFitCost(p,x,y)
v = ENO(p,x);
e = y - v;
% e = (abs(y-v)./y).^2;
end