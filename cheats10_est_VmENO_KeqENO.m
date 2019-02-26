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
p0      = [500 500 0.04 0.5];
lb      = [0 0 0.04 0.5];
ub      = [1000 1000 0.04 0.5];

[p25,resnorm25,residual25,exitflag25,output25,lambda25,jacobian25] = lsqnonlin(@ENOFitCost,p0,lb,ub,options,xdata2,ydata2);
disp('p25');
[p26,resnorm26,residual26,exitflag26,output26,lambda26,jacobian26] = lsqnonlin(@ENOFitCost,p0,lb,ub,options,xdata2,ydata2_1);
disp('p26');
[p27,resnorm27,residual27,exitflag27,output27,lambda27,jacobian27] = lsqnonlin(@ENOFitCost,p0,lb,ub,options,xdata2,ydata2_2);
disp('p27');
[p28,resnorm28,residual28,exitflag28,output28,lambda28,jacobian28] = lsqnonlin(@ENOFitCost,p0,lb,ub,options,xdata2,ydata2_5);
disp('p28');
disp(p25);
disp(p26);
disp(p27);
disp(p28);

% Estimations take 24, 24, 26 22 iterations. Values are better, compared to
% the other couple.

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