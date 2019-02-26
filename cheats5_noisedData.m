% Adding noise to the data. The more the noise, the worse the parameter
% estimate. The number of functions evaluations is not indicative of a good
% fit in this case.
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

ptrue   = [365.806 6.7 0.04 0.5];
p0      = [500 500 500 500];
lb      = [0 0 0 0];
ub      = [1000 1000 1000 1000];
options = optimoptions('lsqnonlin','Display','iter');

    % Optimization
[p5,resnorm5,residual5,exitflag5,output5,lambda5,jacobian5] = lsqnonlin(@ENOFitCost,p0,lb,ub,options,xdata,ydata);
disp('p5');
[p6,resnorm6,residual6,exitflag6,output6,lambda6,jacobian6] = lsqnonlin(@ENOFitCost,p0,lb,ub,options,xdata,ydata_1);
disp('p6');
[p7,resnorm7,residual7,exitflag7,output7,lambda7,jacobian7] = lsqnonlin(@ENOFitCost,p0,lb,ub,options,xdata,ydata_2);
disp('p7');
[p8,resnorm8,residual8,exitflag8,output8,lambda8,jacobian8] = lsqnonlin(@ENOFitCost,p0,lb,ub,options,xdata,ydata_5);
disp('p8');
disp(p5);
disp(p6);
disp(p7);
disp(p8);

% Some estimations already finished, but parameter values are not that
% nice. As expected, the error is greater when more noise there.

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