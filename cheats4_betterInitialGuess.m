% In a similar way to the previous, a closer initial guess should reduce
% the time required as well.
clc, clear, close all

    % Setup
rng default
load('testdata.mat');
xdata = data.x_ENO;
ydata = data.v_ENO;
tdata = data.t_ENO;
% options = optimoptions('lsqnonlin','Display','iter');

ptrue   = [365.806 6.7 0.04 0.5];
p0      = [500 20 0.1 0.2];
lb      = [0 0 0 0];
ub      = [1000 1000 1000 1000];

    % Optimization
options = optimoptions('lsqnonlin','Display','iter','MaxFunEvals',1000);
[p4,resnorm4,residual4,exitflag4,output4,lambda4,jacobian4] = lsqnonlin(@ENOFitCost,p0,lb,ub,options,xdata,ydata);
disp(p4)

% It takes less, 54 iterations and 274 function evaluations.

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