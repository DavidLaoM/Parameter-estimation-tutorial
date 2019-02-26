% Giving the algorithm narrow boundaries. To show that the solution is
% reached faster.
clc, clear, close all

    % Setup
rng default
load('testdata.mat');
xdata = data.x_ENO;
ydata = data.v_ENO;
tdata = data.t_ENO;
% options = optimoptions('lsqnonlin','Display','iter');

ptrue   = [365.806 6.7 0.04 0.5];
p0      = [500 500 500 500];
lb      = [100 1 0 0];
ub      = [500 500 500 500];

    % Optimization
options = optimoptions('lsqnonlin','Display','iter','MaxFunEvals',1000)
[p3,resnorm3,residual3,exitflag3,output3,lambda3,jacobian3] = lsqnonlin(@ENOFitCost,p0,lb,ub,options,xdata,ydata);
disp(p3)

% Reached at 92 iterations, 465 function evaluations.

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
