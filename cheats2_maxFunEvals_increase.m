% THHe number of function evaluations in increased so the algorithm has the
% chance to improve further the results
clc, clear, close all

    % Setup
rng default
load('testdata.mat');
xdata = data.x_ENO;
ydata = data.v_ENO;
tdata = data.t_ENO;
options = optimoptions('lsqnonlin','Display','iter');

ptrue   = [365.806 6.7 0.04 0.5];
p0      = [500 500 500 500];
lb      = [0 0 0 0];
ub      = [1000 1000 1000 1000];

    % Optimization
options = optimoptions('lsqnonlin','Display','iter','MaxFunEvals',1000);
[p2,resnorm2,residual2,exitflag2,output2,lambda2,jacobian2] = lsqnonlin(@ENOFitCost,p0,lb,ub,options,xdata,ydata);
disp(p2)
% v_ENO_sim = ENO(p,xdata);
% v_ENO_sim2 = ENO(p2,xdata);

% Good guess found at 170 iterations, 855 funciton evaluations. 

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
