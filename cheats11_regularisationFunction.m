clc, clear, close all

    % Setup
rng default
load('testdata.mat');
xdata = data.x_ENO;
ydata = data.v_ENO;
tdata = data.t_ENO;

  
% Setup
options = optimoptions('lsqnonlin','Display','iter','MaxFunEvals',1000);
pENO    = [365.806 6.7 0.04 0.5];
p0      = [500 500 500 500];
lb      = [0 0 0 0];
ub      = [1000 1000 1000 1000];

[p29,resnorm29,residual29,exitflag29,output29,lambda29,jacobian29] = lsqnonlin(@ENOFitCost,p0,lb,ub,options,xdata,ydata);
disp('p29');
[p30,resnorm30,residual30,exitflag30,output30,lambda30,jacobian30] = lsqnonlin(@ENOFitCost_Reg,p0,lb,ub,options,xdata,ydata);
disp('p30');
disp(p29);
disp(p30);

% The regularization make 170 iterations go down to 36.

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

% Enolase cost funciton + regularisation over the value of Vmax p(1)
function e = ENOFitCost_Reg(p,x,y)
v = ENO(p,x);
e2 = (p(1) - 365.806) * 0.1;
e = y - v + e2;
% e = (abs(y-v)./y).^2;
end