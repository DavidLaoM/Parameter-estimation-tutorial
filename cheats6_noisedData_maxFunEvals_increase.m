% More iterations or function evaluations allowed, to see the final values
% obtained.
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
% options = optimoptions('lsqnonlin','Display','iter');
options = optimoptions('lsqnonlin','Display','iter','MaxFunEvals',1000)

% Optimization
[p9,resnorm9,residual9,exitflag9,output9,lambda9,jacobian9] = lsqnonlin(@ENOFitCost,p0,lb,ub,options,xdata,ydata);
disp('p9');
[p10,resnorm10,residual10,exitflag10,output10,lambda10,jacobian10] = lsqnonlin(@ENOFitCost,p0,lb,ub,options,xdata,ydata_1);
disp('p10');
[p11,resnorm11,residual11,exitflag11,output11,lambda11,jacobian11] = lsqnonlin(@ENOFitCost,p0,lb,ub,options,xdata,ydata_2);
disp('p11');
[p12,resnorm12,residual12,exitflag12,output12,lambda12,jacobian12] = lsqnonlin(@ENOFitCost,p0,lb,ub,options,xdata,ydata_5);
disp('p12');
disp(p9);
disp(p10);
disp(p11);
disp(p12);

% As expected, the more the  noise, the worse the parameter guess. At 5%
% noise, several parameters are already out of order.

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