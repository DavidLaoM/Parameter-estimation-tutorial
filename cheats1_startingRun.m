% Running parameter eestimation over enolase. Simple case: bounds really
% open and bad initial guess

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

    % Parameter estimation
[p,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(@ENOFitCost,p0,lb,ub,options,xdata,ydata);
disp(p)
v_ENO_sim = ENO(p,xdata);

    % Plotting
figure
plot(tdata,v_ENO_sim)
hold on
plot(tdata, ydata, 'o')
title('Enolase reaction rate fit')
xlabel('time [min]')
ylabel('reaction rate [min-1]')
legend('simulated data', 'mock data')

% The maximum number of funciton evaluations reaches its top. Simulation 
% plot works quite ok, but the minimum point is missed, for example.
% However, see that parameter values are still far away from ideal.

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

