clc, clear, close all

    % Setup
rng default
load('testdata.mat');
xdata = data.x_ENO;
ydata = data.v_ENO;
tdata = data.t_ENO;

x_ENO = xdata;
v_ENO = ydata;

ptrue = [365.806 6.7 0.04 0.5];
[p1range,p2range] = meshgrid(200:0.5:600.1,2:.0125:12);            % More datapoints for pretty plots
[p3range,p4range] = meshgrid(0.01:0.00005:0.11, 0.1:0.0005:1.1);  % More datapoints for pretty plots
% [p1range,p2range] = meshgrid(200:5:600.1,2:.125:12);
% [p3range,p4range] = meshgrid(0.01:0.0005:0.11, 0.1:0.005:1.1);
% [p1range,p2range] = meshgrid(0.01:1:1000.01,0.01:1:1000.01);              % More datapoints for pretty plots
% [p3range,p4range] = meshgrid(0.01:1:1000.01,0.01:1:1000.01);              % More datapoints for pretty plots
% [p1range,p2range] = meshgrid(0.01:10:1000.01,0.01:10:1000.01);
% [p3range,p4range] = meshgrid(0.01:10:1000.01,0.01:10:1000.01);

    V=zeros(length(p1range), length(p2range)); plt = 0;
    h = waitbar(0,'Please wait...');
    for k=1:length(p1range)
        ptemp = ptrue;
        ptemp(1) = p1range(1,k);
        for l=1:length(p2range)
            ptemp(2) = p2range(l);
            nutmp = ENOtest(ptemp,x_ENO);
            nutrue= ENOtest(ptrue,x_ENO);
            e = nutrue-nutmp;
            V(k,l) = e'*e/sqrt(e'*e);%(or take out, changed )%/length(e);
        end
        waitbar(k/length(p1range),h)
    end
    close(h)

    V2=zeros(length(p3range), length(p4range)); plt = 0;
    h = waitbar(0,'Please wait...');
    for k=1:length(p3range)
        ptemp = ptrue;
        ptemp(3) = p3range(1,k);
        for l=1:length(p4range)
            ptemp(4) = p4range(l);
            nutmp = ENOtest(ptemp,x_ENO);
            nutrue= ENOtest(ptrue,x_ENO);
            e2 = nutrue-nutmp;
            V2(k,l) = e2'*e2/sqrt(e2'*e2);%(or take out, changed )%/length(e);
        end
        waitbar(k/length(p3range),h)
    end
    close(h)    
    
    
    figure; subplot(1,2,1);
    meshc(p1range, p2range, log10(V')) %V) %log10(V'))
    set(gca, 'FontSize', 14)
    xlabel('V_m'); ylabel('K_eq'); zlabel('log(V)')
    xlabel('V_m','FontSize', 16); ylabel('K_eq','FontSize', 16)
    subplot(1,2,2);
    meshc(p3range, p4range, log10(V2')) %V) %log10(V'))
    set(gca, 'FontSize', 14)
    xlabel('KmENOP2G'); ylabel('KmENOPEP'); zlabel('log(V2)')
    xlabel('KmENOP2G','FontSize', 16); ylabel('KmENOPEP','FontSize', 16)
    
    
function v = ENOtest(p,x)
    v = (p(1).*(x(:,1) - x(:,2)./p(2)))./(p(3).*(1 + x(:,1)./p(3) + x(:,2)./p(4)));
end
% p(1) = VmENO;     365.806 
% p(2) = KeqENO;    6.7
% p(3) = KmENOP2G;  0.04
% p(4) = KmENOPEP;  0.5
% x(1) = P2G;
% x(2) = PEP;
