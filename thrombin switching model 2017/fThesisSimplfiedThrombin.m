function [t,y] = fThesisSimplfiedThrombin()

y0 = [0 1400 0.2*3400 0]; % [IIa II ATIII II-ATIII];
% y0 = [0 1400 1200 0]; % [IIa II ATIII II-ATIII];
tspan = 0:1:1200;

% K.surf = 3.6e-6;
% K.propagation = 0.36e-4;
% K.in1 = 0.67e-5;
% K.in2 = 1e-7;
% % K.in2 = 1e-5;
% K.IIaThres = 2; % nM

% from PSO
K.surf = 3.6e-6;
K.propagation = 3.976e-5;
K.in1 = 6.13e-6;
K.in2 = 1.688e-7;
% K.in2 = 1e-5;
K.IIaThres = 2; % nM

[t,y] = fThesisSimplfiedThrombinSolve(tspan, y0, K);


fPlotModelComparison(t,y)

end


function fPlotModelComparison(t,y)

load FullModelDataProp.mat;

figure(5324)
% title('Comparison of the simplified and reduced model')
subplot(4,1,1)
plot(t, y(:,1), 'LineWidth',2)
xlabel('Time (sec)', 'FontSize', 14)
% ylabel('Concentration (M)', 'FontSize', 14)
hold on
plot(t, C2*1e09, 'LineWidth',2)
ylabel('Factor II (M)', 'FontSize', 14)
xlabel('Time (s)', 'FontSize', 14)
grid on, box on
h = legend({'Simplified', 'Full'});
set(h,'Location', 'west')

subplot(4,1,2)
% X7: IIa
plot(t, y(:,2), 'LineWidth',2)
hold on
plot(t, C1*1e09, 'LineWidth',2)
ylabel('Thrombin (M)', 'FontSize', 14)
xlabel('Time (s)', 'FontSize', 14)
grid on, box on
h = legend({'Simplified', 'Full'});
set(h,'Location', 'west')

subplot(4,1,3)
% x29: AT
plot(t, y(:,3), 'LineWidth',2)
hold on
plot(t, C3*1e09, 'LineWidth',2)
grid on, box on
ylabel('ATIII (M)', 'FontSize', 14)
xlabel('Time (s)', 'FontSize', 14)
h = legend({'Simplified', 'Full'}); 
set(h,'Location', 'west')

subplot(4,1,4)
% X33: IIa - ATIII
plot(t, y(:,4), 'LineWidth',2)
hold on
plot(t, C4*1e09, 'LineWidth',2)
ylabel('TAT (M)', 'FontSize', 14)
xlabel('Time (s)', 'FontSize', 14)
grid on, box on
h = legend({'Simplified', 'Full'});
set(h,'Location', 'west')

% hLegend = legend('II', 'Thrombin', 'ATIII', 'Thrombin-ATIII', ...
%     'Location','eastoutside');
% set(hLegend,'FontSize',12)
grid on, box on

%%
fName = strcat('Users/jayavelarumugam/Dropbox/J/2017/Thesis/Thesis/figures/ch6SimplifiedModel','ModelComparison');
% saveas(gcf,'fName','png')
% saveas(gcf,'fName','eps')
% saveas(gcf,'fName','pdf')
% saveas(gcf,'fName','tif')

end


function [t,y] = fThesisSimplfiedThrombinSolve(tspan, y0, K)

% [t,y] = ode23(@(t,y) fThesisSimplfiedThrombinRate(t,y,K), tspan, y0);
[t,y] = ode23t(@(t,y) fThesisSimplfiedThrombinRate(t,y,K), tspan, y0);

%
figure(530); plot(t,y(:,1),'-'), hold on
% figure(22); plot(t,y(:,2),'-'), hold on
figure(23); plot(t,y(:,3),'-'), hold on
% figure(23); plot(t,y(:,4),'-'), hold on

figure(531);
subplot(4,1,1), 
plot(t,y(:,1),'-'), hold on, grid on
title('IIa (nM)')

subplot(4,1,2), 
plot(t,y(:,2),'-'), hold on, grid on 
title('II (nM)')

subplot(4,1,3), plot(t,y(:,3),'-'), hold on, grid on
title('ATIII (nM)')

subplot(4,1,4), plot(t,y(:,4),'-'), hold on, grid on
title('IIa-ATIII (nM)')

end



function [ksurf] = fKsurfVal(t, Kksurf)

if(t<= 6000)
   ksurf = Kksurf;
else
   ksurf = 0; 
end

end

function [dy] = fThesisSimplfiedThrombinRate(t,y, K)

Ksurft = fKsurfVal(t, K.surf);

if (y(1) <= K.IIaThres)
    ksurf = Ksurft;
    kin = K.in2;
    kpropagation = 0;
    
elseif ( y(1) > K.IIaThres)
    ksurf = 0;
    kin = K.in1; 
    kpropagation = K.propagation; 
end

dy = zeros(4,1);
% [IIa II AT];
% dy(1) = -kin*y(1)*y(3) + (ksurf + kpropagation)*y(2);
% dy(2) = -(ksurf + kpropagation)*y(2);
% dy(3) = -kin*y(1)*y(3);
% dy(4) = kin*y(1)*y(3);

dy(1) = -kin*y(1)*y(3) + ksurf*y(2) + kpropagation*y(2)*y(1);
dy(2) = -ksurf*y(2) - kpropagation*y(2)*y(1);
dy(3) = -kin*y(1)*y(3);
dy(4) = kin*y(1)*y(3);



end