clearvars
close all

%%

Cinit = zeros(4, 1);
Cinit(1, 1) = 1400;
Cinit(2, 1) = 0;
Cinit(3, 1) = 3400;
Cinit(4, 1) = 0;

load FullModelDataProp.mat;

fun = @(kValues) reducedModelError(kValues, C1, C2, C3, C4, Cinit); 
rng default
nvars = 3; 

% K.propagation = kValues(1);
% K.in1 = kValues(2);
% K.in2 = kValues(3);

% [K.propagation, K.in1, K.in2];

kV0 = [8e-5, 6e-6, 2e-8]; 
lb = [8e-6, 6e-7, 2e-9]; 
ub = [1e-4, 6e-5, 8e-7]; 

options = optimoptions('particleswarm', 'SwarmSize', 20, 'Display','iter');

tic
[kOptimized,fval,exitflag,output] = particleswarm(fun, nvars, lb, ub, options);
toc

% tic
% A = []; 
% b = [];
% Aeq = [];
% beq = []; 
% [kOptimized,fval,exitflag,output] = fmincon(fun, kV0, A, b, Aeq, beq, lb, ub);
% toc

% options = optimoptions('particleswarm', 'SwarmSize', 100, ...
%     'HybridFcn', {@fmincon},'Display','iter', ...
%     'FunctionTolerance',0.05, 'MaxStalliterations',1, ...
%     'MaxIterations',2);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% This was used in thesis
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% lb = [0.36e-5, 0.67e-6, 1e-8];
% % lb = [0.36e-4, 0.67e-5, 1e-7];
% ub = [0.8e-3, 0.67e-4, 1e-6];

% options = optimoptions('particleswarm', 'SwarmSize', 20, ...
%     'HybridFcn', {@fmincon},'Display','iter', 'FunctionTolerance', 100);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% tic
% [kOptimized,fval,exitflag,output] = particleswarm(fun, nvars, lb, ub, options);
% toc


C0 = zeros(4, 1);
C0(1, 1) = Cinit(1);
C0(2, 1) = Cinit(2);
C0(3, 1) = Cinit(3);
C0(4, 1) = Cinit(4);
tRange = 0:1:1200; % s

fun = @(t,y) reducedModelInput(t, y, kOptimized);
options = odeset('AbsTol', 1e-3*ones(1,4)); % in Nano Moles
[T, C] = ode23s(fun, tRange, C0, options);

%%
figure(653); 
plot(T,C)

figure(654)
plot(T,C(:,2))
hold on
plot(T,C2*1e09)

%%
figure(655)
plot( 1e9*C2(2:end) - 1e9*C2(1:end-1),  C(2:end,2) - C(1:end-1,2), 'o')
hold on, plot(linspace(-5,5,100), linspace(-5,5,100))
axis equal