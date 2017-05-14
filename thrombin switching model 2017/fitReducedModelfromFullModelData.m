function fitReducedModelfromFullModelData()

% caseNames = {'meanPhysiologic', 'ACSmean'};
caseNames = {'meanPhysiologic', 'ACSmean', 'hemophilia2c'};

for ii = 1:length(caseNames)
    
    fitWriteReducedModelUsingPSO(caseNames{ii});
    
end

end

function fitWriteReducedModelUsingPSO(caseName)

% Cinit = zeros(4, 1);
% Cinit(1, 1) = 1400;
% Cinit(2, 1) = 0;
% Cinit(3, 1) = 3400;
% Cinit(4, 1) = 0;

fName = strcat(caseName, '_FullModelData.mat');
load(fName);

Cinit = zeros(4, 1);
Cinit(1, 1) = dat.C1(1)*1e09;
Cinit(2, 1) = 0;
Cinit(3, 1) = dat.C3(1)*1e09;
Cinit(4, 1) = 0;


fun = @(kValues) reducedModelError(kValues, dat.C1, dat.C2, dat.C3, dat.C4, Cinit); 
rng default
nvars = 3; 

lb = [8e-6, 6e-7, 2e-9]; 
ub = [1e-4, 6e-5, 8e-7]; 
options = optimoptions('particleswarm', 'SwarmSize', 20, 'Display','iter');

tic
[kOptimized,fval,exitflag,output] = particleswarm(fun, nvars, lb, ub, options);
toc

datRM = solveReducedModel(kOptimized, Cinit, caseName);

datRM.kOptimized = kOptimized;
datRM.fval = fval;
datRM.exitflag = exitflag; 
datRM.output = output;

fOutputName = strcat(caseName, '_ReducedModelData.mat'); 
save(fOutputName, 'datRM')

end


function [datRM] = solveReducedModel(kOptimized, Cinit, caseName)

datRM.C0 = zeros(4, 1);
datRM.C0(1, 1) = Cinit(1);
datRM.C0(2, 1) = Cinit(2);
datRM.C0(3, 1) = Cinit(3);
datRM.C0(4, 1) = Cinit(4);
datRM.tRange = 0:1:1200; % s

fun = @(t,y) reducedModelInput(t, y, kOptimized);
options = odeset('AbsTol', 1e-3*ones(1,4)); % in Nano Moles 
tRange = 0:1:1200; 
[T, C] = ode23s(fun, tRange, datRM.C0, options);
datRM.T = T;
datRM.C = C;

datRM.caseName = caseName;

% figure(653) 
% hold on
% plot(datRM.T,datRM.C)

figure(654)
hold on 
plot(datRM.T,datRM.C(:,2)) 

end