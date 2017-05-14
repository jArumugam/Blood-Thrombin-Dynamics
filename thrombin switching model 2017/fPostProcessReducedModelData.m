function fPostProcessReducedModelData()

% caseNames = {'meanPhysiologic', 'ACSmean'};
caseNames = {'meanPhysiologic', 'ACSmean', 'hemophilia2c'};
caseColors = {[0 0.4470 0.7410], ... 
    [0.8500 0.3250 0.0980], ... 
    [0.9290 0.6940 0.1250]};

for ii=1:length(caseNames)
    fPlotModelData(caseNames{ii}, caseColors{ii});
end

end

function fPlotModelData(caseName, caseColor)

InitTF = 5e-12;
t0 = 0; tf = 2000;
tspan = t0:1:tf;

fRates = @SolverTf.fReaction42Rates2002;
% caseName = caseNames{ii};

[T, C] = SolverTf.fSolveThrombinGenerationNew(tspan, fRates, caseName, InitTF);

figure(40815)
hold on
plot(T,C(:,7)*1e09 + C(:,25)*1e09, '-.', 'color', caseColor);


% fName1 = strcat(caseName, '_FullModelData.mat');
fName2 = strcat(caseName, '_ReducedModelData.mat');
load(fName2)

tRange = 0:1:2000; % s
fun = @(t,y) reducedModelInput(t, y, datRM.kOptimized);
options = odeset('AbsTol', 1e-3*ones(1,4)); % in Nano Moles
[T_R, C_R] = ode23s(fun, tRange, datRM.C0, options);

plot(T_R, C_R(:,2), '-',  'color', caseColor)


end