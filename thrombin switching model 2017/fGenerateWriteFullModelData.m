function fGenerateWriteFullModelData()

caseNames = {'meanPhysiologic', 'ACSmean', 'hemophilia2c'};

for ii=1:length(caseNames)
    
    InitTF = 5e-12; 
    t0 = 0; tf = 1200; 
    tspan = t0:1:tf; 
    
    fRates = @SolverTf.fReaction42Rates2002;
    caseName = caseNames{ii};
    
    [T, C] = SolverTf.fSolveThrombinGenerationNew(tspan, fRates, caseName, InitTF); 
    
    fWriteFullModelData(T,C,caseName, InitTF)
    
end

end

function fWriteFullModelData(T,C,caseName, InitTF)

dat.CaseName = caseName; 

dat.T = T;

dat.C1 = C(:,14); 
dat.C2 = C(:,7) + C(:,25); 
dat.C3 = C(:,29); 
dat.C4 = C(:,33) + C(:,31); 

dat.InitTF = InitTF;


% Plot Full Model Results
% X14: II
figure(4637)
hold on
% X7: IIa
% plot Thrombin
plot(T, dat.C2, 'LineWidth',1)
axis([0 3600 0 0.6e-6])
grid on

fName = strcat(caseName, '_FullModelData.mat');
save(fName, 'dat' )

end