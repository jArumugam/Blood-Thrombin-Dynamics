classdef SolverTf
    %
    % April 21 2015
    % Should the class name necessarily start with Captial letters?
    % tf gives error messages
    
    methods(Static)
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
        function [C0] =  fGetInitialContion(caseName, InitTF)
            % Initialize with mean physiological value
            C0 = zeros(34,1);
            C0(1,1) = InitTF;
            C0(2,1) = 1E-08; % VII
            C0(8,1) = 1.6E-07; % X
            C0(11,1) = 9E-08; % IX
            C0(14,1) = 1.4E-06; % II
            C0(15,1) = 7E-10; % VIII
            C0(21,1) = 2E-08; % V
            C0(26,1) =  2.5E-09; % TFPI
            C0(29,1) = 3.4E-06; % AT
            C0(4,1) = C0(2,1)*0.01; % VIIa
            
            switch caseName
                case 'meanPhysiologic'
                    C0(2,1) = C0(2,1); % VII
                    C0(8,1) = C0(8,1); % X
                    C0(11,1) = C0(11,1); % IX
                    C0(14,1) = C0(14,1); % II
                    C0(15,1) = C0(15,1); % VIII
                    C0(21,1) = C0(21,1); % V
                    C0(26,1) = C0(26,1); % TFPI
                    C0(29,1) = C0(29,1); % AT
                    C0(4,1) = C0(2,1)*0.01; % VIIa
                    
                case 'hemophilia1a' % Models for thrombin generation
                    C0(2,1) = 0.81*C0(2,1); % VII
                    C0(8,1) = 1.02*C0(8,1); % X
                    C0(11,1) = 1.01*C0(11,1); % IX
                    C0(14,1) = 1.08*C0(14,1); % II
                    C0(15,1) = 0.01*C0(15,1); % VIII
                    C0(21,1) = 1.11*C0(21,1); % V
                    C0(26,1) = 0.51*C0(26,1); % TFPI
                    C0(29,1) = 0.97*C0(29,1); % AT
                    C0(4,1) = C0(2,1)*0.01; % VIIa
                    
                case 'ACSmean' % Thrombin generation in ACS and CAD
                    C0(2,1) = 1.16*C0(2,1); % VII
                    C0(8,1) = 1.11*C0(8,1); % X
                    C0(11,1) = 1.17*C0(11,1); % IX
                    C0(14,1) = 1.23*C0(14,1); % II
                    C0(15,1) = 1.62*C0(15,1); % VIII
                    C0(21,1) = 1.21*C0(21,1); % V
                    C0(26,1) = 1.14*C0(26,1); % TFPI
                    C0(29,1) = 0.89*C0(29,1); % AT
                    C0(4,1) = C0(2,1)*0.01; % VIIa
                    
                case 'CADmean' % Thrombin generation in ACS and CAD
                    C0(2,1) = 1.21*C0(2,1); % VII
                    C0(8,1) = 1.14*C0(8,1); % X
                    C0(11,1) = 1.21*C0(11,1); % IX
                    C0(14,1) = 1.03*C0(14,1); % II
                    C0(15,1) = 1.23*C0(15,1); % VIII
                    C0(21,1) = 1.31*C0(21,1); % V
                    C0(26,1) = 1.02*C0(26,1); % TFPI
                    C0(29,1) = 1.13*C0(29,1); % AT
                    C0(4,1) = C0(2,1)*0.01; % VIIa
                    
                case 'hemophilia2c' % Models for thrombin generation
                    C0(2,1) = 1.16*C0(2,1); % VII
                    C0(8,1) = 1.05*C0(8,1); % X
                    C0(11,1) = 1.25*C0(11,1); % IX
                    C0(14,1) = 1.20*C0(14,1); % II
                    C0(15,1) = 0.01*C0(15,1); % VIII
                    C0(21,1) = 1.02*C0(21,1); % V
                    C0(26,1) = 0.71*C0(26,1); % TFPI
                    C0(29,1) = 1.38*C0(29,1); % AT
                    C0(4,1) = C0(2,1)*0.01; % VIIa
                    
                otherwise
                    disp('case not defined \n')
            end
            
        end
        
        
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
        function [T, C] =  fSolveThrombinGenerationNew(tspan, fRates, caseName, InitTF)
            
            [C0] =  SolverTf.fGetInitialContion(caseName, InitTF);
            %             figure(9945), hold on, plot(C0)
            options = odeset('AbsTol', 1e-15);
            [T, C] = ode23t(fRates, tspan, C0, options);
            C(C < 0) = 0;
        end % fSolveThrombinGeneration Ends
        
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
        function [T, C] =  fSolveThrombinGeneration(tspan, fRates)
            
            
            options = odeset('AbsTol', 1e-15);
            [T, C] = ode23t(fRates,tspan, C0, options);
            C(C < 0) = 0;
        end % fSolveThrombinGeneration Ends
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
        % fSolverForRAWDATA
        function [] = fSolveAndWriteData(Trigger, CaseName, NoOfDataPoints, fRates,...
                t0, tf, dt, flagRAW, NIntervals)
            
            
            NConst = dlmread('data/NConst.txt');
            InitTF = Trigger*1e-12;
            tspan = t0:dt:tf;
            
            % Solve and write solution for all data points
            hWaitBar = waitbar(0,'Please wait...');
            steps = NoOfDataPoints;
            
            for i = 1:NoOfDataPoints
                
                step = i;
                waitbar(step / steps);
                
                SampleNo = i;
                
                C0 = PreProsTf.fGetInitialConcentration(InitTF, CaseName, SampleNo);
                [T, C] = SolverTf.fSolveThrombinGeneration(C0,  tspan, fRates);
                
                %~~~~~~~~~~~~~~~~~%
                % writes RAW DATA %
                %~~~~~~~~~~~~~~~~~%
                % if flagRAW SolverTf.fWriteSolution(C, fNamePrefix, varNames); end
                
                % write thrombin parameters
                [Time2nM] = PostProsTf.fGetT2nM(T, C(:, 7));
                [MaxL, TMaxL, MaxR, TMaxR, AUC] = ...
                    PostProsTf.fGetSummaryParameters(T, C(:, 7));
                ThrombinData = [Time2nM, MaxL, TMaxL, MaxR, TMaxR, AUC];
                PostProsTf.fAppendRowVecToTextFile(ThrombinData, CaseName, ...
                    Trigger, 'IIa');
                
                % write Active Thrombin parameters
                [Time2nM] = PostProsTf.fGetT2nM(T, C(:, 7) + 1.2*C(:,25));
                [MaxL, TMaxL, MaxR, TMaxR, AUC] = ...
                    PostProsTf.fGetSummaryParameters(T, C(:, 7)+ 1.2*C(:,25));
                ActiveThrombinData = [Time2nM, MaxL, TMaxL, MaxR, TMaxR, AUC];
                PostProsTf.fAppendRowVecToTextFile(ActiveThrombinData, CaseName, ...
                    Trigger, 'ActiveThrombin');
                
                % write Xa paramaters
                [MaxL, TMaxL, MaxR, TMaxR, AUC] = ...
                    PostProsTf.fGetSummaryParameters(T, C(:, 6));
                XaData = [MaxL, TMaxL, MaxR, TMaxR, AUC];
                PostProsTf.fAppendRowVecToTextFile(XaData, CaseName, Trigger, 'Xa');
                
                % write Coeff data
                C = PostProsTf.fNormalizeCColumns(C, NConst);
                % CoeffRowVector = PostProsTf.fFindPchipCoeffsWithBadReshape(T, C', NIntervals, t0, tf);
                CoeffRowVector = PostProsTf.fFindPchipCoeffsWithBetterReshape(T, C, NIntervals, t0, tf);
                PostProsTf.fAppendRowVecToTextFile(CoeffRowVector, CaseName, Trigger, 'pchipCoeffs');
                
            end
            
            % Convert txt files to excel files
            %             PostProsTf.fReadTxtFileWriteXlsFile(CaseName, Trigger, 'Xa');
            %             PostProsTf.fReadTxtFileWriteXlsFile(CaseName, Trigger, 'IIa');
            %             PostProsTf.fReadTxtFileWriteXlsFile(CaseName, Trigger, 'ActiveThrombin');
            %             PostProsTf.fWritePchipDataToXlsFile(NIntervals, CaseName, Trigger);
            
            close(hWaitBar);
            
        end % fSolveAndWriteData Ends
        
        
        
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
        function [] = fWriteSolution(C, fNamePrefix, varNames)
            % sWriteSolution(C, fNamePrefix)
            % Writes each variable solution to a separate file
            % i = 7;
            for i = 1:34
                
                fName = strcat(fNamePrefix,'_Species_',varNames{1}{i},'.txt');
                dlmwrite(fName, ...
                    C(:,i)', ...
                    'delimiter','\t', '-append',  'newline', 'pc');
                
            end
            
            % fileID = fopen('SpeciesNames.txt');
            % varNames = textscan(fileID,'%s');
            % fclose(fileID);
            % celldisp(varNames);
            
            % strcat('rgtehte_',varNames{1}{34})
            
        end
        
        
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
        function dC = fReaction42Rates2002(t,C)
            
            dC = zeros(34,1);
            
            k1 = 3.1E-03;
            k2 = 3.2E+06;
            k3 = 3.1E-03;
            k4 = 2.3E+07;
            k5 = 4.4E+05;
            
            k6 = 1.3E+07;
            k7 = 2.3E+04;
            
            k8 = 1.05E+00;
            k9 = 2.5E+07;
            k10 = 6.0E+00;
            k11 = 19.0E+00;
            k12 = 2.2E+07;
            k13 = 2.4E+00;
            k14 = 1.0E+07;
            k15 = 1.8E+00;
            k16 = 7.5E+03;
            k17 = 2.0E+07;
            k18 = 5.0E-03;
            k19 = 1.0E+07;
            k20 = 1.0E-03;
            k21 = 1.0E+08;
            k22 = 8.2E+00;
            k23 = 2.2E+04;
            k24 = 6.0E-03;
            k25 = 1.0E-03;
            k26 = 2.0E+07;
            k27 = 0.2E+00;
            k28 = 4.0E+08;
            k29 = 103.0E+00;
            k30 = 1.0E+08;
            k31 = 63.5E+00;
            k32 = 1.5E+07;
            k33 = 3.6E-04;
            k34 = 9.0E+05;
            k35 = 1.1E-04;
            k36 = 3.2E+08;
            k37 = 5.0E+07;
            k38 = 1.5E+03;
            k39 = 7.1E+03;
            k40 = 4.9E+02;
            k41 = 7.1E+03;
            k42 = 2.3E+02;
            
            % From Danforth's Blood Coagulation Sensitivity paper
            k43 = 0;
            k44 = 0;
            
            dC(1,1) = -k2*C(1)*C(2) + k1*C(3) - k4*C(1)*C(4) + k3*C(5);
            dC(2,1) = -k2*C(1)*C(2) + k1*C(3) - k5*C(5)*C(2) - k6*C(6)*C(2) - k7*C(7)*C(2);
            dC(3,1) = -k1*C(3) + k2*C(1)*C(2);
            dC(4,1) = -k4*C(1)*C(4) + k3*C(5) + k5*C(5)*C(2) + k6*C(6)*C(2) + k7*C(7)*C(2);
            dC(5,1) = -k3*C(5) + k4*C(1)*C(4) - k9*C(5)*C(8) + k8*C(9)...
                -k12*C(5)*C(6) + k11*C(10) - k14*C(5)*C(11) + k13*C(12)...
                +k15*C(12) - k37*C(5)*C(27) - k42*C(5)*C(29);
            dC(6,1) = -k12*C(5)*C(6) + k11*C(10) + k22*C(18) - k28*C(6)*C(22)...
                +k27*C(23) - k34*C(6)*C(26) + k33*C(27) - k38*C(6)*C(29)...
                +k43*C(13)*C(8);
            dC(7,1) = k16*C(6)*C(14) + k32*C(25)*C(23) - k41*C(7)*C(29);
            dC(8,1) = -k9*C(5)*C(8) + k8*C(9) - k21*C(17)*C(8) + k20*C(18)...
                +k25*C(18) - k43*C(13)*C(8);
            dC(9,1) = k9*C(5)*C(8) - k10*C(9) - k8*C(9);
            dC(10,1) = k10*C(9) + k12*C(5)*C(6) - k11*C(10) - k36*C(10)*C(26)...
                +k35*C(28);
            dC(11,1) = -k14*C(5)*C(11) + k13*C(12);
            dC(12,1) = k14*C(5)*C(11) - k13*C(12) -k15*C(12);
            dC(13,1) = k15*C(12) - k19*C(16)*C(13) + k18*C(17) + k25*C(18) ...
                +k25*C(17) - k40*C(13)*C(29);
            dC(14,1) = -k16*C(6)*C(14) - k30*C(23)*C(14) + k29*C(24);
            dC(15,1) = -k17*C(7)*C(15);
            dC(16,1) = k17*C(7)*C(15) - k19*C(16)*C(13) + k18*C(17) - k24*C(16)...
                +k23*C(19)*C(20);
            dC(17,1) = k19*C(16)*C(13) - k18*C(17) - k21*C(17)*C(8) + k20*C(18)...
                +k22*C(18) - k25*C(17);
            dC(18,1) = k21*C(17)*C(8) - k20*C(18) - k22*C(18) - k25*C(18);
            dC(19,1) = k24*C(16) + k25*C(18) + k25*C(17) - k23*C(19)*C(20);
            dC(20,1) = k24*C(16) + k25*C(18) + k25*C(17) - k23*C(19)*C(20);
            dC(21,1) = -k26*C(7)*C(21) - k44*C(25)*C(21);
            dC(22,1) = k26*C(7)*C(21) - k28*C(6)*C(22) + k27*C(23) + k44*C(25)*C(21);
            dC(23,1) = k28*C(6)*C(22) - k27*C(23) - k30*C(23)*C(14) + k29*C(24)...
                +k31*C(24);
            dC(24,1) = k30*C(23)*C(14) - k29*C(24) - k31*C(24);
            dC(25,1) = k31*C(24) - k32*C(25)*C(23) - k39*C(25)*C(29);
            dC(26,1) = -k34*C(6)*C(26) + k33*C(27) - k36*C(10)*C(26) + k35*C(28);
            dC(27,1) = k34*C(6)*C(26) - k33*C(27) - k37*C(5)*C(27);
            dC(28,1) = k36*C(10)*C(26) - k35*C(28) + k37*C(5)*C(27);
            dC(29,1) = -k38*C(6)*C(29) - k39*C(25)*C(29) - k40*C(13)*C(29) ...
                -k41*C(7)*C(29) - k42*C(5)*C(29);
            dC(30,1) = k38*C(6)*C(29);
            dC(31,1) = k39*C(25)*C(29);
            dC(32,1) = k40*C(13)*C(29);
            dC(33,1) = k41*C(7)*C(29);
            dC(34,1) = k42*C(5)*C(29);
            
        end
        
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
        function dC = fReaction42Rates2002_JacEst(C)
            
            dC = zeros(34,1);
            
            k1 = 3.1E-03;
            k2 = 3.2E+06;
            k3 = 3.1E-03;
            k4 = 2.3E+07;
            k5 = 4.4E+05;
            
            k6 = 1.3E+07;
            k7 = 2.3E+04;
            
            k8 = 1.05E+00;
            k9 = 2.5E+07;
            k10 = 6.0E+00;
            k11 = 19.0E+00;
            k12 = 2.2E+07;
            k13 = 2.4E+00;
            k14 = 1.0E+07;
            k15 = 1.8E+00;
            k16 = 7.5E+03;
            k17 = 2.0E+07;
            k18 = 5.0E-03;
            k19 = 1.0E+07;
            k20 = 1.0E-03;
            k21 = 1.0E+08;
            k22 = 8.2E+00;
            k23 = 2.2E+04;
            k24 = 6.0E-03;
            k25 = 1.0E-03;
            k26 = 2.0E+07;
            k27 = 0.2E+00;
            k28 = 4.0E+08;
            k29 = 103.0E+00;
            k30 = 1.0E+08;
            k31 = 63.5E+00;
            k32 = 1.5E+07;
            k33 = 3.6E-04;
            k34 = 9.0E+05;
            k35 = 1.1E-04;
            k36 = 3.2E+08;
            k37 = 5.0E+07;
            k38 = 1.5E+03;
            k39 = 7.1E+03;
            k40 = 4.9E+02;
            k41 = 7.1E+03;
            k42 = 2.3E+02;
            
            % From Danforth's Blood Coagulation Sensitivity paper
            k43 = 0;
            k44 = 0;
            
            dC(1,1) = -k2*C(1)*C(2) + k1*C(3) - k4*C(1)*C(4) + k3*C(5);
            dC(2,1) = -k2*C(1)*C(2) + k1*C(3) - k5*C(5)*C(2) - k6*C(6)*C(2) - k7*C(7)*C(2);
            dC(3,1) = -k1*C(3) + k2*C(1)*C(2);
            dC(4,1) = -k4*C(1)*C(4) + k3*C(5) + k5*C(5)*C(2) + k6*C(6)*C(2) + k7*C(7)*C(2);
            dC(5,1) = -k3*C(5) + k4*C(1)*C(4) - k9*C(5)*C(8) + k8*C(9)...
                -k12*C(5)*C(6) + k11*C(10) - k14*C(5)*C(11) + k13*C(12)...
                +k15*C(12) - k37*C(5)*C(27) - k42*C(5)*C(29);
            dC(6,1) = -k12*C(5)*C(6) + k11*C(10) + k22*C(18) - k28*C(6)*C(22)...
                +k27*C(23) - k34*C(6)*C(26) + k33*C(27) - k38*C(6)*C(29)...
                +k43*C(13)*C(8);
            dC(7,1) = k16*C(6)*C(14) + k32*C(25)*C(23) - k41*C(7)*C(29);
            dC(8,1) = -k9*C(5)*C(8) + k8*C(9) - k21*C(17)*C(8) + k20*C(18)...
                +k25*C(18) - k43*C(13)*C(8);
            dC(9,1) = k9*C(5)*C(8) - k10*C(9) - k8*C(9);
            dC(10,1) = k10*C(9) + k12*C(5)*C(6) - k11*C(10) - k36*C(10)*C(26)...
                +k35*C(28);
            dC(11,1) = -k14*C(5)*C(11) + k13*C(12);
            dC(12,1) = k14*C(5)*C(11) - k13*C(12) -k15*C(12);
            dC(13,1) = k15*C(12) - k19*C(16)*C(13) + k18*C(17) + k25*C(18) ...
                +k25*C(17) - k40*C(13)*C(29);
            dC(14,1) = -k16*C(6)*C(14) - k30*C(23)*C(14) + k29*C(24);
            dC(15,1) = -k17*C(7)*C(15);
            dC(16,1) = k17*C(7)*C(15) - k19*C(16)*C(13) + k18*C(17) - k24*C(16)...
                +k23*C(19)*C(20);
            dC(17,1) = k19*C(16)*C(13) - k18*C(17) - k21*C(17)*C(8) + k20*C(18)...
                +k22*C(18) - k25*C(17);
            dC(18,1) = k21*C(17)*C(8) - k20*C(18) - k22*C(18) - k25*C(18);
            dC(19,1) = k24*C(16) + k25*C(18) + k25*C(17) - k23*C(19)*C(20);
            dC(20,1) = k24*C(16) + k25*C(18) + k25*C(17) - k23*C(19)*C(20);
            dC(21,1) = -k26*C(7)*C(21) - k44*C(25)*C(21);
            dC(22,1) = k26*C(7)*C(21) - k28*C(6)*C(22) + k27*C(23) + k44*C(25)*C(21);
            dC(23,1) = k28*C(6)*C(22) - k27*C(23) - k30*C(23)*C(14) + k29*C(24)...
                +k31*C(24);
            dC(24,1) = k30*C(23)*C(14) - k29*C(24) - k31*C(24);
            dC(25,1) = k31*C(24) - k32*C(25)*C(23) - k39*C(25)*C(29);
            dC(26,1) = -k34*C(6)*C(26) + k33*C(27) - k36*C(10)*C(26) + k35*C(28);
            dC(27,1) = k34*C(6)*C(26) - k33*C(27) - k37*C(5)*C(27);
            dC(28,1) = k36*C(10)*C(26) - k35*C(28) + k37*C(5)*C(27);
            dC(29,1) = -k38*C(6)*C(29) - k39*C(25)*C(29) - k40*C(13)*C(29) ...
                -k41*C(7)*C(29) - k42*C(5)*C(29);
            dC(30,1) = k38*C(6)*C(29);
            dC(31,1) = k39*C(25)*C(29);
            dC(32,1) = k40*C(13)*C(29);
            dC(33,1) = k41*C(7)*C(29);
            dC(34,1) = k42*C(5)*C(29);
            
        end
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
        function dC = fReaction44Rates2009(t,C)
            
            dC = zeros(34,1);
            
            k1 = 3.1E-03;
            k2 = 3.2E+06;
            k3 = 3.1E-03;
            k4 = 2.3E+07;
            k5 = 4.4E+05;
            
            k6 = 1.3E+07;
            k7 = 2.3E+04;
            
            k8 = 1.05E+00;
            k9 = 2.5E+07;
            k10 = 6.0E+00;
            k11 = 19.0E+00;
            k12 = 2.2E+07;
            k13 = 2.4E+00;
            k14 = 1.0E+07;
            k15 = 1.8E+00;
            k16 = 7.5E+03;
            k17 = 2.0E+07;
            k18 = 5.0E-03;
            k19 = 1.0E+07;
            k20 = 1.0E-03;
            k21 = 1.0E+08;
            k22 = 8.2E+00;
            k23 = 2.2E+04;
            k24 = 6.0E-03;
            k25 = 1.0E-03;
            k26 = 2.0E+07;
            k27 = 0.2E+00;
            k28 = 4.0E+08;
            k29 = 103.0E+00;
            k30 = 1.0E+08;
            k31 = 63.5E+00;
            k32 = 1.5E+07;
            k33 = 3.6E-04;
            k34 = 9.0E+05;
            k35 = 1.1E-04;
            k36 = 3.2E+08;
            k37 = 5.0E+07;
            k38 = 1.5E+03;
            k39 = 7.1E+03;
            k40 = 4.9E+02;
            k41 = 7.1E+03;
            k42 = 2.3E+02;
            
            % From Danforth's Blood Coagulation Sensitivity paper
            k43 = 5.7E+03;
            k44 = 3.0E+06;
            
            dC(1,1) = -k2*C(1)*C(2) + k1*C(3) - k4*C(1)*C(4) + k3*C(5);
            dC(2,1) = -k2*C(1)*C(2) + k1*C(3) - k5*C(5)*C(2) - k6*C(6)*C(2) - k7*C(7)*C(2);
            dC(3,1) = -k1*C(3) + k2*C(1)*C(2);
            dC(4,1) = -k4*C(1)*C(4) + k3*C(5) + k5*C(5)*C(2) + k6*C(6)*C(2) + k7*C(7)*C(2);
            dC(5,1) = -k3*C(5) + k4*C(1)*C(4) - k9*C(5)*C(8) + k8*C(9)...
                -k12*C(5)*C(6) + k11*C(10) - k14*C(5)*C(11) + k13*C(12)...
                +k15*C(12) - k37*C(5)*C(27) - k42*C(5)*C(29);
            dC(6,1) = -k12*C(5)*C(6) + k11*C(10) + k22*C(18) - k28*C(6)*C(22)...
                +k27*C(23) - k34*C(6)*C(26) + k33*C(27) - k38*C(6)*C(29)...
                +k43*C(13)*C(8);
            dC(7,1) = k16*C(6)*C(14) + k32*C(25)*C(23) - k41*C(7)*C(29);
            dC(8,1) = -k9*C(5)*C(8) + k8*C(9) - k21*C(17)*C(8) + k20*C(18)...
                +k25*C(18) - k43*C(13)*C(8);
            dC(9,1) = k9*C(5)*C(8) - k10*C(9) - k8*C(9);
            dC(10,1) = k10*C(9) + k12*C(5)*C(6) - k11*C(10) - k36*C(10)*C(26)...
                +k35*C(28);
            dC(11,1) = -k14*C(5)*C(11) + k13*C(12);
            dC(12,1) = k14*C(5)*C(11) - k13*C(12) -k15*C(12);
            dC(13,1) = k15*C(12) - k19*C(16)*C(13) + k18*C(17) + k25*C(18) ...
                +k25*C(17) - k40*C(13)*C(29);
            dC(14,1) = -k16*C(6)*C(14) - k30*C(23)*C(14) + k29*C(24);
            dC(15,1) = -k17*C(7)*C(15);
            dC(16,1) = k17*C(7)*C(15) - k19*C(16)*C(13) + k18*C(17) - k24*C(16)...
                +k23*C(19)*C(20);
            dC(17,1) = k19*C(16)*C(13) - k18*C(17) - k21*C(17)*C(8) + k20*C(18)...
                +k22*C(18) - k25*C(17);
            dC(18,1) = k21*C(17)*C(8) - k20*C(18) - k22*C(18) - k25*C(18);
            dC(19,1) = k24*C(16) + k25*C(18) + k25*C(17) - k23*C(19)*C(20);
            dC(20,1) = k24*C(16) + k25*C(18) + k25*C(17) - k23*C(19)*C(20);
            dC(21,1) = -k26*C(7)*C(21) - k44*C(25)*C(21);
            dC(22,1) = k26*C(7)*C(21) - k28*C(6)*C(22) + k27*C(23) + k44*C(25)*C(21);
            dC(23,1) = k28*C(6)*C(22) - k27*C(23) - k30*C(23)*C(14) + k29*C(24)...
                +k31*C(24);
            dC(24,1) = k30*C(23)*C(14) - k29*C(24) - k31*C(24);
            dC(25,1) = k31*C(24) - k32*C(25)*C(23) - k39*C(25)*C(29);
            dC(26,1) = -k34*C(6)*C(26) + k33*C(27) - k36*C(10)*C(26) + k35*C(28);
            dC(27,1) = k34*C(6)*C(26) - k33*C(27) - k37*C(5)*C(27);
            dC(28,1) = k36*C(10)*C(26) - k35*C(28) + k37*C(5)*C(27);
            dC(29,1) = -k38*C(6)*C(29) - k39*C(25)*C(29) - k40*C(13)*C(29) ...
                -k41*C(7)*C(29) - k42*C(5)*C(29);
            dC(30,1) = k38*C(6)*C(29);
            dC(31,1) = k39*C(25)*C(29);
            dC(32,1) = k40*C(13)*C(29);
            dC(33,1) = k41*C(7)*C(29);
            dC(34,1) = k42*C(5)*C(29);
            
        end
        
        
        
        function B = fUntestedGenerateMassConservationMatrix()
            % copied from sGenerateMassConservationMatrix.m
            % Mass conservation Matrix
            B = zeros(9,34);
            B(1,1) = 1; B(1,3) = 1; B(1,5) = 1; B(1,9) = 1; B(1,10) = 1;
            B(1,12) = 1; B(1,28) = 1; B(1,34) = 1;
            B(2,2) = 1; B(2,3) = 1; B(2,4) = 1; B(2,5) = 1; B(2,9) = 1;
            B(2,10) = 1; B(2,12) = 1; B(2,28) = 1; B(2,34) = 1;
            B(3,6) = 1; B(3,8) = 1; B(3,9) = 1; B(3,10) = 1; B(3,18) = 1;
            B(3,23) = 1; B(3,24) = 1; B(3,27) = 1; B(3,28) = 1; B(3,30) = 1;
            B(4,11) = 1; B(4,13) = 1; B(4,17) = 1; B(4,18) = 1; B(4,32) = 1; B(4,12) = 1;
            B(5,7) = 1; B(5,14) = 1; B(5,24) = 1; B(5,33) = 1; B(5,31) = 1; B(5,25) = 1;
            B(6,15) = 1; B(6,16) = 1; B(6,17) = 1; B(6,18) = 1; B(6,19) = 1;
            % B(6,15) = 1; B(6,16) = 1; B(6,17) = 1; B(6,18) = 1; B(6,19) = 1; B(6,20) = 1;
            B(7,21) = 1; B(7,22) = 1; B(7,23) = 1; B(7,24) = 1;
            B(8,26) = 1; B(8,27) = 1; B(8,28) = 1;
            B(9,29) = 1; B(9,30) = 1; B(9,31) = 1; B(9,32) = 1; B(9,33) = 1; B(9,34) = 1;
            
        end
        
        
        
        function [sm,pr] = sumprod(A)
            % returns the sum and product of the elements of a vector.
            
            % Other functions needs class name for calling
            sm = Tf.getsum(A);
            pr = Tf.getprod(A);
        end
        
        
        
        % functions are named normally
        function sm = getsum(A)
            sm = sum(A);
        end
        
        
        
        function pr = getprod(A)
            pr = prod(A);
        end
        
    end
    
end