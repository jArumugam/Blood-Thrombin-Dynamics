function [sumErr] = reducedModelError(kValues, C1, C2, C3, C4, Cinit)
% Enter C1 C2 C3 C4 in this order
% II IIa AT IIaAT
% Jayavel Arumugam
% [sumErr] = reducedModelError(kS, kA, kI)
% 11/8/16

C0 = zeros(4, 1);
C0(1, 1) = Cinit(1);
C0(2, 1) = Cinit(2);
C0(3, 1) = Cinit(3);
C0(4, 1) = Cinit(4);
tRange = 0:1:1200; % s

%% Computation
fun = @(t,y) reducedModelInput(t, y, kValues);
% options = odeset('AbsTol', 1e-3*ones(1,4));
% [~, C] = ode23t(fun, tRange, C0, options);
[~, C] = ode23s(fun, tRange, C0);

s1 = (C1*1e09 - C(:,1))/C0(1);
s2 = (C2*1e09 - C(:,2))/C0(1);
s3 = (C3*1e09 - C(:,3))/C0(3);
s4 = (C4*1e09 - C(:,4))/C0(3);

mTimePoints = length(C1); 
sumErr = 1/mTimePoints*sum( s1.^2 ) ...
        + 1/mTimePoints*sum( s2.^2 ) ...
        + 1/mTimePoints*sum( s3.^2 ) ... 
        + 1/mTimePoints*sum( s4.^2 ); 

% sumErr = 1/mTimePoints*sqrt(sum( s1.^2 )) ...
%         + 1/mTimePoints*sqrt(sum( s2.^2 )) ...
%         + 1/mTimePoints*sqrt(sum( s3.^2 )) ... 
%         + 1/mTimePoints*sqrt(sum( s4.^2 )); 

% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% % This was used in thesis
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% % Normalized Error
% errC1 = norm(C1*1e09 - C(:,1))/C0(1);
% errC2 = norm(C2*1e09 - C(:,2))/C0(1);
% errC3 = norm(C3*1e09 - C(:,3))/C0(3);
% errC4 = norm(C4*1e09 - C(:,4))/C0(3);
% % % 
% sumErr = errC1 + errC2 + errC3 + errC4;


end

