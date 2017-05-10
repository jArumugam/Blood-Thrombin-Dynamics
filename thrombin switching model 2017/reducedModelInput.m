function dy = reducedModelInput(t, y, kValues)
% dyVec = inputFunc(t, C)
%   This function is used to compartmentalize a system of ODEs so that
%   ode23t() can handle it.

% K.surf = kValues(1);
% K.propagation = kValues(2);
% K.in1 = kValues(3);
% K.in2 = kValues(4);

K.surf = 3.6e-06;

K.propagation = kValues(1);
K.in1 = kValues(2);
K.in2 = kValues(3);

K.IIaThres = 2; % nM

Ksurft = K.surf; 

if (y(2) <= K.IIaThres)
    ksurf = Ksurft;
    kin = K.in2;
    kpropagation = 0;
    
elseif ( y(2) > K.IIaThres)
    ksurf = 0;
    kin = K.in1; 
    kpropagation = K.propagation; 
end

dy = zeros(4,1);
% [II IIa ATIII];

dy(1) = -ksurf*y(1) - kpropagation*y(1)*y(2);
dy(2) = -kin*y(2)*y(3) + ksurf*y(1) + kpropagation*y(1)*y(2);
dy(3) = -kin*y(2)*y(3);
dy(4) = kin*y(2)*y(3);



end