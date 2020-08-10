function xdot = modelNL_filtro(x, u)
% This function is used to generate states 
    % State vector -> x = [Gp Gi Q I P1 P2 D]
 
    %u = C_input(t);
    
    % Parameters
    pG = 0.006;         % Patient endogenous glucose removal [min^-1]
    Si = 0.0002;         % From former model since it is identified online [L/mU/min]
    aG = 0.0154;        % Saturation of insulin-stimulated glucose [L/mU]
    EGPb = 1.16;        % Basal endogenous glucose production [mmol/min]
    CNS = 0.3;          % Insulin independent central nervous system glucose uptake [mmol/min]
    VG = 13.3;          % Glucose distribution volume [L]
    VI = 3.15;          % Insulin distribution volume [L]
    aI = 0.0017;        % Saturation of plasma insulin disappearance [L/mU]
    nI = 0.003;         % Transcapillary difussion rate [min^-1]
    nC = nI;            % [min^-1]
    nK = 0.0542;        % Kidney clearance [min^-1]
    nL = 0.1578;        % Patient specific liver clearance [min^-1]
    xL = 0.67;          % First pass endogenous insulin hepatic uptake []
    d1 = 0.0347;        % Transport rate [min^-1]
    d2 = 0.0069;        % Transport rate [min^-1]
    Pmax = 6.11;        % Saturation value of P2 [mmol/min]
    k1 = 45.7;          % Base rate for endogenous insulin production [mU/min]
    k2 = 1.5;           % Generic constant for exponential supression []
    k3 = 1000;          % Generic constant for exponential supression []
    beta1=0.09;
    beta2=beta1;

    P = min(d2*x(6),Pmax) + u(2);
    uen = k1*exp(-x(4)^(k2/k3));
    
    xdot = zeros(6, 1);
        
    xdot(1) = -pG.*x(1) - Si.*x(1).*x(3)./(1 + aG.*x(3)) + (P + EGPb - CNS)./VG;
    xdot(2) = beta1.*x(1)-beta2.*x(2);
    xdot(3) = nI.*(x(4) - x(3)) - nC.*x(3)./(1 + aG.*x(3));
    xdot(4) = -nK.*x(4) - nL.*x(4)./(1 + aI.*x(4)) - nI.*(x(4) - x(3)) + u(1)./VI + (1 - xL).*uen./VI;
    xdot(5) = -d1.*x(5) + u(3);
    xdot(6) = -min(d2.*x(6),Pmax) + d1.*x(5);
    
end