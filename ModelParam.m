
%**************************************************************************
%     The Intensive Control Insulin-Nutrition-Glucose (ICING) Model
%     Lin, J. et al (2011)
%     Last modified: 31/05/16
%     Author: Jose Fernando Garcia and Estefania Aguirre Zapata
%**************************************************************************

%%
%Symbolic

syms pG Si aG VG k1 k2 k3 d1 d2 aI VI EGPb CNS nI nC nK nL xL Pmax Qsp Isp BGsp Gi P1sp P2sp uexsp Dsp PNsp uensp Psp beta1 beta2

%%
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

%%
P=d2*P2sp+PNsp ;
Psp = 0.7672;
uensp = k1*exp(-Isp^(k2/k3));
           % Assumed by lack of info
        % No dextrose is added intravenously


%%
%Definitions 

%Functions 
F1= -pG*BGsp - Si*BGsp*Qsp/(1 + aG*Qsp) + (P + EGPb - CNS)/VG;
F2=beta1*BGsp-beta2*Gi;
F3=nI*(Isp - Qsp) - nC*Qsp/(1 + aG*Qsp);
F4=-nK*Isp - nL*Isp/(1 + aI*Isp)- nI*(Isp - Qsp) + uexsp/VI + (1 - xL)*uensp/VI;
F5=-d1*P1sp+Dsp;
F6=-d2*P2sp+d1*P1sp;


%Operating points
%xss=[5;20;10.7655;22;111];
%uss=[7.5933;5;0];
xss=[5 5 10.7655 20 22 111];
uss=[7.5933 0];
wss=0.7672;
yss=5;
xssn=[xss wss];

%System definition
system=[F1;F2;F3;F4;F5;F6];


%States and inputs definition 
states=[BGsp Gi Qsp Isp P1sp P2sp];
inputs= [uexsp PNsp];
dist=Dsp;


%%
%Linealization 

%Jacobian Matrix
Ac=jacobian(system,states);
Bc3= jacobian(system,inputs);
Bc1= jacobian(system,dist);
Cc=[0 1 0 0 0 0];
Dc=0;

% Output Matrix
C=Cc;
D=Dc;

%Replace the operating points
A=subs(Ac,[states inputs],[xss uss]);
A=double(A);

B3= subs(Bc3,[states inputs dist],[xss uss wss]);
B3=double(B3); %Matriz asociada a la variable manipulada

B1= subs(Bc1,[states inputs dist],[xss uss wss]);
B1=double(B1); %Matriz asociada a la perturbaci?n desconocida 
B=[B3 B1];


%%
%DISCRETIZACI?N - Tiempo de muestre de 1segundo

%State Space
Sysc= ss(A,B,C,D);

%Discretizaci?n 
Sysd= c2d(Sysc,1); 

%Asignaci?n 
Ad= Sysd.a;
Bd= Sysd.b;
Cd= Sysd.c;
Dd= Sysd.d;

Bd3=Bd(:,1);
Bd1=Bd(:,3);


global aux;
aux=[];

global auxfiltro;
auxfiltro=[];

%--------------------------------------------------------------------------
%Matriz A ampliada 
%-------------------------------------------------------------------------


% An=[Ad Bd1; 0 0 0 0 0 0 1];
% Bn=[Bd3;0];
% Cn=[0 1 0 0 0 0];

Ann=[Ac Bc1; 0 0 0 0 0 0 1];
Bnn=[Bc3(:,1);0];
Cnn=[0 1 0 0 0 0 0];


O=obsv(Ad,Cd);
Ob=rank(O);
C=ctrb(Ad,Bd);
Con=rank(C);
