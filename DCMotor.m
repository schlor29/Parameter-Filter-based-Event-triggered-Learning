function [Ad, Bd] = DCMotor(JL)

% parameters
LS = 1.0;
dS = 0.02;
JM = 0.5;
betaM = 0.1;
R = 20;
KT = 10;
rho = 20;
kTheta = 1280.2;
JHatL = 20 * JM;
betaL = 25;
Ts = 0.1;
% uncertain parameter (load)
% JL = JHatL; % 10 * JM <= JL <= 30 * JM


A_cont = [   0                       1                   0                           0;
            -kTheta/JL      -betaL/JL   kTheta/(rho*JL)     0;
            0                       0                   0                           1;
            kTheta/(rho*JM) 0                   -kTheta/(rho^2*JM)  -(betaM+KT^2/R)/JM];


B_cont = [0;0;0;KT/R*JM];

n = size(A_cont,1);
% continuous ss
SS_cont = ss(A_cont, B_cont, eye(n), 0);

% discrete ss
sysd = c2d(SS_cont,Ts); 
Ad = sysd.A;
Bd = sysd.B;

end