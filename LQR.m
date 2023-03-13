function K = LQR(Q,R,mod)
% discrete-time LQR

% Discrete-time Algebraic Riccati Equation
[P,~,~] = dare(mod.A,mod.B,Q,R);

% Discrete-time LQR controller
K = -1*inv(R + mod.B.'*P*mod.B)*mod.B.'*P*mod.A;
end