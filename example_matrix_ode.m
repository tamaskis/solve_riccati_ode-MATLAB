%% example_matrix_ode.m
% ODE Solver Toolbox
%
% Example for solving a matrix-valued ODE (the Riccati differential
% equation).
%
% Copyright © 2021 Tamas Kis
% Last Update: 2021-12-22
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com



%% SCRIPT SETUP

% clears Workspace and Command Window, closes all figures
clear; clc; close all;

% adds path to all "ODE Solver Toolbox" functions
addpath(genpath('..'));



%% DEFINING SYSTEM

% state matrix
A = [1   1;
     2   1];

% input matrix
B = [1;
     1];

% state weighting matrix
Q = [2   1;
     1   1];

% input weighting matrix
R = 1;

% cross-coupling weighting matrix
S = [0  0;
     0  0];

% final condition
PT = [1   1;
      1   1];

% final time
T = 5;



%% SOLVING RICCATI DIFFERENTIAL EQUATION USING ODE SOLVER

% defines the Riccati differential equation (a matrix-valued ODE)
F = @(t,P) -(A.'*P+P*A-(P*B+S)/R*(B.'*P+S.')+Q);

% converts the matrix-valued ODE to a vector-valued ODE
f = odefun_mat2vec(F);

% final condition
yT = odeIC_mat2vec(PT);

% solves vector-valued ODE using a step size of h = 0.001
[t,y] = RK4(f,[T,0],yT,0.001);

% transforms solution matrix for vector-valued ODE into solution array for
% matrix-valued ODE
P = odesol_vec2mat(y);

% solution for P0 (will be at end of array since P solved for backwards in
% time)
P0 = P(:,:,end)



% %% SOLVING RICCATI DIFFERENTIAL EQUATION USING ONE-STEP PROPAGATION
% 
% % time vector between t = 5 and t = 0 with a spacing of h = 0.001.
% h = -0.001;
% t = (5:h:0)';
% 
% % preallocate vector to store solution
% P = zeros(2,2,length(t));
% 
% % store initial condition
% P(:,:,1) = PT;
% 
% % solving using "RK4_step"
% for i = 1:(length(t)-1)
%     P(:,:,i+1) = RK4_step(F,t(i),P(:,:,i),h);
% end
% 
% % solution for P0 using one-step propagation
% P0_step = P(:,:,end);
% 
% % maximum absolute error between the two results (should be 0)
% max(abs(P0-P0_step),[],'all')


%% test
Pinf = icare(A,B,Q,R);
P_norm = zeros(size(t));
for i = 1:length(t)
    P_norm(i) = norm(P(:,:,i),'fro');
end

figure;
hold on;
plot(t,P_norm,'LineWidth',1.5);
plot(t,norm(Pinf,'fro')*ones(size(t)),'k--','LineWidth',1.5);
hold off;
grid on;
xlabel('$t$','interpreter','latex','fontsize',18);
ylabel('$\|\mathbf{P}\|_{\mathrm{F}}$','interpreter','latex','fontsize',18);
legend('$\mathbf{P}(t)$','$\mathbf{P}_{\infty}$','interpreter','latex','fontsize',14,'location','southeast');
hold off;