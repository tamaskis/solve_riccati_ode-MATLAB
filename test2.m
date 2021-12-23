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



%% DEFINING SYSTEM

% state matrix
A = [1   1;
     2   1];

% input matrix
B = [1;
     1];

% cross-coupling weighting matrix
S = [0  0;
     0  0];

% state weighting matrix
Q = [2   1;
     1   1];

% input weighting matrix
R = 1;

% final condition
PT = [1   1;
      1   1];

% final time
T = 5;



%% SOLVING RICCATI DIFFERENTIAL EQUATION USING ODE SOLVER

[t,P] = odericcati(A,B,Q,R,[],PT,[0,T]);
P(:,:,1)
Pf = icare(A,B,Q,R)