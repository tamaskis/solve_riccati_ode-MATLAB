clear;clc;close all;

% state space system
A = [ 0   1   0   0;
      0   0   1   0;
     -3   1   2   3;
      2   1   0   0];
B = [0   0;
     0   0;
     1   2;
     0   2];

n = size(A,1);
m = size(B,2);

Q = eye(n);
R = eye(m);
N = zeros(n,m);
PT = ones(n);

%P = odericcati(A,B,Q,R,[],P0);

tspan = [0 10];

[t,P] = odericcati(A,B,Q,R,[],'final',PT,tspan);
P