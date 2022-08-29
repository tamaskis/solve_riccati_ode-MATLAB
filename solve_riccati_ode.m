%==========================================================================
%
% solve_riccati_ode  Solves the Riccati differential equation for the 
% finite-horizon linear quadratic regulator.
%
%   [t,P] = solve_riccati_ode(A,B,Q,R,[],PT,tspan)
%   [t,P] = solve_riccati_ode(A,B,Q,R,S,PT,tspan)
%
% Copyright © 2021 Tamas Kis
% Last Update: 2022-08-28
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com
%
% TECHNICAL DOCUMENTATION:
% https://tamaskis.github.io/files/Riccati_Differential_Equation.pdf
%
% DEPENDENCIES:
%   • IVP Solver Toolbox (https://www.mathworks.com/matlabcentral/fileexchange/103975-ivp-solver-toolbox)
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   A       - (n×n double) state/system/dynamics matrix
%   B       - (n×m double) input/control matrix
%   Q       - (n×n double) state weighting matrix
%   R       - (m×m double) input/control weighting matrix
%   S       - (OPTIONAL) (n×m double) cross-coupling weighting matrix
%               --> defaults to 0
%   PT      - (n×n double) terminal condition at time t = T
%   tspan   - (1×2 or 1×(N+1) double) time interval to solve over
%               --> To specify the interval only and allow MATLAB to choose
%                   the intermediate times, use the 1×2 vector [t0 T]
%               --> To specify specific times you want the solution at, use
%                   the 1×(N+1) vector: [t0 t1 ... tN]
%
% -------
% OUTPUT:
% -------
%   t       - ((N+1)×1 double) time vector
%   P       - (n×n×(N+1) double) solution of Riccati differential equation
%
% -----------
% CONDITIONS:
% -----------
%   1. M ⪰ 0 (positive semidefinite). If S = 0, this conditions reduces to
%      the following two conditions:
%       (a) Q ⪰ 0 (positive semidefinite)
%       (b) R ≻ 0 (positive definite)
%   2. PT ⪰ 0 (positive semidefinite)
%   3. (A,B) stabilizable
%   4. (A-BR^(-1)S^T,Q-SR^(-1)S^T) detectable
%       • if S = 0, this condition reduces to (A,Q^(1/2)) detectable
%
%==========================================================================
function [t,P] = solve_riccati_ode(A,B,Q,R,S,PT,tspan)
    
    % ----------------------
    % Determines dimensions.
    % ----------------------
    
    % state dimension
    n = size(A,1);
    
    % input dimension
    m = size(B,2);
    
    % ----------------------------------------------------
    % Sets unspecified parameters to their default values.
    % ----------------------------------------------------
    
    % defaults cross-coupling matrix to 0
    if isempty(S)
        S = zeros(n,m);
    end
    
    % ---------
    % Solution.
    % ---------
    
    % flips tspan so the Riccati ODE is solved backwards in time
    tspan = fliplr(tspan);
    
    % defines Riccati ODE has a matrix-valued ODE
    dPdt = @(t,P) -(A.'*P+P*A-(P*B+S)/R*(B.'*P+S.')+Q);
    
    % converts matrix-valued ODE to vector-valued ODE
    dydt = mat2vec_ode(dPdt);
    
    % converts matrix initial condition to vector initial condition
    yT = mat2vec_IC(PT);
    
    % solves Riccati ODE
    [t,y] = ode45(dydt,tspan,yT);
    
    % transforms solution matrix for vector-valued ODE into solution array
    % for matrix-valued ODE
    P = vec2mat_sol(y);
    
    % reorders t so that time is increasing
    t = flipud(t);
    
    % reorders solution for P
    P_reordered = zeros(size(P));
    for k = 1:length(t)
        P_reordered(:,:,k) = P(:,:,length(t)-k+1);
    end
    P = P_reordered;
    
end