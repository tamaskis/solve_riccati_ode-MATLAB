%==========================================================================
%
% odericcati  Solves the Riccati differential equation for the 
% finite-horizon linear quadratic regulator.
%
%   P = odericcati(A,B,Q,R,N,'initial',P0,tspan)
%   P = odericcati(A,B,Q,R,N,'final',PT,tspan)
%
% Copyright © 2021 Tamas Kis
% Last Update: 2021-12-23
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com
%
% TECHNICAL DOCUMENTATION:
% https://tamaskis.github.io/documentation/Riccati_Differential_Equation.pdf
%
% REFERENCES:
%   [1] https://www.mathworks.com/matlabcentral/answers/94722-how-can-i-solve-the-matrix-riccati-differential-equation-within-matlab
%   [2] https://en.wikipedia.org/wiki/Linear-quadratic_regulator
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
%   N       - (OPTIONAL) (n×m double) cross-coupling weighting matrix
%               --> defaults to 0
%   type    - (char) 'final' or 'initial'; describes condition type on P
%   Pc      - (n×n double) initial condition (i.e. P0) if type = 'initial', 
%             final condition (i.e. PT) if type = 'final'
%   tspan   - (1×2 or 1×(N+1) double) time interval to solve over
%               --> if input as the 1×2 vector [t0 T], TODO
%               --> if input as the 1×(N+1) vector [t0 t1 ... tN], TODO
%
% -------
% OUTPUT:
% -------
%   t       - ((N+1)×1 double) time vector
%   P       - (n×n×(N+1) double) solution of Riccati differential equation
%
% -----
% NOTE:
% -----
%   --> N+1 = length of time vector
%
%==========================================================================
function [t,P] = odericcati(A,B,Q,R,N,type,Pc,tspan)
    
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
    if (nargin < 5) || isempty(N)
        N = zeros(n,m);
    end

    % ---------
    % Solution.
    % ---------

    % flips tspan if given final condition (i.e. solve backwards in time)
    if strcmpi(type,'final')
        tspan = fliplr(tspan);
    end

    % defines Riccati ODE has a matrix-valued ODE
    dPdt = @(t,P) -(A.'*P+P*A-(P*B+N)/R*(B.'*P+N.')+Q);

    % converts matrix-valued ODE to vector-valued ODE
    dydt = odefun_mat2vec(dPdt);

    % converts matrix IC to vector IC
    y0 = odeIC_mat2vec(Pc);
    
    % solves Riccati ODE
    [t,y] = ode45(dydt,tspan,y0);

    % transforms solution matrix for vector-valued ODE into solution array
    % for matrix-valued ODE
    P = odesol_vec2mat(y);

    % reorders t and P if solved backwards in time
    if strcmpi(type,'final')

        % reorders t
        t = flipud(t);

        % reorders solution for P
        P_reordered = zeros(size(P));
        for k = 1:length(t)
            P_reordered(:,:,k) = P(:,:,length(t)-k+1);
        end
        P = P_reordered;

    end
    
end