%==========================================================================
%
% odericcati  Solves the Riccati differential equation.
%
%   P = odericcati(A,B,Q,R,N,'initial',P0,tspan)
%   P = odericcati(A,B,Q,R,N,'final',PT,tspan)
%
% Copyright © 2021 Tamas Kis
% Last Update: 2021-11-10
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com
%
% TECHNICAL DOCUMENTATION:
% https://tamaskis.github.io/documentation/Riccati_Differential_Equation.pdf
%
% REFERENCES:
%   [1] https://www.mathworks.com/matlabcentral/answers/94722-how-can-i-solve-the-matrix-riccati-differential-equation-within-matlab
%   [2] https://en.wikipedia.org/wiki/Linear%E2%80%93quadratic_regulator#Finite-horizon,_continuous-time_LQR
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   A       - (n×n double) system/dynamics matrix
%   B       - (n×m double) input/control matrix
%   Q       - (n×n double) state weighting matrix
%   R       - (m×m double) control weighting matrix
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
    
    % ---------------------
    % Determines dimension.
    % ---------------------

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
    
    % solves Riccati ODE
    [t,P] = ode45(@(t,P)riccati(t,P,A,B,Q,R,N),tspan,Pc);

    % reshape (N+1)×n^2 solution matrix to n×n×(N+1)
    P_reshaped = zeros(n,n,length(t));
    for i = 1:length(t)
        P_reshaped(:,:,i) = reshape(P(i,:),size(A));
    end
    P = P_reshaped;

    % flips t if solved backwards in time
    if strcmpi(type,'final')
        t = flipud(t);
    end

    % reorders solution for P if solved backwards in time
    P_reordered = zeros(size(P));
    for i = 1:length(t)
        P_reordered(:,:,i) = P(:,:,length(t)-i+1);
    end
    P = P_reordered;

    %----------------------------------------------------------------------
    % riccati
    %
    % Riccati differential equation.
    %----------------------------------------------------------------------
    %
    % INPUT:
    %   t       - (1×1 double) time (needed for ODE solver convention)
    %   P       - (n^2×1 double) column vectors of P(t) stacked into a
    %             single vector
    %   A       - (n×n double) system/dynamics matrix
    %   B       - (n×m double) input/control matrix
    %   Q       - (n×n double) state weighting matrix
    %   R       - (m×m double) control weighting matrix
    %   N       - (n×m double) cross-coupling weighting matrix
    %
    % OUTPUT:
    %   dPdt    - (n^2×1 double) column vectors of dP/dt at time t stacked
    %             into a single vector
    %
    %----------------------------------------------------------------------
    function dPdt = riccati(t,P,A,B,Q,R,N)

        % reshapes n^2×1 vector into n×n matrix
        P = reshape(P,size(A));
        
        % evaluates Riccati differential equation
        dPdt = -(A.'*P+P*A-(P*B+N)/R*(B.'*P+N.')+Q);
        
        % reshapes n×n matrix into n^2×1 vector
        dPdt = dPdt(:);

    end
    
end