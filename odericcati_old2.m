%==========================================================================
%
% odericcati  Solves the Riccati differential equation.
%
%   P = odericcati(A,B,Q,R,N,'initial',P0,tspan)
%   P = odericcati(A,B,Q,R,N,'final',PT,tspan)
%
% Copyright © 2021 Tamas Kis
% Last Update: 2021-12-13
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
function [t,P] = odericcati_old2(A,B,Q,R,N,type,Pc,tspan)
    
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

    % defines Riccati ODE has a matrix-valued ODE
    dPdt = @(t,P) -(A.'*P+P*A-(P*B+N)/R*(B.'*P+N.')+Q);

    % converts matrix-valued ODE to vector-valued ODE
    dydt = @(t,y) odefun_mat2vec(dPdt,t,y);

    % converts matrix IC to vector IC
    y0 = Pc(:);
    
    % solves Riccati ODE
    [t,y] = ode45(dydt,tspan,y0);

    % TODO
    P = odesol_vec2mat(y);

    % flips t if solved backwards in time
    if strcmpi(type,'final')
        t = flipud(t);
    end

    % reorders solution for P if solved backwards in time
    P_reordered = zeros(size(P));
    for k = 1:length(t)
        P_reordered(:,:,k) = P(:,:,length(t)-k+1);
    end
    P = P_reordered;



    function ydot = odefun_mat2vec(F,t,y,p)
    
        % determine state dimension if not input (assuming a Y is square)
        if nargin < 4
            p = sqrt(length(y));
        end
    
        % determines q, where Y is a p×q matrix
        q = length(y)/p;
        
        % reshapes pq×1 state vector into p×q state matrix
        Y = reshape(y,[p,q]);
    
        % evaluates matrix ODE
        Ydot = F(t,Y);
        
        % reshapes p×q state matrix derivative into pq×1 state vector deriv.
        ydot = Ydot(:);
        
    end

    function Y = odesol_vec2mat(y,p)
    
        % state vector dimension
        pq = size(y,2);
    
        % determine state dimension if not input (assuming a Y is square)
        if nargin < 2
            p = sqrt(pq);
        end
    
        % determines q, where Y is a p×q matrix
        q = pq/p;
        
        % determines N, where the ODE solution is given at N+1 points in time
        N = size(y,1)-1;
        
        % preallocate array to store time history of state matrix
        Y = zeros(p,q,N+1);
    
        % populates state matrix
        for i = 1:(N+1)
            Y(:,:,i) = reshape(y(i,:),[p,q]);
        end
        
    end
    
end