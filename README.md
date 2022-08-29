# `solve_riccati_ode` [![View Solution of the Riccati Differential Equation (odericcati) on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/104030-solution-of-the-riccati-differential-equation-odericcati)

Solves the Riccati differential equation for the finite-horizon linear quadratic regulator.

**NOTE:** This function requires the [IVP Solver Toolbox](https://www.mathworks.com/matlabcentral/fileexchange/103975-ivp-solver-toolbox).


## Syntax

`[t,P] = solve_riccati_ode(A,B,Q,R,[],PT,tspan)`\
`[t,P] = solve_riccati_ode(A,B,Q,R,S,PT,tspan)`


## Description

`[t,P] = solve_riccati_ode(A,B,Q,R,[],PT,tspan)` solves the Riccati differential equation <img src="https://latex.codecogs.com/svg.latex?\inline&space;\dot{\mathbf{P}}=-\left[\mathbf{A}^{T}\mathbf{P}+\mathbf{P}\mathbf{A}-(\mathbf{P}\mathbf{B}+\mathbf{S})\mathbf{R}^{-1}(\mathbf{B}^{T}\mathbf{P}+\mathbf{S}^{T})+\mathbf{Q}\right]" title="" /> for <img src="https://latex.codecogs.com/svg.latex?\inline&space;\mathbf{P}(t)\in\mathbb{R}^{{n}\times{n}}" title="" />, given the state matrix <img src="https://latex.codecogs.com/svg.latex?\inline&space;\mathbf{A}\in\mathbb{R}^{{n}\times{n}}" title="" />, input matrix <img src="https://latex.codecogs.com/svg.latex?\inline&space;\mathbf{B}\in\mathbb{R}^{{n}\times{m}}" title="" />, state weighting matrix <img src="https://latex.codecogs.com/svg.latex?\inline&space;\mathbf{Q}\in\mathbb{R}^{{n}\times{n}}" title="" />, input weighting matrix <img src="https://latex.codecogs.com/svg.latex?\inline&space;\mathbf{R}\in\mathbb{R}^{{m}\times{m}}" title="" />, terminal condition <img src="https://latex.codecogs.com/svg.latex?\inline&space;\mathbf{P}_{T}\in\mathbb{R}^{{n}\times{n}}" title="" />, and the time span `tspan` over which to solve. `tspan` can be specified either as the 1×2 double `[t0,T]` where <img src="https://latex.codecogs.com/svg.latex?\inline&space;t=t_{0}" title="" /> is the initial time and <img src="https://latex.codecogs.com/svg.latex?\inline&space;t=T" title="" /> is the final time, or as a 1×(N+1) vector of times `[t0,t1,...,tNminus1,T]` at which to return the solution for <img src="https://latex.codecogs.com/svg.latex?\inline&space;\mathbf{P}(t)" title="" />. It is assumed that the cross-coupling weighting matrix is <img src="https://latex.codecogs.com/svg.latex?\inline&space;\mathbf{S}=\mathbf{0}_{{n}\times{m}}" title="" />.

`[t,P] = solve_riccati_ode(A,B,Q,R,S,PT,tspan)` does the same as the syntax above, but this time the cross-coupling weighting matrix <img src="https://latex.codecogs.com/svg.latex?\inline&space;\mathbf{S}\in\mathbb{R}^{{n}\times{m}}" title="" /> *is* specified.


## Time Vector and Solution Array

The time vector, <img src="https://latex.codecogs.com/svg.latex?\inline&space;\mathbf{t}\in\mathbb{R}^{{(N+1)}\times{1}}" title="" />, is defined as

<img src="https://latex.codecogs.com/svg.latex?\inline&space;\mathbf{t}=\begin{bmatrix}t_{0}\\\vdots\\t_{N}\end{bmatrix}" title="" />

The ith "layer" of `P` (i.e. `P(:,:,i)`) stores <img src="https://latex.codecogs.com/svg.latex?\inline&space;\mathbf{P}(t_{i})" title="" />, where <img src="https://latex.codecogs.com/svg.latex?\inline&space;t_{i}" title="" /> is the time stored in the ith element of the time vector, <img src="https://latex.codecogs.com/svg.latex?\inline&space;\mathbf{t}" title="" />.


## Examples and Additional Documentation

   - See "EXAMPLES.mlx" or the "Examples" tab on the File Exchange page for examples. 
   - See [Riccati_Differential_Equation.pdf](https://tamaskis.github.io/documentation/Riccati_Differential_Equation.pdf) (also included with download) for additional documentation.
