# `odericcati`

Solves the Riccati differential equation for the finite-horizon linear quadratic regulator.


## Syntax

`P = odericcati(A,B,Q,R,N,'initial',P0,tspan)`\
`P = odericcati(A,B,Q,R,N,'final',PT,tspan)`


## Description
`P = odericcati(A,B,Q,R,N,'initial',P0,tspan)` solves the Riccati differential equation (<img src="https://latex.codecogs.com/svg.latex?\inline&space;\dot{\mathbf{P}}=-\left[\mathbf{A}^{T}\mathbf{P}+\mathbf{P}\mathbf{A}-(\mathbf{P}\mathbf{B}+\mathbf{N})\mathbf{R}^{-1}(\mathbf{B}^{T}\mathbf{P}+\mathbf{N}^{T})+\mathbf{Q}\right]" title="" />) for <img src="https://latex.codecogs.com/svg.latex?\inline&space;\mathbf{P}\in" title="" />.


## Examples and Additional Documentation

   - See "EXAMPLES.mlx" or the "Examples" tab on the File Exchange page for examples. 
   - See [Riccati_Differential_Equation.pdf](https://tamaskis.github.io/documentation/Riccati_Differential_Equation.pdf) (also included with download) for additional documentation.
