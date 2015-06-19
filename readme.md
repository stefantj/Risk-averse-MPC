# Getting started

This library uses the julia language (http://julialang.org/), with the packages Convex, SCS, and Debug. To add these packages, use 
```
   Pkg.add("Convex")
   Pkg.add("SCS")
   Pkg.add("Debug")
```

To run the simulation, run julia in the directory with these files and run
```
   include("MPC_sim.jl")
```

The program will output a csv file with fields

|  Status data                   | State/action data                             |  Performance data     |
| ------------------------------ | --------------------------------------------- | --------------------- |
| iteration, feasible, rand_ind, | X(1), ... , X(n_state), u(1), ... ,u(n_input),| r, c_i, gamma_i, E(c) |

This data can be easily interpreted by running the plotter.m file in Matlab.

# Verification of code
Each function that calls the solver has an optional final argument that is a test-mode flag.

When the flag is set, then the function will print out diagnostics on each of the constraints.

When testing positive (semi-)definiteness of a matrix, the minimum eigen value is returned. This should be within tolerance (typically 1e-5) of zero.

When testing state dynamics, the 2-norm of the error is printed.
=======
# ToDos
1. Extend to constrained case
2. Extend to robust case
3. Extend to N=3,4,5
4. Add different risk-polytopes
