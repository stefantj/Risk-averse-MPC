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

There (will be) a suite of tests to verify that the code is correct.

# ToDos
1. Write a test suite
	- Verify each constraint (positive definiteness etc)
2. AFTER everything else clears, extend to constrained case
3. AFTER everything else clears, extend to robust case

