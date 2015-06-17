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

| First Header  | Second Header |
| ------------- | ------------- |
| Content Cell  | Content Cell  |
| Content Cell  | Content Cell  |
