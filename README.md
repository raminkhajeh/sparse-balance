This code implements the model of:
R Khajeh, F Fumarola & LF Abbott. Sparse balance: excitatory-inhibitory networks with small bias currents and broadly distributed synaptic weights. PLOS Comp Bio (2022).

REQUIREMENTS:
This code is written in Julia (julialang.org) and tested on Julia 1.3.0. Required packages are listed at the top of the code (Using [package name]).
## installation
install all dependencies using ```julia install.jl```

USE:
To run the model, within Julia REPL, enter the command:
>>> include("model.jl"); D,F = v.runsim(); 
where D and F are dictionaries, containing observables and trajectories.
The notation used in the code is similar to that in the paper.

CONTACT:
Ramin Khajeh
rk2899@columbia.edu

