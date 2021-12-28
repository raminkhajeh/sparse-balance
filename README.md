This code implements the model of:
R Khajeh, F Fumarola & LF Abbott. Sparse balance: excitatory-inhibitory networks with small bias currents and broadly distributed synaptic weights.

## Installation
This code is written in Julia (julialang.org) and tested on Julia 1.3.0.
Install all dependencies using ```julia install.jl```
We have attempted to keep variable names consistent with those that appear in the manuscript. 

## Use:
To run the model, within Julia REPL, enter the command:
```include("model.jl"); D,F = v.runsim();``` 
where D and F are dictionaries, containing observables and trajectories.
The notation used in the code is similar to that in the paper.

## Contact:
Ramin Khajeh
