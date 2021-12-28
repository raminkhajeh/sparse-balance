This code implements the model of:
R Khajeh, F Fumarola & LF Abbott. Sparse balance: excitatory-inhibitory networks with small bias currents and broadly distributed synaptic weights.

## Installation
This code is written in Julia (julialang.org) and tested on Julia 1.3.0.
Install all dependencies using ```julia install.jl```

## Use
To run the model, within Julia REPL, enter the command:
```include("model.jl"); D,F = v.runsim();``` 
where D and F are dictionaries, containing observables and trajectories.
The variable names in the code is consistent with those that appear in the manuscript.

## Contact
Ramin Khajeh
