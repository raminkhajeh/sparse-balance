This code implements the model of:
R Khajeh, F Fumarola & LF Abbott. Sparse balance: excitatory-inhibitory networks with small bias currents and broadly distributed synaptic weights.

## Installation
This code is written in Julia (julialang.org) and tested on Julia 1.3.0.
Install all dependencies using ```julia install.jl```

## Use
To directly run ```model.jl```, inside Julia REPL, enter the command
```include("src/model.jl"); D,F = v.runsim();```  
where ```D``` and ```F``` are dictionaries, containing observables and trajectories.
To run the model using ```run.jl```, enter the command
```julia run.jl```
which calls ```model.jl```, performs a sweep over the in-degree K, and plots example responses of the networks. Data and plots are saved under ```plots/``` and ```data/```, respectively.

The variable names in the code are as consistent as possible with those that appear in the manuscript.

## Contact
Ramin Khajeh
rk2899 at Columbia dot edu
