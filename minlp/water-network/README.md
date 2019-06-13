# Optimizing water network with Knitro

This folder implements some methods to optimize flows on a water-network
with Knitro. The problem formulates naturally as a Non Linear Program (NLP).
We aim at studying the impact of arcs removals by using the MINLP solver
implemented in Knitro.

The tutorial was first developed for the first
[Journée Julia et optimisation](https://julialang.univ-nantes.fr/journee-julia-et-optimisation/),
held June 17th, 2019 at Nantes University. The author thanks the organizers
for having invited him.


## Contents

The tutorial implements:

- `data.jl`: define a small water network;
- `utils.jl`: some plotting utilities;
- `network.jl`: the tutorial itself.

An extension for large-scale networks is implemented with:
- `large_scale.jl`: a script to build a large-scale network;
- `large_network.jl`: optimize the large-scale network with Knitro and JuMP.

A HTML version of the tutorial can be built locally by using Weave.jl:

```julia
using Weave
weave("network.jl")

```

## Installation

We know focus on installing the tutorial from scratch. We recommend
to install Julia `1.1` before.


### Installing Knitro

First, we need to install Knitro locally. Knitro is proprietary software
and require a license for use.

1- Create an account on the [Artelys website](https://www.artelys.com/fr/espace-client/).
2- Download Knitro for your specific platform (Windows, Mac or Linux).
3- Install Knitro on your machine.
4- Run `get_machine_ID` script to get your machine ID.
4- Download a valid license corresponding to your machine ID (licenses are provided for the tutorial).

### Setting up the tutorial

The tutorial runs with Julia 1.1. It uses `PyPlot` to plot the different
results. Note that even if `PyPlot` is not required, it is recommended
to install it to visualize the result.

Activate the tutorial environment inside Julia with:

```julia
] activate .
] instantiate

```

