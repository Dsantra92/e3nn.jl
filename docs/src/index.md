# E3nn: E(3) Equivariant Neural Networks

[E3nn.jl](https://github.com/Dsantra92/e3nn.jl) provides a pure julia implementation of the [e3nn framework](https://e3nn.org/).
It aims to provide fast and extendable APIs for working with $\mathbb{E}^3$ Equivariant Neural Networks and other performing related operations.
It is built on top of [GraphNeuralNetworks.jl]() and [Flux.jl]().
The APIs are generally consistent with e3nn's [PyTorch](https://github.com/e3nn/e3nn/) and [JAX](https://github.com/e3nn/e3nn-jax/) libraries but there might be certain differences.

## Installation

The package isn't available through the Julia Package Manager (yet).
If you want to use it, you can add it directly from GitHub.

```julia-repl
julia> using Pkg

julia> Pkg.add("https://github.com/Dsantra92/e3nn.jl.git")
```

!!! note
The package is still in under early stages of development and the APIs might change in the future.

## Getting started

We highly recommend starting with [Introduction to Irreps]() to gain some understanding of the primary components of the network.
If you are looking for a bit more hands on experience, have a look at the [examples]().
