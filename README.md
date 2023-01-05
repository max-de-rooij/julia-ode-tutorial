# ODE Modelling and Parameter Estimation with Julia
This repository contains sample files explaining the basics of Julia, required to get started with ODE modelling, as well as some more advanced tutorials on simulating systems of ordinary differential equations, performing sensitivity analysis, parameter estimation and parameter identifiability analysis. In short, this repository contains all basic components needed to get started with using Julia.

## Getting Started with Julia
To start using Julia, you can download the latest version from [julialang.org](https://julialang.org) to be able to start using it. To be able to use it in VS Code, you also require the Julia VS Code extension.

### The REPL
The command line in Julia is called the REPL (Read-Eval-Print-Loop). When opening Julia, this is what launches. From here, we can execute julia code, install packages, manage projects or run code files. The output of julia code will be printed, unless ended with a semi-colon (;). For example, the following code will print the value '3'

```julia
julia> 1 + 2
```


### The Julia package manager (Pkg)
Installing packages happens through `Pkg`, the julia package manager. This also includes functionality for environments. When opening Julia, we can instantiate the package manager by pressing "`]`". We can press backspace to return to the regular Julia REPL. After activating `Pkg`, you can install a package by typing `add <package name>`. The convenient plotting package `Plots.jl` can be installed using

```julia
pkg> add Plots
```

This will install the Plots library in the global environment. To create a new environment, use

```julia
pkg> activate <.path/to/working/directory>
```

Once this has been activated, you can add packages in a similar way. When you start adding packages, notice that a `Project.toml` file and sometimes a `Manifest.toml` file will be created. These files contain the project dependencies and allow others to easily set up an environment with the same project dependencies as yours using
```julia
pkg> activate <.path/to/working/directory>
pkg> instantiate
```

In a similar way, you can use the `Project.toml` and `Manifest.toml` files in this project to set up the right environment for running all code. After activating the right environment, you can get started with the first notebook: [1-the-basics.ipynb](https://github.com/max-de-rooij/julia-ode-tutorial/blob/d3f9e86674ce46118cdf536ac63714f45bceb2f9/1-the-basics.ipynb).