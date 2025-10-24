# IncompleteACAforH2-Matrices

This repository contains all simulation scripts required to reproduce the results presented in the paper:  

**“The Incomplete Adaptive Cross Approximation for the Fast Construction of H²-Matrices and Its Application to the Electric Field Integral Equation for Electrically Small Problems”**  
*Joshua Tetzner, Simon Adrian*  

---

## 1. Overview

This repository provides the simulation scripts used to generate all figures and tables in the paper.  
All implementations are written in **Julia**, and the workflow is fully script-driven.  

At the top level, a **bash script** is provided that runs all simulations sequentially using 16 threads. The script demonstrates how each simulation can be launched individually from the command line.

---

## 2. Setup

### 2.1 Install Julia
The code was developed and tested with **Julia ≥ 1.10**.  
You can download Julia from [https://julialang.org/downloads](https://julialang.org/downloads).

### 2.2 Clone the repositories
Clone the simulation repository:
```
git clone https://github.com/JoshuaTetzner/IncompleteACAforH2-Matrices.git
```
### 2.3 Dependencies
Start a Julia REPEL and run
```julia
using Pkg
Pkg.insantiate()
```
This command installs all Julia dependencies defined in the `Project.toml` and `Manifest.toml` files.
After installation, the environment is fully configured and ready to run.

To verify the installation, you can run:
```julia
using NestedCrossApproximation
```
This should load the package without throwing any errors.

## 3. External Requirements
The geometries used in the simulations are meshed using [Gmsh](https://gmsh.info/).
Please install Gmsh and ensure that the executable gmsh is available in your system path.
You can verify this by running:
```
gmsh --version
```
Meshes are generated via the Julia package [CompScienceMeshes.jl](https://github.com/krcools/CompScienceMeshes.jl), which internally interfaces with Gmsh.
If Gmsh is not found in your path, mesh generation will fail.

On Linux, Gmsh can typically be installed via your package manager, e.g.:
```
sudo apt install gmsh
```
On Windows, download the precompiled binary from the [website](https://gmsh.info/#Download). 
