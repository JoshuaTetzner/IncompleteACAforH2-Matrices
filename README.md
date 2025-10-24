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

## 4. Hardware Requirements
Running all simulations requires a system with at least 200 GB of RAM in the peak, depending on the problem size and number of threads used.

For reproducibility of timing results, use a single thread.
While multithreading is supported, perfect scaling cannot be expected.

## 5. Running the Simulations
### 5.1 Run all simulations
A bash script is provided at the top level to execute all simulations sequentially using 16 threads:
```bash
bash runall.sh
```
This script runs every simulation in sequence and stores the generated data in the corresponding results directory.

### 5.2 Run individual simulations
Each subdirectory inside `simulations/ contains a Julia script corresponding to a specific numerical experiment from the paper.
```julia
# Fig. 6
julia --project=. --threads=16 simulations/rectangle/rectanglepivoting.jl

# Fig. 7
julia --project=. --threads=16 simulations/typhoon/typhoonpivoting.jl

# Fig. 8
julia --project=. --threads=16 simulations/sphere/spheresweep.jl

# Fig. 9
julia --project=. --threads=16 simulations/typhoon/typhoonsweep.jl

# Fig. 10
julia --project=. --threads=16 simulations/sphere/spherefull.jl

# Fig. 11
julia --project=. --threads=16 simulations/typhoon/typhoonfull.jl

# Fig. 12
julia --project=. --threads=16 simulations/sphere/spheresweep_lf.jl

# Fig. 13
julia --project=. --threads=16 simulations/sphere/spherefull_lf.jl
```
Fig. 3 and Fig. 4 show the pivots selected in `simulations/rectangle/rectanglepivoting.jl`.
For Table I the sourcecode of the IACA in [NestedCrossApproximation.jl](https://github.com/JoshuaTetzner/NestedCrossApproximation.jl/blob/paper/IACAforH2-Matrices/src/incompletefactorization/incompleteaca.jl) must me modified to track the time required for the pivoting and for the entiere IACA. 

## 7. Reproducibility Statement

The source code that support the findings of this paper is openly available at [NestedCrossApproximaiton.jl](https://github.com/JoshuaTetzner/NestedCrossApproximation.jl/tree/paper/IACAforH2-Matrices) and Zenodo.
This is a priliminary version of the NCA implementations, the main branch of the package might have been updated in the meantime. 




