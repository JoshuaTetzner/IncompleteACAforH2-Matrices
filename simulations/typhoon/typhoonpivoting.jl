using BEAST
using CompScienceMeshes
using ClusterTrees
using FastBEAST
using NestedCrossApproximation
using LinearAlgebra
using Random
using DelimitedFiles

include(pwd() * "/src/pivoting.jl")

λ = 20
k = 2 * π / λ
run(`gmsh simulations/typhoon/geometry/typhoon.geo -2 -clmax 0.2  -format
msh2 -o simulations/typhoon/geometry/typhoon.msh`)
Γ = CompScienceMeshes.read_gmsh_mesh(pwd() * "/simulations/typhoon/geometry/typhoon.msh")
op = Maxwell3D.singlelayer(; wavenumber=k)
space = raviartthomas(Γ)

A = assemble(op, space, space)

Random.seed!(1)
tree = create_tree(space.pos, KMeansTreeOptions(; nmin=400))
blktree = ClusterTrees.BlockTrees.BlockTree(tree, tree)
_, fars = FastBEAST.computeinteractions(blktree; η=0.9)
fars = NestedCrossApproximation.testfars(length(tree.nodes), fars)[43]

##

idx = 43
tidcs = value(tree, idx)
sidcs = value(tree, fars)
blk = A[tidcs, sidcs]

fprows, fpcols, fperr = fullpivoting(blk)
pprows, ppcols, pperr = partialpivoting(blk)
pivstrat = IACAPivoting(space.pos[sidcs])
pivstrat = pivstrat(
    Vector(1:length(space.pos[sidcs])); ref=sum(space.pos[tidcs]) / length(tidcs)
)
mprows, mpcols, mperr = mimicrypivoting(blk, pivstrat)
##
res = hcat(Vector(1:length(fperr)), fperr)
res = hcat(res, pperr)
res = hcat(res, mperr)

##

writedlm(pwd() * "/results/typhoon/typhoonpivoting.txt", ["r\t fullpivoting\t partialpivoting\t mimicrypivoting"])
# append data
open(pwd() * "/results/typhoon/typhoonpivoting.txt", "a") do io
    writedlm(io, res)
end
