using BEAST
using CompScienceMeshes
using FastBEAST
using Random

include(pwd() * "/src/simsweep.jl")
##

tols = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10]
filename = pwd() * "/results/sphere/spheresweep.txt"
results = file = open(filename, "r")
oldresults = read(file, String)
close(file)
file = open(filename, "w")
results = "tol\tlrberr_topdown\t lrberr_bottomup\tlrberr_hmat\n"
write(file, oldresults * "\n" * results)
close(file)


λ = 2.0
k = 2 * π / λ
Γ = meshicosphere(28, 1.0)
op = Maxwell3D.singlelayer(; wavenumber=k)
space = raviartthomas(Γ)
A = assemble(op, space, space)
for tol in tols
    Random.seed!(3)
    tree = create_tree(space.pos, KMeansTreeOptions(; nmin=200, maxlevel=100))
    @time simsweep(A, op, space, tree, filename; tol=tol, maxrank=200, η=1.0, multithreading=true)
end
