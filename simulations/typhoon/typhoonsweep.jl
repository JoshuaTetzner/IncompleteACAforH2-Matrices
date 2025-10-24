using BEAST
using CompScienceMeshes
using FastBEAST
using Random

include(pwd() * "/src/simsweep.jl")
##
tols = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8]
filename = pwd() * "/results/typhoon/typhoonsweep.txt"
results = file = open(filename, "r")
oldresults = read(file, String)
close(file)
file = open(filename, "w")
results = "tol\tlrberr_topdown\t lrberr_bottomup\tlrberr_hmat\n"
write(file, oldresults * "\n" * results)
close(file)
##
λ = 15
k = 2 * π / λ
run(`gmsh simulations/typhoon/geometry/typhoon.geo -2 -clmax 0.122 -format
msh2 -o simulations/typhoon/geometry/typhoon.msh`)
Γ = CompScienceMeshes.read_gmsh_mesh(
    pwd() * "/simulations/typhoon/geometry/typhoon.msh"
)
op = Maxwell3D.singlelayer(; wavenumber=k)
space = raviartthomas(Γ)
A = assemble(op, space, space)
for tol in tols
    Random.seed!(1)
    tree = create_tree(space.pos, KMeansTreeOptions(; nmin=250, maxlevel=100))
    @time simsweep(A, op, space, tree, filename; tol=tol, maxrank=400, η=1.0, multithreading=true)
end
