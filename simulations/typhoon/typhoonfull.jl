using BEAST
using CompScienceMeshes
using FastBEAST
using Random

include(pwd() * "/src/simfull.jl")
##
filename = pwd() * "/results/typhoon/typhoonfull.txt"
file = open(filename, "r")
oldresults = read(file, String)
close(file)
file = open(filename, "w")
results =
    "N\t t_topdown\t t_bottomup\t t_hmat\t s_topdown\t s_bottomup\t" *
    "s_hmat\t lrberr_topdown\t lrberr_bottomup\t lrberr_hmat\n"
write(file, oldresults * "\n" * results)
close(file)

for h in [0.25, 0.1735, 0.122, 0.0875, 0.06, 0.0425, 0.03, 0.021, 0.015]
    λ = 15
    k = 2 * π / λ
    run(`gmsh simulations/typhoon/geometry/typhoon.geo -2 -clmax $h  -format
    msh2 -o simulations/typhoon/geometry/typhoon.msh`)
    Γ = CompScienceMeshes.read_gmsh_mesh(
        pwd() * "/simulations/typhoon/geometry/typhoon.msh"
    )

    op = Maxwell3D.singlelayer(; wavenumber=k)
    space = raviartthomas(Γ)
    Random.seed!(2)
    tree = create_tree(space.pos, KMeansTreeOptions(; nmin=100, maxlevel=100))

    @time simfull(op, space, tree, filename; tol=1e-3, maxrank=50, η=1.0, multithreading=true)
end
