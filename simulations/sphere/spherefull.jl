using BEAST
using CompScienceMeshes
using FastBEAST
using Random

include(pwd() * "/src/simfull.jl")
##
filename = pwd() * "/results/sphere/spherefull.txt"
file = open(filename, "r")
oldresults = read(file, String)
close(file)
file = open(filename, "w")
results =
    "N\t t_topdown\t t_bottomup\t t_hmat\t s_topdown\t s_bottomup\t" *
    "s_hmat\t lrberr_topdown\t lrberr_bottomup\t lrberr_hmat\n"
write(file, oldresults * "\n" * results)
close(file)

for h in [20, 28, 40, 57, 80, 113, 160, 226]
    println("h: ", h)
    λ = 5.0
    k = 2 * π / λ
    Γ = meshicosphere(h, 1.0)

    op = Maxwell3D.singlelayer(; wavenumber=k)
    space = raviartthomas(Γ)
    Random.seed!(3)
    tree = create_tree(space.pos, KMeansTreeOptions(; nmin=100, maxlevel=100))

    @time simfull(op, space, tree, filename; tol=1e-3, maxrank=50, η=1.0, multithreading=true)
end
