using BEAST
using CompScienceMeshes
using FastBEAST
using Random

include(pwd() * "/src/simfull.jl")
##
filename = pwd() * "/results/sphere/spherefull_lf.txt"

#T_Φ
file = open(filename, "r")
oldresults = read(file, String)
close(file)
file = open(filename, "w")
results =
    "T_Φ\nN\t t_topdown\t t_bottomup\t t_hmat\t s_topdown\t s_bottomup\t" *
    "s_hmat\t lrberr_topdown\t lrberr_bottomup\t lrberr_hmat\n"
write(file, oldresults * "\n" * results)
close(file)

for h in [20, 28, 40, 57, 80, 113, 160]
    println("h: ", h)
    λ = 5.0
    k = 2 * π / λ
    Γ = meshicosphere(h, 1.0)
    gamma = im * k
    op = Maxwell3D.singlelayer(; gamma=gamma, alpha=0.0, beta=1.0)
    space = raviartthomas(Γ)
    Random.seed!(3)
    tree = create_tree(space.pos, KMeansTreeOptions(; nmin=100, maxlevel=100))

    @time simfull(op, space, tree, filename; tol=1e-3, maxrank=50, η=1.0, multithreading=true)
end

#T_A
file = open(filename, "r")
oldresults = read(file, String)
close(file)
file = open(filename, "w")
results =
    "T_A\nN\t t_topdown\t t_bottomup\t t_hmat\t s_topdown\t s_bottomup\t" *
    "s_hmat\t lrberr_topdown\t lrberr_bottomup\t lrberr_hmat\n"
write(file, oldresults * "\n" * results)
close(file)

for h in [20, 28, 40, 57, 80, 113, 160]
    println("h: ", h)
    λ = 5.0
    k = 2 * π / λ
    Γ = meshicosphere(h, 1.0)
    gamma = im * k
    op = Maxwell3D.singlelayer(; gamma=gamma, alpha=1.0, beta=0.0)
    space = raviartthomas(Γ)
    Random.seed!(3)
    tree = create_tree(space.pos, KMeansTreeOptions(; nmin=100, maxlevel=100))

    @time simfull(op, space, tree, filename; tol=1e-3, maxrank=50, η=1.0, multithreading=true)
end
