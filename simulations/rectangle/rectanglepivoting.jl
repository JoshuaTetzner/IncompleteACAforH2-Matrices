using BEAST
using CompScienceMeshes
using NestedCrossApproximation
using LinearAlgebra
using StaticArrays
using DelimitedFiles

include(pwd() * "/src/pivoting.jl")

λ = 10
k = 2 * π / λ
Γt = meshrectangle(1.0, 1.0, 0.1)
Γs = CompScienceMeshes.translate(meshrectangle(2.0, 4.5, 0.1), SVector(2.5, 0.0, 0.0))
Γs2 = CompScienceMeshes.translate(meshrectangle(2.5, 2.0, 0.1), SVector(0.0, 2.5, 0.0))
Γs = weld(Γs, Γs2)
op = Maxwell3D.singlelayer(; wavenumber=k)
t = raviartthomas(Γt)
Ft = raviartthomas(Γs)
A = assemble(op, t, Ft)

fprows, fpcols, fperr = fullpivoting(A)
pprows, ppcols, pperr = partialpivoting(A)
pivstrat = IACAPivoting(Ft.pos)
pivstrat = pivstrat(Vector(1:length(Ft.pos)); ref=SVector(0.5, 0.5, 0.5))
mprows, mpcols, mperr = mimicrypivoting(A, pivstrat)

##

res = hcat(Vector(1:length(fperr)), fperr)
res = hcat(res, pperr)
res = hcat(res, mperr)

##

writedlm(pwd() * "/results/rectangle/rectanglepivoting.txt", ["r\t fullpivoting\t partialpivoting\t mimicrypivoting"])
open(pwd() * "/results/rectangle/rectanglepivoting.txt", "a") do io
    writedlm(io, res)
end
