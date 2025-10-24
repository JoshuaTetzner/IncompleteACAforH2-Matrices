julia --project=. --threads=16 simulations/sphere/spherefull.jl
julia --project=. --threads=16 simulations/sphere/spherefull_lf.jl
julia --project=. --threads=16 simulations/sphere/spheresweep.jl
julia --project=. --threads=16 simulations/sphere/spheresweep_lf.jl
julia --project=. --threads=16 simulations/typhoon/typhoonfull.jl
julia --project=. --threads=16 simulations/typhoon/typhoonsweep.jl
julia --project=. --threads=16 simulations/typhoon/typhoonpivoting.jl
julia --project=. --threads=16 simulations/rectangle/rectanglepivoting.jl