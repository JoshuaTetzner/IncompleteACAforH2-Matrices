function fullpivoting(M::Matrix{F}; maxiter=rank(M) - 2) where {F}
    R = copy(M)
    norms = Float64[]
    rows = Int[]
    cols = Int[]
    push!(norms, 1.0)

    for i in 1:maxiter
        ij = argmax(abs.(R))
        R -= R[:, ij[2]] * R[ij[1], ij[2]]^-1 * transpose(R[ij[1], :])
        push!(rows, ij[1])
        push!(cols, ij[2])
        push!(norms, norm(R) / norm(M))
    end
    return rows, cols, norms
end

function partialpivoting(R::Matrix{K}; maxiter=rank(R) - 2) where {K}
    M = copy(R)
    err = Float64[]
    usedr = zeros(Bool, size(M, 1))
    usedc = zeros(Bool, size(M, 2))
    refnorm = norm(M)
    push!(err, 1.0)
    c1 = 1
    usedc[c1] = true
    r1 = argmax(abs.(M[:, c1]) .* (.!usedr))
    usedr[r1] = true
    c = argmax(abs.(M[r1, :]) .* (.!usedc))
    usedc[c] = true
    M -= M[:, c1] * M[r1, c1]^-1 * transpose(M[r1, :])
    push!(err, norm(M))

    for i in 1:(maxiter-1)
        r = argmax(abs.(M[:, c]) .* (.!usedr))
        usedr[r] = true
        cn = argmax(abs.(M[r, :]) .* (.!usedc))
        usedc[cn] = true
        M -= M[:, c] * M[r, c]^-1 * transpose(M[r, :])
        c = cn

        push!(err, norm(M) / refnorm)
    end
    rows = findall(x -> x, usedr)
    cols = findall(x -> x, usedc)
    return rows, cols, err
end

function mimicrypivoting(A::Matrix{K}, pivstrat; maxiter=rank(A) - 2) where {K}
    rows = Int[]
    cols = Int[]
    err = [1.0]
    Ac = copy(A)
    nextcol = pivstrat()
    nextrow = argmax(abs.(A[:, nextcol]))
    push!(rows, nextrow)
    push!(cols, nextcol)
    for i in 1:maxiter
        Ac =
            Ac -
            Ac[:, cols[end]] * 1 / Ac[rows[end], cols[end]] * transpose(Ac[rows[end], :])
        push!(err, norm(Ac) / norm(A))
        nextcol = pivstrat(i + 1)
        nextrow = argmax(abs.(Ac[:, nextcol]))
        push!(rows, nextrow)
        push!(cols, nextcol)
    end
    return rows, cols, err
end
