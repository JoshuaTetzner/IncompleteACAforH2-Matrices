using FastBEAST
using NestedCrossApproximation
using ClusterTrees
using LinearAlgebra

function simsweep(
    A, op, space, tree, filename; tol=1e-3, maxrank=200, multithreading=true, η=1.0
)
    println("top-down NCA")
    testcomp = TopDownCompressor(IACA(space.pos), nothing)
    trialcomp = TopDownCompressor(
        IACA(
            space.pos;
            rowpivoting=MimicryPivoting(space.pos),
            columnpivoting=FastBEAST.LRF.MaximumValue(),
        ),
        nothing,
    )
    topdown = PetrovGalerkinNCA(
        op,
        space,
        space;
        testtree=tree,
        trialtree=tree,
        testcompressor=testcomp,
        trialcompressor=trialcomp,
        maxrank=maxrank,
        tol=tol,
        multithreading=multithreading,
        η=η,
    )

    println("bottom-up NCA")
    testcomp = BottomUpCompressor(
        IACA(
            LRF.MaximumValue(Bool[]),
            TreeMimicryPivoting(tree, space.pos),
            IncompleteNormEstimator(Float64[], Float64(0.0)),
        ),
        nothing,
    )
    trialcomp = BottomUpCompressor(
        IACA(
            TreeMimicryPivoting(tree, space.pos),
            LRF.MaximumValue(Bool[]),
            IncompleteNormEstimator(Float64[], Float64(0.0)),
        ),
        nothing,
    )
    bottomup = PetrovGalerkinNCA(
        op,
        space,
        space;
        testtree=tree,
        trialtree=tree,
        testcompressor=testcomp,
        trialcompressor=trialcomp,
        maxrank=maxrank,
        multithreading=multithreading,
        tol=tol,
        η=η,
    )
    println("H-Matrix")
    hmat = HM.assemble(
        op,
        space,
        space;
        testtree=tree,
        trialtree=tree,
        compressor=FastBEAST.ACAOptions(; tol=tol, maxrank=maxrank),
        multithreading=multithreading,
        η=η,
    )

    lrb_td = NestedCrossApproximation.lrbmat(topdown)
    lrb_bu = NestedCrossApproximation.lrbmat(bottomup)
    lrb_hmat = HM.lrbmat(hmat)
    lrb_fm = NestedCrossApproximation.lrbmat(A, topdown)

    rel_lrb_td = norm(lrb_td - lrb_fm) / norm(lrb_fm)
    rel_lrb_bu = norm(lrb_bu - lrb_fm) / norm(lrb_fm)
    rel_lrb_hmat = norm(lrb_hmat - lrb_fm) / norm(lrb_fm)
    #---------------------------------------
    # Write data
    #---------------------------------------
    file = open(filename, "r")
    oldresults = read(file, String)
    close(file)
    file = open(filename, "w")
    results =
        oldresults *
        string(tol) *
        "\t" *
        string(rel_lrb_td) *
        "\t" *
        string(rel_lrb_bu) *
        "\t" *
        string(rel_lrb_hmat) *
        "\n"
    write(file, results)
    return close(file)
    #--------------------------------------
    # Finished data
    #--------------------------------------
end
