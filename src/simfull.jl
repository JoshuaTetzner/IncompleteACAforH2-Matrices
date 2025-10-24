using Dates
using FastBEAST
using NestedCrossApproximation

function simfull(op, space, tree, filename; tol=1e-3, maxrank=50, multithreading=true, η=1.0)

    println("This Work")
    testcomp = TopDownCompressor(IACA(space.pos), nothing)
    trialcomp = TopDownCompressor(
        IACA(
            space.pos;
            rowpivoting=MimicryPivoting(space.pos),
            columnpivoting=FastBEAST.LRF.MaximumValue(),
        ),
        nothing,
    )
    t_td = @elapsed topdown = PetrovGalerkinNCA(
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

    println("ButtomUp")
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
    t_bu = @elapsed bottomup = PetrovGalerkinNCA(
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

    t_hmat = @elapsed hmat = HM.assemble(
        op,
        space,
        space;
        testtree=tree,
        trialtree=tree,
        compressor=FastBEAST.ACAOptions(; tol=tol, maxrank=maxrank),
        multithreading=multithreading,
        η=η,
    )

    refmat = HM.assemble(
        op,
        space,
        space;
        testtree=tree,
        trialtree=tree,
        compressor=FastBEAST.ACAOptions(; tol=1e-5, maxrank=100),
        multithreading=multithreading,
        η=η,
    )

    s_td = NestedCrossApproximation.storage(topdown)
    s_bu = NestedCrossApproximation.storage(bottomup)
    s_hmat = HM.storage(hmat)

    lrb_td = NestedCrossApproximation.lrbh2mat(topdown)
    lrb_bu = NestedCrossApproximation.lrbh2mat(bottomup)
    lrb_hmat = HM.lrbhmat(hmat)
    lrbrefmat = HM.lrbhmat(refmat)
    lrberr_td = estimate_reldifference(lrb_td, lrbrefmat)
    lrberr_bu = estimate_reldifference(lrb_bu, lrbrefmat)
    lrberr_hmat = estimate_reldifference(lrb_hmat, lrbrefmat)

    #---------------------------------------
    # Write data
    #---------------------------------------
    file = open(filename, "r")
    oldresults = read(file, String)
    close(file)
    file = open(filename, "w")
    results =
        oldresults *
        string(length(space.pos)) *
        "\t" *
        string(t_td) *
        "\t" *
        string(t_bu) *
        "\t" *
        string(t_hmat) *
        "\t" *
        string(s_td[1]) *
        "\t" *
        string(s_bu[1]) *
        "\t" *
        string(s_hmat[1]) *
        "\t" *
        string(lrberr_td) *
        "\t" *
        string(lrberr_bu) *
        "\t" *
        string(lrberr_hmat) *
        "\n"
    write(file, results)
    return close(file)
    #--------------------------------------
    # Finished data
    #--------------------------------------
end
