using FastBEAST
using NestedCrossApproximation
using BEAST

function simsingle(
    op,
    space,
    tree,
    filename;
    tol=1e-3,
    maxrank=50,
    multithreading=true,
    η=1.0,
    type="topdown",
)
    if type == "topdown"
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
        time = @elapsed app = PetrovGalerkinNCA(
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
    elseif type == "bottomup"
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

        time = @elapsed app = PetrovGalerkinNCA(
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
    elseif type == "H-Matrix"
        time = @elapsed app = HM.assemble(
            op,
            space,
            space;
            testtree=tree,
            trialtree=tree,
            compressor=FastBEAST.ACAOptions(; tol=tol, maxrank=maxrank),
            multithreading=multithreading,
            η=η,
        )
    end

    #---------------------------------------
    # Write data
    #---------------------------------------
    file = open(pwd() * filename, "r")
    oldresults = read(file, String)
    close(file)
    file = open(pwd() * filename, "w")
    results =
        oldresults *
        string(length(space.pos)) *
        "\t" *
        string(length(tree.levels)) *
        "\t" *
        string(time) *
        "\n"
    write(file, results)
    return close(file)
    #--------------------------------------
    # Finished data
    #--------------------------------------
end
