using DrWatson
@quickactivate "PDE"
include(srcdir("Parameters.jl"));
include(srcdir("Run_zebrafish.jl"));
include(srcdir("Visualisations.jl"));


ParamsDict = Dict(
    "L₁" => [1.5],   # x-length ## NB: the domain in the x-direction is [-L₁, L₁]
    "L₂" => [1.5],   # y-length ## NB: the domain in the y-direction is [-L₂, L₂]
    "N₁" => [120],    # number of gridpoints (x-axis) ## NB: This value is 0.5x the value of N_bin
    "N₂" => [120],    # number of gridpoints (y-axis) ## NB: This value is 0.5x the value of N_bin
    "ϵ"  => [0],     # cross-diffusion strength
    "ν"  => [20],     # output frequency

    ## Unscaled potentials
    # "α"  => [[1.0, 0, 0, 0]],  # unscaled melanophore potential
    # "α"  => [[0, 0, 0, 1.0]],  # unscaled xanthophore potential

    ## MelStripe (a.k.a CentredBlackStripe) Initial Condition
    # "α"  => [[0.3993934507468565, 0, 0, 0]], # value found from nonlinear least squares for CentredBlackStripe, 1000 runs, 60 x 60 x 150 array
    # "α"  => [[0.5816611520542333, 0, 0, 0]], # value found from nonlinear least squares for CentredBlackStripe, 1000 runs, 120 x 120 x 150 array
    # "α"  => [[0.6331931583429926, 0, 0, 0]], # value found from nonlinear least squares for CentredBlackStripe, 1000 runs, 240 x 240 x 150 array
    # "α"  => [[, 0, 0, 0]], #
    # "α"  => [[0.637318325759627, 0, 0, 0]], # value found from nonlinear least squares for CentredBlackStripe, 1000 runs, 240 x 240 x 150 array, FIT-TO-60X60-MESH-DATA

    ## MelBox Initial Condition
    # "α"  => [[0.40771925127808584, 0, 0, 0]], # value found from nonlinear least squares for MelBox, 10000 runs, 60 x 60 x 150 array
    # "α"  => [[0.5932776267571058, 0, 0, 0]], # value found from nonlinear least squares for MelBox, 10000 runs, 120 x 120 x 150 array
    "α"  => [[0.6420080713784326, 0, 0, 0]], # value found from nonlinear least squares for MelBox, 10000 runs, 240 x 240 x 150 array
    # "α"  => [[0.6574173594493339, 0, 0, 0]], # value found from nonlinear least squares for MelBox, 10000 runs, 480 x 480 x 150 array
    # "α"  => [[, 0, 0, 0]], #
    # "α"  => [[0.6543938884792129, 0, 0, 0]], # value found from nonlinear least squares for MelBox, 10000 runs, 240 x 240 x 150 array, FIT-TO-60X60-MESH-DATA
    # "α"  => [[0.39665231948564755, 0, 0, 0]], # value found from L∞ metric for MelBox support, 10000 runs, 240 x 240 x 150 array

    ## OffsetRectanglesMelOnly Initial Condition
    # "α"  => [[0.38055781650597237, 0, 0, 0]], # value found from nonlinear least squares for OffsetRectanglesMelOnly, 1000 runs, 60 x 60 x 150 array
    # "α"  => [[0.5587302670916284, 0, 0, 0]], # value found from nonlinear least squares for OffsetRectanglesMelOnly, 1000 runs, 120 x 120 x 150 array
    # "α"  => [[0.6040842079291687, 0, 0, 0]], # value found from nonlinear least squares for OffsetRectanglesMelOnly, 1000 runs, 240 x 240 x 150 array
    # "α"  => [[, 0, 0, 0]], #
    # "α"  => [[0.617959619336824, 0, 0, 0]], # value found from nonlinear least squares for OffsetRectanglesMelOnly, 1000 runs, 240 x 240 x 150 array, FIT-TO-60X60-MESH-DATA

    ## TwoStripesMelOnly Initial Condition
    # "α"  => [[0.39965357219023745, 0, 0, 0]], # value found from nonlinear least squares for TwoStripesMelOnly, 1000 runs, 60 x 60 x 150 array
    # "α"  => [[0.6252370291710823, 0, 0, 0]], # value found from nonlinear least squares for TwoStripesMelOnly, 1000 runs, 120 x 120 x 150 array
    # "α"  => [[0.6386999960117437, 0, 0, 0]], # value found from nonlinear least squares for TwoStripesMelOnly, 1000 runs, 240 x 240 x 150 array
    # "α"  => [[, 0, 0, 0]], #
    # "α"  => [[0.6543470799941368, 0, 0, 0]], # value found from nonlinear least squares for TwoStripesMelOnly, 1000 runs, 240 x 240 x 150 array, FIT-TO-60X60-MESH-DATA




    ## XanBox Initial Condition
    # "α"  => [[0, 0, 0, 0.13611520908912572]], # value found from nonlinear least squares for XanBox, 5000 runs, 60 x 60 x 150 array
    # "α"  => [[0, 0, 0, 0.35023053731672094]], # value found from nonlinear least squares for XanBox, 5000 runs, 120 x 120 x 150 array
    # "α"  => [[0, 0, 0, 0.476620798374944]], # value found from nonlinear least squares for XanBox, 5000 runs, 240 x 240 x 150 array
    # "α"  => [[0, 0, 0, 0.5055427719256165]], # value found from nonlinear least squares for XanBox, 5000 runs, 480 x 480 x 150 array
    # "α"  => [[, 0, 0, 0]], #
    # "α"  => [[0, 0, 0, 0.5054483806532654]], # value found from nonlinear least squares for XanBox, 5000 runs, 240 x 240 x 150 array, FIT-TO-60X60-MESH-DATA
    # "α"  => [[0, 0, 0, 0.2748972776308001]], # value found from L∞ metric for XanBox support, 5000 runs, 240 x 240 x 150 array

    ## XanStripe Initial Condition
    # "α"  => [[0, 0, 0, 0.14464598602408615]], # value found from nonlinear least squares for XanStripe, 5000 runs, 60 x 60 x 150 array
    # "α"  => [[0, 0, 0, 0.37888261801295753]], # value found from nonlinear least squares for XanStripe, 5000 runs, 120 x 120 x 150 array
    # "α"  => [[0, 0, 0, 0.5024141240648899]], # value found from nonlinear least squares for XanStripe, 5000 runs, 240 x 240 x 150 array
    # "α"  => [[0, 0, 0, ]], #
    # "α"  => [[0, 0, 0, 0.5263402011472122]], # value found from nonlinear least squares for XanStripe, 5000 runs, 240 x 240 x 150 array, FIT-TO-60X60-MESH-DATA

    ## TwoStripesXanOnly Initial Condition
    # "α"  => [[0, 0, 0, 0.15171662829282737]], # value found from nonlinear least squares for TwoStripesXanOnly, 1000 runs, 60 x 60 x 150 array
    # "α"  => [[0, 0, 0, 0.3731606365354493]], # value found from nonlinear least squares for TwoStripesXanOnly, 1000 runs, 120 x 120 x 150 array
    # "α"  => [[0, 0, 0, 0.501454101499117]], # value found from nonlinear least squares for TwoStripesXanOnly, 1000 runs, 240 x 240 x 150 array
    # "α"  => [[0, 0, 0, ]], #
    # "α"  => [[0, 0, 0, 0.5342028699730562]], # value found from nonlinear least squares for TwoStripesXanOnly, 1000 runs, 240 x 240 x 150 array, FIT-TO-60X60-MESH-DATA

    ## OffsetRectanglesXanOnly Initial Condition
    # "α"  => [[0, 0, 0, 0.1383150325886836]], # value found from nonlinear least squares for OffsetRectanglesXanOnly, 1000 runs, 60 x 60 x 150 array
    # "α"  => [[0, 0, 0, 0.35557192336448035]], # value found from nonlinear least squares for OffsetRectanglesXanOnly, 1000 runs, 120 x 120 x 150 array
    # "α"  => [[0, 0, 0, 0.48200836751174053]], # value found from nonlinear least squares for OffsetRectanglesXanOnly, 1000 runs, 240 x 240 x 150 array
    # "α"  => [[0, 0, 0, ]], #
    # "α"  => [[0, 0, 0, 0.5130293928685054]], # value found from nonlinear least squares for OffsetRectanglesXanOnly, 1000 runs, 240 x 240 x 150 array, FIT-TO-60X60-MESH-DATA

    "Δt" => [0.05], # time step
    "T"  => [150],  # final time
    "IC" => ["MelBox"] #  "MelBox", "CentredBlackStripe", "ABMData", "XanBox", "XanStripe", "TwoStripesMelOnly", "OffsetRectanglesMelOnly", "TwoStripesXanOnly", "OffsetRectanglesXanOnly"
);
ParamList = dict_list(ParamsDict);
params = ParamList[1];

make_melano_movie_or_run(params)
