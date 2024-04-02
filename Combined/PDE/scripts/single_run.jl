# cd("/Users/markusschmidtchen/Dropbox/Zebrafish_MicroMacroModelling/Code/Combined/Macro") # Markus's path
# cd("/Users/guest91/Dropbox/Zebrafish_MicroMacroModelling/Code/ImageStats")# Markus uni



using DrWatson
@quickactivate "Macro"

using Base.Threads: @threads, @spawn
using BenchmarkTools

include(srcdir("Parameters.jl"));
include(srcdir("Run_zebrafish.jl"));
include(srcdir("Visualisations.jl"));


# ParamsDict = Dict(
#       "L₁"     => [1.5],   # x-length
#       "L₂"     => [1.5],   # y-length
#       "N₁"     => [120],    # number of gridpoints (x-axis)
#       "N₂"     => [120],    # number of gridpoints (y-axis)
#       "ϵ"      => [0],     # cross-diffusion strength
#       # "α"      => [[0.6276, 0.6401440037117885, 0.4655416234321407, 0.4809]],
#       "α"      => [[0.6420080713784326, 0, 0, 0]], # value found from nonlinear least squares for MelBox Migration model, 10000 runs, 240 x 240 x 150 array
#       # "c₊"     => [], # Maximum value of
#       # "γ"      => [],
#       # "N_dart" => [],
#       "ν"      => [20],     # output frequency
#       "Δt"     => [0.05], # time step
#       "T"      => [150],  # final time
#       "IC"     => ["MelBox"], # "uniform", "MelBox", "CentredBlackStripe", "TwoStripes", "UniformPerturbed", "ABMData", "XanBox", "XanStripe", "OffsetRectangles", "ReplicationExptUniformPerturbed", "ReplicationExptPerturbedStripes", "UniformPerturbedWholeDomainMX", "CheckerboardBio", "TwoStripesMelOnly", "OffsetRectanglesMelOnly", "TwoStripesXanOnly", "OffsetRectanglesXanOnly"
# );

# ParamList = dict_list(ParamsDict);
# params = ParamList[1]

# The input α is vector-valued, so set
# [αₘₘ, γₘ, c₊] = [1.906055716924618, 0.0663099011024432, 6.628139504174085] # Value found from minimising L2 Error from MelBox IC
# [αₘₘ, γₘ, c₊] = [1.8438393660228698, 0.06838572421250995, 6.632025591319713] # Value found from minimising L2 Error from OffsetRectanglesMelOnly IC
αₘₘ = 1.906055716924618;
αₘₓ = 0;
αₓₘ = 0;
αₓₓ = 0;
γₘ  = 0.0663099011024432;
c₊  = 6.628139504174085;
Ndart = 150;
# αs = [p[1], 0, 0, 0]; # MM, MX, XM, XX
αs = [αₘₘ, αₘₓ, αₓₘ, αₓₓ];
# Set the value of α for the PDE params variable:
params = LoadCombinedModelParameters(αs, γₘ, 0.0, c₊, 2.0, Ndart);

@time run_zebrafish(params; parallel=false)
# make_melano_movie_or_run(params; type="state")





# plot_initial_data(params)
# plot_timeseries(params)
# make_melano_movie_or_run(params)
# plot_potential(params, "MelMel")

# plot_initial_data(params)
