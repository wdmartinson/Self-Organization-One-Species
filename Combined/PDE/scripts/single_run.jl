using DrWatson
@quickactivate "PDE"

using Base.Threads: @threads, @spawn
using BenchmarkTools

include(srcdir("Parameters.jl"));
include(srcdir("Run_zebrafish.jl"));
include(srcdir("Visualisations.jl"));


ParamsDict = Dict(
      "L₁"     => [1.5],    # x-length ## NB: the domain in the x-direction is [-L₁, L₁]
      "L₂"     => [1.5],    # y-length ## NB: the domain in the y-direction is [-L₂, L₂]
      "N₁"     => [120],    # number of gridpoints (x-axis) ## NB: This value is 0.5x the value of N_bin
      "N₂"     => [120],    # number of gridpoints (y-axis) ## NB: This value is 0.5x the value of N_bin
      "ϵ"      => [0],      # cross-diffusion strength
      "α"      => [[0.6420080713784326, 0, 0, 0]], # value found from nonlinear least squares for MelBox Migration model, 10000 runs, 240 x 240 x 150 array
      "R₀"     => [0.075],  # Radius of the ball that cells use to evaluate proliferation criteria
      "Rl"     => [0.1],    # Inner radius of annulus that cells use to evaluate proliferation criteria
      "Ru"     => [0.2],    # Outer radius of annulus that cells use to evaluate proliferation criteria
      "c₋"     => [1.0],    # Theoretical value of c₋
      "c₊"     => [7.0686], # Theoretical value of c₊ = 400.0 * pi * R₀^2
      "γ"      => [0.1221], # Growth rate for cell proliferation when N_bir = 1
      "N_bir"  => [150],    # Number of points evaluated per day for cell birth
      "ν"      => [20],     # output frequency
      "Δt"     => [0.05],   # time step
      "T"      => [150],    # final time
      "IC"     => ["MelBox"], #  "MelBox", "OffsetRectanglesMelOnly", "XanBox", "OffsetRectanglesXanOnly"
);

ParamList = dict_list(ParamsDict);
params = ParamList[1]


@time run_zebrafish(params; parallel=false)
