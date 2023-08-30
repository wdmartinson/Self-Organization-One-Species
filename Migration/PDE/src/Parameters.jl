function LoadParameters()
    ParamsDict = Dict(
    "L₁" => [1.5],   # x-length
    "L₂" => [1.5],   # y-length
    "N₁" => [120],    # number of gridpoints (x-axis)
    "N₂" => [120],    # number of gridpoints (y-axis)
    "ϵ"  => [0],     # cross-diffusion strength
    "ν"  => [20],     # output fequency
    "α"  => [0.6420080713784326, 0, 0, 0],  # scaling of amplitude of interact. potentials
    "Δt" => [0.05], # time step
    "T"  => [100],  # final time
    "IC" => ["MelBox"], # "MelBox", "CentredBlackStripe", "ABMData", "XanBox", "XanStripe", "TwoStripesMelOnly", "OffsetRectanglesMelOnly", "TwoStripesXanOnly", "OffsetRectanglesXanOnly"
    );
    return ParamsDict;
end

# New function for nonlinear least squares:
function LoadParametersWithAlpha(alphas)
    ParamsDict = Dict(
    "L₁" => 1.5,   # x-length
    "L₂" => 1.5,   # y-length
    "N₁" => 120,    # number of gridpoints (x-axis)
    "N₂" => 120,    # number of gridpoints (y-axis)
    "ϵ"  => 0,     # cross-diffusion strength
    "ν"  => 20,     # output fequency
    "α"  => alphas, # scaling of amplitude of interact. potentials
    "Δt" => 0.05, # time step
    "T"  => 150,  # final time
    "IC" => "MelBox", # "MelBox", "CentredBlackStripe", "ABMData", "XanBox", "XanStripe", "TwoStripesMelOnly", "OffsetRectanglesMelOnly", "TwoStripesXanOnly", "OffsetRectanglesXanOnly"
    );
    return ParamsDict;
end

# New function for nonlinear least squares:
function LoadParametersWithAlpha(αₘₘ, αₘₓ, αₓₘ, αₓₓ)
    ParamsDict = Dict(
    "L₁" => 1.5,   # x-length
    "L₂" => 1.5,   # y-length
    "N₁" => 120,    # number of gridpoints (x-axis)
    "N₂" => 120,    # number of gridpoints (y-axis)
    "ϵ"  => 0,     # cross-diffusion strength
    "ν"  => 20,     # output fequency
    "α"  => [αₘₘ, αₘₓ, αₓₘ, αₓₓ],# scaling of amplitude of interact. potentials
    "Δt" => 0.05, # time step
    "T"  => 150,  # final time
    "IC" => "MelBox", # "MelBox", "CentredBlackStripe", "ABMData", "XanBox", "XanStripe", "TwoStripesMelOnly", "OffsetRectanglesMelOnly", "TwoStripesXanOnly", "OffsetRectanglesXanOnly"
    );
    return ParamsDict;
end

function LoadParametersWithAlphaAndIC(alphas, IC)
    ParamsDict = Dict(
    "L₁" => 1.5,   # x-length
    "L₂" => 1.5,   # y-length
    "N₁" => 120,    # number of gridpoints (x-axis)
    "N₂" => 120,    # number of gridpoints (y-axis)
    "ϵ"  => 0,     # cross-diffusion strength
    "ν"  => 20,     # output fequency
    "α"  => alphas, # scaling of amplitude of interact. potentials
    "Δt" => 0.05, # time step
    "T"  => 150,  # final time
    "IC" => IC, # "MelBox", "CentredBlackStripe", "ABMData", "XanBox", "XanStripe", "TwoStripesMelOnly", "OffsetRectanglesMelOnly", "TwoStripesXanOnly", "OffsetRectanglesXanOnly"
    );



    return ParamsDict;
end
