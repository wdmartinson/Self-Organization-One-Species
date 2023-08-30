using PyCall

function initialise(XMID, YMID; type = "MelBox")
    # compute grid size params.
    Δx = XMID[1, 2] - XMID[1, 1]
    Δy = YMID[2, 1] - YMID[1, 1]

    # uniform distribution of all cell types
    if type == "MelBox"
        Mel = (XMID .< 0.5) .* (XMID .> -0.5) .* (YMID .> -0.5) .* (YMID .< 0.5)
        Mel = 400 * Mel / (Δx * Δy * sum(Mel)) # with melanophore density of 400 cells/mm^2
        Xan = zeros(size(YMID))

        return Mel, Xan


    elseif type == "CentredBlackStripe"
        Mel = (YMID .< 0.25) .* (YMID .> -0.25) # This is where the error is: for certain mesh sizes, fractions of voxels can be covered by cells and this function does not capture that
        Mel = 600 * Mel / (Δx * Δy * sum(Mel)) # with number of melanophores equal to 600 (density = 400 cells/mm^2)
        Xan = zeros(size(YMID))

        return Mel, Xan


    elseif type == "ABMData"
        sio = pyimport("scipy.io");
        ABM_Data = sio.loadmat("../../ABM/data/sims/MC/ABM_Histograms/XXXXX.mat");
        Mel = ABM_Data["Chr"]; # only care about the matrix (Array{Float64,3} type)
        Mel = Mel[:,:,1]';
        Xan = zeros(size(Mel)); # can include the xanthophore data from the ABM similar to the above

        return Mel, Xan


    elseif type == "XanBox"
        Xan = (XMID .< 0.5) .* (XMID .> -0.5) .* (YMID .> -0.5) .* (YMID .< 0.5)
        Xan = 772 * Xan / (Δx * Δy * sum(Xan)) # with xanthophore density of 772 cells/mm^2 (number density)
        Mel = zeros(size(YMID))

        return Mel, Xan


    elseif type == "XanStripe"
        Xan = (YMID .< 0.25) .* (YMID .> -0.25) # This is where the error is: for certain mesh sizes, fractions of voxels can be covered by cells and this function does not capture that
        Xan = 772 * 0.5 * 3* Xan / (Δx * Δy * sum(Xan)) # with Xanthophore density of 772 cells/mm^2 (number density)
        Mel = zeros(size(YMID))

        return Mel, Xan


    elseif type == "TwoStripesMelOnly"
        Mel1 = (YMID .< 0.5) .* (YMID .> 0)
        Mel1 = 600 * Mel1 / (Δx * Δy * sum(Mel1)) # with number of melanophores equal to 600 (density = 400 cells/mm^2)

        Mel2 = (YMID .< -0.082) .* (YMID .> -(0.082 + 0.5))
        Mel2 = 600 * Mel2 / (Δx * Δy * sum(Mel2)) # with number of xanthophores equal to 834 (density = 772 cells/mm^2)

        Mel = Mel1 .+ Mel2;
        Xan = zeros(size(YMID))

        return Mel, Xan


    elseif type == "OffsetRectanglesMelOnly"
        Mel1 = (XMID .< -0.1) .* (XMID .> -0.6) .* (YMID .> -0.75) .* (YMID .< 0.25)
        Mel1 = 200 * Mel1 / (Δx * Δy * sum(Mel1))

        Mel2 = (XMID .< 0.6) .* (XMID .> 0.1) .* (YMID .> -0.25) .* (YMID .< 0.75)
        Mel2 = 200 * Mel2 / (Δx * Δy * sum(Mel2))
        Mel  = Mel1 .+ Mel2;

        Xan  = zeros(size(YMID));

        return Mel, Xan


    elseif type == "TwoStripesXanOnly"
        Mel  = zeros(size(YMID));

        Xan1 = (YMID .< 0.36) .* (YMID .> 0)
        Xan1 = 834 * Xan1 / (Δx * Δy * sum(Xan1)) # with number of melanophores equal to 600 (density = 400 cells/mm^2)

        Xan2 = (YMID .< -0.082) .* (YMID .> -(0.082 + 0.36))
        Xan2 = 834 * Xan2 / (Δx * Δy * sum(Xan2)) # with number of xanthophores equal to 834 (density = 772 cells/mm^2)
        Xan  = Xan1 .+ Xan2;

        return Mel, Xan


    elseif type == "OffsetRectanglesXanOnly"
        Mel  = zeros(size(YMID));

        Xan1 = (XMID .< -0.1) .* (XMID .> -0.6) .* (YMID .> -0.75) .* (YMID .< 0.25)
        Xan1 = 772/2 * Xan1 / (Δx * Δy * sum(Xan1))

        Xan2 = (XMID .< 0.6) .* (XMID .> 0.1) .* (YMID .> -0.25) .* (YMID .< 0.75)
        Xan2 = 772/2 * Xan2 / (Δx * Δy * sum(Xan2))
        Xan  = Xan1 .+ Xan2;

        return Mel, Xan
    end
end

function PlaceHorizontalStripe(lower, upper, XMID, YMID, type)
    k, l = size(YMID)

    Δx = XMID[1, 2] - XMID[1, 1]
    Δy = YMID[2, 1] - YMID[1, 1]

    temp = (YMID .< upper) .* (YMID .> lower);
    # Without noise:
    temp = max.(temp + temp .* 0 .* (rand(k, l) .- 0.5), 0)
    # With noise:
    # temp = max.(temp + temp .* 0.1 * (rand(k, l) .- 0.5), 0)

    if type == "Mel"
        temp = 3 * 400 * (upper - lower) * temp / (Δx * Δy * sum(temp)) ;
    else
        temp = 3 * 772 * (upper - lower) * temp / (Δx * Δy * sum(temp)) ;
    end

    return temp
end

function PlaceBox(Lx, Ly, Ux, Uy, XMID, YMID, type)
    k, l = size(YMID)

    Δx = XMID[1, 2] - XMID[1, 1];
    Δy = YMID[2, 1] - YMID[1, 1];

    temp = (YMID .< Uy) .* (YMID .> Ly) .* (XMID .< Ux) .* (XMID .> Lx);

    #temp = max.(temp + temp .* 1 * (rand(k, l) .- 0.5), 0)

    if type == "Mel"
        temp = 400 * (Uy - Ly) * (Ux - Lx) * temp / (Δx * Δy * sum(temp)) ;
    else
        temp = 772 * (Uy - Ly) * (Ux - Lx) * temp / (Δx * Δy * sum(temp)) ;
    end

    return temp
end
