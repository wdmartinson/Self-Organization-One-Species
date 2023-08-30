function initialise(XMID, YMID; type = "uniform")
    # compute grid size params.
    Δx = XMID[1, 2] - XMID[1, 1]
    Δy = YMID[2, 1] - YMID[1, 1]

    # uniform distribution of all cell types
    if type == "MelBox"
        Mel = (XMID .< 0.5) .* (XMID .> -0.5) .* (YMID .> -0.5) .* (YMID .< 0.5)
        Mel = 400 * Mel / (Δx * Δy * sum(Mel)) # with melanophore density of 400 cells/mm^2
        Xan = zeros(size(YMID))

        return Mel, Xan

    elseif type == "XanBox"
        Xan = (XMID .< 0.5) .* (XMID .> -0.5) .* (YMID .> -0.5) .* (YMID .< 0.5)
        Xan = 772 * Xan / (Δx * Δy * sum(Xan)) # with xanthophore density of 772 cells/mm^2 (number density)
        Mel = zeros(size(YMID))

        return Mel, Xan


    elseif type == "OffsetRectanglesMelOnly"
        Mel1 = (XMID .< -0.1) .* (XMID .> -0.6) .* (YMID .> -0.75) .* (YMID .< 0.25)
        Mel1 = 200 * Mel1 / (Δx * Δy * sum(Mel1))

        Mel2 = (XMID .< 0.6) .* (XMID .> 0.1) .* (YMID .> -0.25) .* (YMID .< 0.75)
        Mel2 = 200 * Mel2 / (Δx * Δy * sum(Mel2))
        Mel  = Mel1 .+ Mel2;

        Xan  = zeros(size(YMID));

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
