using PyCall

function initialise(XMID, YMID; type = "uniform")
    # compute grid size params.
    Δx = XMID[1, 2] - XMID[1, 1]
    Δy = YMID[2, 1] - YMID[1, 1]

    # uniform distribution of all cell types
    if type == "uniform"
        XanD = 200 * (XMID .< 1) .* (XMID .> -1) .* (YMID .> -1) .* (YMID .< 1)
        Mel = 200 * (XMID .< 1)  .* (XMID .> -1) .* (YMID .> -1) .* (YMID .< 1)
        IriD = 200 * (XMID .< 1) .* (XMID .> -1) .* (YMID .> -1) .* (YMID .< 1)
        IriL = 200 * (XMID .< 1) .* (XMID .> -1) .* (YMID .> -1) .* (YMID .< 1)
        return XanD, Mel, IriD, IriL

        # random fluctuations in the initial cell densities
    elseif type == "MelBox"
        Mel = (XMID .< 0.5) .* (XMID .> -0.5) .* (YMID .> -0.5) .* (YMID .< 0.5)
        Mel = 400 * Mel / (Δx * Δy * sum(Mel)) # with melanophore density of 400 cells/mm^2
        # Mel = Mel / (Δx * Δy * sum(Mel)) # normalise to 1
        Xan = zeros(size(YMID))

        return Mel, Xan
    elseif type == "CentredBlackStripe"
        Mel = (YMID .< 0.25) .* (YMID .> -0.25) # This is where the error is: for certain mesh sizes, fractions of voxels can be covered by cells and this function does not capture that
        # Mel = 400 * Mel / (Δx * Δy * sum(Mel)) # with melanophore density of 400cells/mm^2
        # Mel = Mel / (Δx * Δy * sum(Mel)) # normalise to 1.
        Mel = 600 * Mel / (Δx * Δy * sum(Mel)) # with number of melanophores equal to 600 (density = 400 cells/mm^2)
        Xan = zeros(size(YMID))

        return Mel, Xan
    elseif type == "TwoStripes"
        Mel = (YMID .< 0.5) .* (YMID .> 0)
        # Mel = 400 * Mel / (Δx * Δy * sum(Mel)) # with melanophore density of 400cells/mm^2
        Mel = 600 * Mel / (Δx * Δy * sum(Mel)) # with number of melanophores equal to 600 (density = 400 cells/mm^2)

        Xan = (YMID .< -0.082) .* (YMID .> -(0.082 + 0.36))
        # Xan = 772 * Xan / (Δx * Δy * sum(Xan)) # with xanthophore density of 772cells/mm^2
        Xan = 834 * Xan / (Δx * Δy * sum(Xan)) # with number of xanthophores equal to 834 (density = 772 cells/mm^2)
        return Mel, Xan
    elseif type == "UniformPerturbed"
        Mel = ones(size(XMID))
        k, l = size(XMID)
        Mel = Mel + 0.5 * (rand(k, l) .- 0.5)
        Mel = Mel .* rand(k, l)
        Mel = 400 * Mel / (Δx * Δy * sum(Mel))

        Xan = ones(size(XMID))
        Xan = Xan + 0.5 * (rand(k, l) .- 0.5)
        Xan = Xan .* rand(k, l)
        Xan = 400 * Xan / (Δx * Δy * sum(Xan))

        return Mel, Xan
    elseif type == "UniformPerturbed"
        Mel = (XMID .< 1) .* (XMID .> 0.5) .* (YMID .< -0.5) .* (YMID .> -1)
        Mel = Mel / (Δx * Δy * sum(Mel))
        return Mel, Xan
    elseif type == "ReplicationExptUniformPerturbed"
        k, l = size(XMID)

        Mel = ones(size(XMID))
        Mel = Mel + 0.05 * (rand(k, l) .- 0.5)
        # Mel = Mel .* rand(k, l)
        Mel = 400 * 9 * Mel / (Δx * Δy * sum(Mel))

        Xan = ones(size(XMID))
        Xan = Xan + 0.05 * (rand(k, l) .- 0.5)
        # Xan = Xan .* rand(k, l)
        Xan = 772 * 9 * Xan / (Δx * Δy * sum(Xan))

        return Mel, Xan

    elseif type == "ReplicationExptPerturbedStripes"
        # The interstripe distance is 0.082 (Alexandria's code, Ref.)
        # Typical width of a Melanophore stripe is 0.5 (Alexandria's code, Ref.)
        # Typical width of a Xanthophore stripe is 0.36 (Alexandria's code, Ref)
        Mel = zeros(size(XMID))
        Xan = zeros(size(YMID))
        for i = 1:3
            lower_Xan = -1.5 + (i-1)*(0.5 + 2*0.082 + 0.36)
            upper_Xan = -1.5 + 0.36 + (i-1)*(2*0.082 + 0.36 + 0.5)
            Xan .+= PlaceHorizontalStripe(lower_Xan, upper_Xan, XMID, YMID, "Xan")

            lower_Mel = -1.5 + 0.36 + 0.082 + (i-1)*(0.5 + 2*0.082 + 0.36)
            upper_Mel = -1.5 + 0.36 + 0.082 + 0.5 + (i-1)*(2*0.082 + 0.36 + 0.5)
            Mel .+= PlaceHorizontalStripe(lower_Mel, upper_Mel, XMID, YMID, "Mel")
        end

        return Mel, Xan

    elseif type == "UniformPerturbedWholeDomainMX"
        k, l = size(XMID)
        Mel = 1 .+ 0 * XMID
        Mel = max.( Mel .+ 10 * (rand(k, l) .- 0.5), 0)
        Mel = 400 * 3*3 * Mel / (Δx * Δy * sum(Mel))


        Xan = 1 .+ 0 * XMID
        Xan = max.( Xan .+ 10 * (rand(k, l) .- 0.5), 0)
        Xan = 772 * 3*3 * Xan / (Δx * Δy * sum(Xan))

        return Mel, Xan
    elseif type == "CheckerboardBio"
        Mel = zeros(size(XMID))
        Xan = zeros(size(YMID))
        dx = (XMID[1,2] - XMID[1,1]);
        dy = (YMID[2,1] - XMID[1,1]);

        dG  = 0.082; # The interstripe distance is 0.082 (Alexandria's code, Ref.) (gap)
        dM  = 0.5;# Typical width of a Melanophore stripe is 0.5 (Alexandria's code, Ref.)
        dX  = 0.36;# Typical width of a Xanthophore stripe is 0.36 (Alexandria's code, Ref)
        ofs = -0.983;#
        MelDensity = 400;
        XanDensity = 772;

        for j = 1:2 # num patch (columns - either mel or xan)
            for i=1:2
                Mel += PlaceBox(ofs + (i-1)*(dX+dG) + (j-1)*(2*dG+dX+dM), ofs + (i-1)*(dX+dG),                  ofs + dM + (i-1)*(dG+dX) + (j-1)*(2*dG+dX+dM), ofs + (i-1)*(dX+dG) + dM, XMID, YMID, "Mel")
                Mel += PlaceBox(ofs + (i-1)*(dX+dG) + (j-1)*(2*dG+dX+dM), ofs + (i-1)*(dX+dG) + 2*dG + dX + dM, ofs + dM + (i-1)*(dG+dX) + (j-1)*(2*dG+dX+dM), ofs + (i-1)*(dX+dG) + 2*(dG+dM) + dX, XMID, YMID, "Mel")

                Xan += PlaceBox(ofs - (i-1)*(dM+dG) + (j-1)*(2*dG+dM+dX) + (dM+dG), ofs + (i-1)*(dM+dG),                ofs  +(dM+dG) - (i-1)*(dM+dG) + dX + (j-1)*(2*dG+dM+dX), ofs + dX + (i-1)*(dM+dG), XMID, YMID, "Xan")
                Xan += PlaceBox(ofs - (i-1)*(dM+dG) + (j-1)*(2*dG+dM+dX) + (dM+dG), ofs + (i-1)*(dM+dG)+1*dM + 2*dG+dX, ofs + (dM+dG) - (i-1)*(dG+dM) + dX + (j-1)*(2*dG+dM+dX), ofs + dX + (i-1)*(dM+dG)+(dM+2*dG+dX), XMID, YMID, "Xan")
            end
        end

        Mel = Mel .> 1;
        Xan = Xan .> 1;

        MelArea = dx*dy*sum(Mel)
        XanArea = dx*dy*sum(Xan)

        Mel = Mel / (dx*dy*sum(Mel)) * MelArea * MelDensity
        Xan = Xan / (dx*dy*sum(Xan)) * XanArea * XanDensity

        return Mel, Xan
    elseif type == "CheckerboardEqui"
        Mel = zeros(size(XMID))
        Xan = zeros(size(YMID))
        dx = (XMID[1,2] - XMID[1,1]);
        dy = (YMID[2,1] - XMID[1,1]);

        dG  = 0.082; # The interstripe distance is 0.082 (Alexandria's code, Ref.) (gap)
        dM  = 0.5;# Typical width of a Melanophore stripe is 0.5 (Alexandria's code, Ref.)
        dX  = 0.36;# Typical width of a Xanthophore stripe is 0.36 (Alexandria's code, Ref.)
        dh = 3 / 8;
        for i = 0:3
            for j = 0:3
                x_i = -1.5 + (2i+1)*dh;
                y_j = -1.5 + (2j+1)*dh;
                if mod(i+j,2) == 0 # melano
                    #println(x_i,y_i)
                    Mel = Mel .+ PlaceBox(x_i - dM/2, y_j - dM/2, x_i + dM/2, y_j + dM/2, XMID, YMID, "Mel")
                else
                    Xan = Xan .+ PlaceBox(x_i - dX/2, y_j - dX/2, x_i + dX/2, y_j + dX/2, XMID, YMID, "Xan")
                end
            end
        end
        return Mel, Xan
    elseif type == "SmallCheckerboardEqui"
        Mel = zeros(size(XMID))
        Xan = zeros(size(YMID))
        dx = (XMID[1,2] - XMID[1,1]);
        dy = (YMID[2,1] - XMID[1,1]);

        dG  = 0.082; # The interstripe distance is 0.082 (Alexandria's code, Ref.) (gap)
        dM  = 0.5;# Typical width of a Melanophore stripe is 0.5 (Alexandria's code, Ref.)
        dX  = 0.36;# Typical width of a Xanthophore stripe is 0.36 (Alexandria's code, Ref.)

        a = 0.5
        dh = 3 / 8 / 2;
        for i = 0:7
            for j = 0:7
                x_i = -1.5 + (2i+1) * dh;
                y_j = -1.5 + (2j+1) * dh;
                if mod(i+j,2) == 0 # melano
                    #println(x_i,y_i)
                    Mel = Mel .+ PlaceBox(x_i - dM/2*a, y_j - dM/2*a, x_i + dM/2*a, y_j + dM/2*a, XMID, YMID, "Mel")
                else
                    Xan = Xan .+ PlaceBox(x_i - dX/2*a, y_j - dX/2*a, x_i + dX/2*a, y_j + dX/2*a, XMID, YMID, "Xan")
                end
            end
        end
        return Mel, Xan
    elseif type == "ABMData"
        sio = pyimport("scipy.io");
        # ABM_Data = sio.loadmat("../../ZF_Migration_ABM/data/sims/MC/MonteCarlo_1Stripe_500_runs_T=150.mat");
        # ABM_Data = sio.loadmat("../../ZF_Migration_ABM/data/sims/MC/MonteCarlo_1Stripe_5000_runs_T=150.mat");
        ABM_Data = sio.loadmat("../../ZF_Migration_ABM/data/sims/MC/ABM_Histogram_Square1_1000_runs_NumBins=120_T=150.mat");
        Mel = ABM_Data["Mel"]; # only care about the matrix (Array{Float64,3} type)
        Mel = Mel[:,:,1]';
        Xan = zeros(size(Mel)); # can include the xanthophore data from the ABM later
        return Mel, Xan
    elseif type == "XanBox"
        Xan = (XMID .< 0.5) .* (XMID .> -0.5) .* (YMID .> -0.5) .* (YMID .< 0.5)
        Xan = 772 * Xan / (Δx * Δy * sum(Xan)) # with xanthophore density of 772 cells/mm^2 (number density)
        # Xan = Xan / (Δx * Δy * sum(Xan)) # normalise to 1 (probability density)
        Mel = zeros(size(YMID))

        return Mel, Xan
    elseif type == "XanStripe"
        Xan = (YMID .< 0.25) .* (YMID .> -0.25) # This is where the error is: for certain mesh sizes, fractions of voxels can be covered by cells and this function does not capture that
        Xan = 772 * 0.5 * 3* Xan / (Δx * Δy * sum(Xan)) # with Xanthophore density of 772 cells/mm^2 (number density)
        # Xan = Xan / (Δx * Δy * sum(Xan)) # normalise to 1 (probability density)
        Mel = zeros(size(YMID))

        return Mel, Xan
    elseif type == "OffsetRectangles"
        Mel = (XMID .< -0.1) .* (XMID .> -0.6) .* (YMID .> -0.75) .* (YMID .< 0.25)
        Mel = 200 * Mel / (Δx * Δy * sum(Mel))

        Xan = (XMID .< 0.6) .* (XMID .> 0.1) .* (YMID .> -0.25) .* (YMID .< 0.75)
        Xan = 772/2 * Xan / (Δx * Δy * sum(Xan))
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
