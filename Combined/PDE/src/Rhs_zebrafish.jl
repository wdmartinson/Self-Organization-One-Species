function rhs(t, y, p; parallel=false)
      # XanD = kron(ones(3,3), XanD)
      Mel = y[1]
      # Xan = y[2]

      # params = p[5]

      @unpack α, Δt, T, IC, L₁, L₂, N₁, N₂, ν, ϵ, CmaxMel, CminMel, max_cpd, GrowthRateMel, GrowthRateXan = p[5]
      Δx = L₁ / N₁
      Δy = L₂ / N₂

      sMel, ~ = size(Mel)
      k₁, ~ = size(p[1])
      #k₁ = 2 * N₁ + 1;
      k = Int((k₁ - 1) / 2)

      # second argument, k, is the number of rows
      # we want to extend the density matrix by
      # this is determined by the cut_off_radius
      # in the main file
      Mel_extended = extend(Mel, k)
      # Xan_extended = extend(Xan, k)

      # computing local cell counts for nonlocal growth dynamics a la Shekar & Al.
      IndB = p[6]
      # TODO: (1) compute IndB, IndA;
      LMC = Threads.@spawn begin
            buffer = conv(IndB, Mel_extended)
            Δx * Δy * buffer[2*k+1:end-2*k, 2*k+1:end-2*k]
      end

      # LXC = @spawn begin
      #       buffer = conv(IndB, Xan_extended)
      #       Δx * Δy * buffer[2*k+1:end-2*k, 2*k+1:end-2*k]
      # end

      # AMC = @spawn begin
      #       buffer = conv(IndA, Mel_extended)
      #       Δx * Δy * buffer[2*k+1:end-2*k, 2*k+1:end-2*k]
      # end

      # AXC = @spawn begin
      #       buffer = conv(IndA, Xan_extended)
      #       Δx * Δy * buffer[2*k+1:end-2*k, 2*k+1:end-2*k]
      # end





      resMel = Threads.@spawn begin
            #dMel = begin
            ###### MELANOPHORES ONLY #######

            WMM = p[1]
            # WMX = p[2]

            # computing convolutions for melano migration
            Wₘₘ = conv(WMM, Mel_extended)
            # Wₘₓ = @spawn conv(WMX, Xan_extended)


            # truncate convolutions back to domain
            Wₘₘ = Δx * Δy * Wₘₘ[2*k+1:end-2*k, 2*k+1:end-2*k]
            # Wₘₓ = Δx * Δy * fetch(Wₘₓ)[2*k+1:end-2*k, 2*k+1:end-2*k]

            # ξ_Mel = ϵ * (Mel + Xan) + Wₘₘ + Wₘₓ
            ξ_Mel = ϵ * (Mel) + Wₘₘ
            ∂₁ξ_Mel = zeros(2 * N₁ + 1, 2 * N₂)
            ∂₁ξ_Mel[2:2*N₁, :] = -(ξ_Mel[2:end, :] - ξ_Mel[1:end-1, :]) / Δx
            F₁_Mel = zeros(2 * N₁ + 1, 2 * N₂)

            F₁_Mel[2:2*N₁, :] =
                  max.(∂₁ξ_Mel[2:2*N₁, :], 0) .* Mel[1:end-1, :] .+
                  min.(∂₁ξ_Mel[2:2*N₁, :], 0) .* Mel[2:end, :]
            F₁_Mel[1, :] =
                  max.(∂₁ξ_Mel[2, :], 0) .* Mel[end, :] .+ min.(∂₁ξ_Mel[2*N₂, :], 0) .* Mel[1, :]
            F₁_Mel[2*N₂+1, :] = F₁_Mel[1, :]

            ∂₂ξ_Mel = zeros(2 * N₁, 2 * N₂ + 1)
            ∂₂ξ_Mel[:, 2:2*N₂] = -(ξ_Mel[:, 2:end] - ξ_Mel[:, 1:end-1]) / Δy

            F₂_Mel = zeros(2 * N₁, 2 * N₂ + 1)

            F₂_Mel[:, 2:2*N₂] =
                  max.(∂₂ξ_Mel[:, 2:2*N₂], 0) .* Mel[:, 1:end-1] +
                  min.(∂₂ξ_Mel[:, 2:2*N₂], 0) .* Mel[:, 2:end]
            F₂_Mel[:, 1] =
                  max.(∂₂ξ_Mel[:, 2], 0) .* Mel[:, end] + min.(∂₂ξ_Mel[:, 2*N₂], 0) .* Mel[:, 1]
            F₂_Mel[:, 2*N₂+1] = F₂_Mel[:, 1]

            ∂₁F₁_Mel = (F₁_Mel[2:end, :] - F₁_Mel[1:end-1, :]) / Δx
            ∂₂F₂_Mel = (F₂_Mel[:, 2:end] - F₂_Mel[:, 1:end-1]) / Δy

            dMel = -(∂₁F₁_Mel + ∂₂F₂_Mel)


            # # adding growth set_up_domain
            # # TODO: define r, N, cmax, cmin for Xan and Mel -


            # NMel = 5

            # r₂Mel = 0 # short-range death rate constant, Melanophores
            # r₃Mel = 0 # long-range death rate constant, Melanophores
            # ξ = 1.2 #as in Al's note but tbc


            # Growth Dynamics for biological case
            # LocalMelCount = fetch(LMC)
            # LocalXanCount = fetch(LXC)
            # AnnulMelCount = fetch(AMC)
            # AnnulXanCount = fetch(AXC)

            # birth_rate_Mel =
            #       rMel * NMel * (LocalMelCount .< CmaxMel) .* (LocalMelCount .> CminMel) .*
            #       (Mel .< 400.0)
            # death_rate_Mel_short_range = r₂Mel * (LocalXanCount .> LocalMelCount)
            # death_rate_Mel_long_range = r₃Mel * (ξ * AnnulXanCount .< AnnulMelCount)
            # dMel = birth_rate_Mel - death_rate_Mel_short_range - death_rate_Mel_long_range # make this +=

            # # Growth dynamics for test case (1 species paper)
            # #R₀ = 0.075 # From Shekar's code, TODO: Define R₀ in Parameters.jl
            # #CmaxMel = 400.0 * pi * R₀^2 # taken from Shekar
            # # CmaxMel = 8.8494 # Fitted from Reaction-Only PDE data
            # CmaxMel = 7.0686 # 400.0 * pi * R₀^2 -- Maximum value from ABM(? TODO: check with Alexandria)
            # # CmaxMel = 7.0
            # CminMel = 1.0 # taken from Shekar
            # max_cpd = 150 # numbers of cells nucleated per day (taken from Shekar)
            # GrowthRateMel = 0.1221 # Value fitted to 2D Reaction only ABM code # 0.27 # taken from Shekar's code that he sent us, fit to 1D model(?)
            # # GrowthRateMel = 0.27 # taken from Shekar
            # # GrowthRateMel = 0.16 # From 1D code, this doesn't seem to introduce too much error...

            IndB = p[6]
            LocalMelCount = fetch(LMC)


            # dMel += GrowthRateMel * max_cpd * (CminMel .< LocalMelCount) .* (LocalMelCount .< CmaxMel)
            # The following is the same line as above, but you use tanh functions to smooth the indicator functions...
            dMel += GrowthRateMel * max_cpd * (0.5 .+ 0.5 .* tanh.( (LocalMelCount .- CminMel)./0.01 ) ) .* (0.5 .+ 0.5 .* tanh.( (CmaxMel .- LocalMelCount)./0.01 ))

      end

      # resXan = Threads.@spawn begin
      #       #dXan = begin
      #       ###### XANTHOPHORES ONLY #######
      #
      #       WXM = p[3]
      #       WXX = p[4]
      #
      #       # computing convolutions for xantho migration
      #       Wₓₘ = @spawn conv(WXM, Mel_extended)
      #       Wₓₓ = conv(WXX, Xan_extended)
      #
      #       # truncate convolutions back to domain
      #       Wₓₓ = Δx * Δy * Wₓₓ[2*k+1:end-2*k, 2*k+1:end-2*k]
      #       Wₓₘ = Δx * Δy * fetch(Wₓₘ)[2*k+1:end-2*k, 2*k+1:end-2*k]
      #
      #
      #
      #
      #
      #
      #       ξ_Xan = ϵ * (Mel + Xan) + Wₓₓ + Wₓₘ
      #       ∂₁ξ_Xan = zeros(2 * N₁ + 1, 2 * N₂)
      #       ∂₁ξ_Xan[2:2*N₁, :] = -(ξ_Xan[2:end, :] - ξ_Xan[1:end-1, :]) / Δx
      #
      #       F₁_Xan = zeros(2 * N₁ + 1, 2 * N₂)
      #       F₁_Xan[2:2*N₁, :] =
      #             max.(∂₁ξ_Xan[2:2*N₁, :], 0) .* Xan[1:end-1, :] .+
      #             min.(∂₁ξ_Xan[2:2*N₁, :], 0) .* Xan[2:end, :]
      #       F₁_Xan[1, :] =
      #             max.(∂₁ξ_Xan[2, :], 0) .* Xan[end, :] .+ min.(∂₁ξ_Xan[2*N₂, :], 0) .* Xan[1, :]
      #       F₁_Xan[2*N₂+1, :] = F₁_Xan[1, :]
      #
      #       ∂₂ξ_Xan = zeros(2 * N₁, 2 * N₂ + 1)
      #       ∂₂ξ_Xan[:, 2:2*N₂] = -(ξ_Xan[:, 2:end] - ξ_Xan[:, 1:end-1]) / Δy
      #
      #       F₂_Xan = zeros(2 * N₁, 2 * N₂ + 1)
      #
      #       F₂_Xan[:, 2:2*N₂] =
      #             max.(∂₂ξ_Xan[:, 2:2*N₂], 0) .* Xan[:, 1:end-1] +
      #             min.(∂₂ξ_Xan[:, 2:2*N₂], 0) .* Xan[:, 2:end]
      #       F₂_Xan[:, 1] =
      #             max.(∂₂ξ_Xan[:, 2], 0) .* Xan[:, end] + min.(∂₂ξ_Xan[:, 2*N₂], 0) .* Xan[:, 1]
      #       F₂_Xan[:, 2*N₂+1] = F₂_Xan[:, 1]
      #
      #
      #       ∂₁F₁_Xan = (F₁_Xan[2:end, :] - F₁_Xan[1:end-1, :]) / Δx
      #       ∂₂F₂_Xan = (F₂_Xan[:, 2:end] - F₂_Xan[:, 1:end-1]) / Δy
      #
      #
      #       dXan = -(∂₁F₁_Xan + ∂₂F₂_Xan)
      #
      #
      #
      #
      #       # LocalMelCount = fetch(LMC)
      #
      #       # AnnulMelCount = fetch(AMC)
      #       # AnnulXanCount = fetch(AXC)
      #
      #
      #       # NXan = 5
      #       # rXan = 0.27 # set equal to the melanophores
      #       # CmaxXan = 772.0 * pi * R₀^2 # taken from Shekar
      #       # CminXan = 1.0 # taken from Shekar
      #       # r₂Xan = 0 # short-range death rate constant, Xanthophores
      #
      #       # birth_rate_Xan =
      #       #       rXan * NXan * (LocalXanCount .< CmaxXan) .* (LocalXanCount .> CminXan) .*
      #       #       (Xan .< 772.0)
      #       # death_rate_Xan_short_range = r₂Xan * (LocalXanCount .< LocalMelCount)
      #
      #       # dXan = birth_rate_Xan - death_rate_Xan_short_range # make this +=
      #
      #
      #       # # Growth dynamics for test case (1 species paper)
      #       # #R₀ = 0.075 # From Shekar's code, TODO: Define R₀ in Parameters.jl
      #       # #CmaxXan = 400.0 * pi * R₀^2 # taken from Shekar
      #       # CmaxXan = 8.8494
      #       # CminXan = 1.0 # taken from Shekar
      #       # max_cpd = 150 # numbers of cells nucleated per day (taken from Shekar)
      #       # GrowthRateXan = 0.1221 # 0.27 # taken from Shekar
      #       # IndB = p[6]
      #
      #       LocalXanCount = fetch(LXC)
      #
      #
      #       dXan += GrowthRateXan * max_cpd * (CminXan .< LocalXanCount) .* (LocalXanCount .< CmaxXan)
      # end


      # return [fetch(resMel), fetch(resXan)]
      return [fetch(resMel), zeros(2 * N₁, 2 * N₂)]

end

##### REACTION Part
#
#
# XanD = kron(ones(3,3), XanD)
# Mel  = kron(ones(3,3), Mel)
# IriD = kron(ones(3,3), IriD)
#
#
#
#
# MelNonLoc_1 = conv(ANNULUS, XanD + IriD)
# MelNonLoc_1 = Δx * Δy * MelNonLoc_1[3*N₁+1:5*N₁, 3*N₂+1:5*N₂]
#
# MelNonLoc_2 = conv(ANNULUS, Mel)
# MelNonLoc_2 = Δx * Δy * MelNonLoc_2[3*N₁+1:5*N₁, 3*N₂+1:5*N₂]
#
# IriDNonLoc  = conv(ANNULUS, XanD)
# IriDNonLoc  = Δx * Δy * IriDNonLoc[3*N₁+1:5*N₁, 3*N₂+1:5*N₂]
#
#
# Mel   =  Mel[2*N₁+1:4*N₁, 2*N₂+1:4*N₂]
# XanD  = XanD[2*N₁+1:4*N₁, 2*N₂+1:4*N₂]
# IriD  = IriD[2*N₁+1:4*N₁, 2*N₂+1:4*N₂]
#
# Mel   = Mel  + Δt * (M_ϵ*(1 - Mel) +  (1 - Mel) .* MelNonLoc_1  -  Mel .* (XanD + MelNonLoc_2))
# XanD  = XanD + Δt * ((1 - XanD) .* (XanD + IriD) - 1 * XanD .* Mel)
# IriD  = IriD + Δt * ((1 - IriD) .* (XanD + IriD) - IriD .* (Mel + IriDNonLoc))


# XanD = max(XanD, 0)
# Mel = max(Mel, 0)

# t = t + Δt;
