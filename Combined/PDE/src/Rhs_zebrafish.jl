function rhs(t, y, p; parallel=false)
      Mel = y[1]
      Xan = y[2]

      params = p[5]

      @unpack α, Δt, T, IC, L₁, L₂, N₁, N₂, ν, ϵ, c₊, c₋, γ, N_bir = params
      Δx = L₁ / N₁
      Δy = L₂ / N₂

      sMel, ~ = size(Mel)
      k₁, ~ = size(p[1])
      k = Int((k₁ - 1) / 2)

      # second argument, k, is the number of rows
      # we want to extend the density matrix by
      # this is determined by the cut_off_radius
      # in the main file
      Mel_extended = extend(Mel, k)
      Xan_extended = extend(Xan, k)

      # computing local cell counts for nonlocal growth dynamics
      IndB = p[6]

      ## STEP (1): compute cell number at each voxel with IndB:
      LMC = @spawn begin
            buffer = conv(IndB, Mel_extended)
            Δx * Δy * buffer[2*k+1:end-2*k, 2*k+1:end-2*k]
      end

      LXC = @spawn begin
            buffer = conv(IndB, Xan_extended)
            Δx * Δy * buffer[2*k+1:end-2*k, 2*k+1:end-2*k]
      end






      resMel = @spawn begin
            ###### MELANOPHORES ONLY #######

            WMM = p[1]

            # computing convolutions for melano migration
            Wₘₘ = conv(WMM, Mel_extended)

            # truncate convolutions back to domain
            Wₘₘ = Δx * Δy * Wₘₘ[2*k+1:end-2*k, 2*k+1:end-2*k]

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



            # Growth Dynamics for biological case
            LocalMelCount = fetch(LMC)
            # LocalXanCount = fetch(LXC)
            # AnnulMelCount = fetch(AMC)
            # AnnulXanCount = fetch(AXC)

            # Growth dynamics for 1 species
            dMel += γ * N_bir * (c₋ .< LocalMelCount) .* (LocalMelCount .< c₊)

      end

      resXan = @spawn begin
            ###### XANTHOPHORES ONLY #######

            WXX = p[4]

            # computing convolutions for xantho migration
            Wₓₓ = conv(WXX, Xan_extended)

            # truncate convolutions back to domain
            Wₓₓ = Δx * Δy * Wₓₓ[2*k+1:end-2*k, 2*k+1:end-2*k]






            ξ_Xan = ϵ * ( Xan ) + Wₓₓ
            ∂₁ξ_Xan = zeros(2 * N₁ + 1, 2 * N₂)
            ∂₁ξ_Xan[2:2*N₁, :] = -(ξ_Xan[2:end, :] - ξ_Xan[1:end-1, :]) / Δx

            F₁_Xan = zeros(2 * N₁ + 1, 2 * N₂)
            F₁_Xan[2:2*N₁, :] =
                  max.(∂₁ξ_Xan[2:2*N₁, :], 0) .* Xan[1:end-1, :] .+
                  min.(∂₁ξ_Xan[2:2*N₁, :], 0) .* Xan[2:end, :]
            F₁_Xan[1, :] =
                  max.(∂₁ξ_Xan[2, :], 0) .* Xan[end, :] .+ min.(∂₁ξ_Xan[2*N₂, :], 0) .* Xan[1, :]
            F₁_Xan[2*N₂+1, :] = F₁_Xan[1, :]

            ∂₂ξ_Xan = zeros(2 * N₁, 2 * N₂ + 1)
            ∂₂ξ_Xan[:, 2:2*N₂] = -(ξ_Xan[:, 2:end] - ξ_Xan[:, 1:end-1]) / Δy

            F₂_Xan = zeros(2 * N₁, 2 * N₂ + 1)

            F₂_Xan[:, 2:2*N₂] =
                  max.(∂₂ξ_Xan[:, 2:2*N₂], 0) .* Xan[:, 1:end-1] +
                  min.(∂₂ξ_Xan[:, 2:2*N₂], 0) .* Xan[:, 2:end]
            F₂_Xan[:, 1] =
                  max.(∂₂ξ_Xan[:, 2], 0) .* Xan[:, end] + min.(∂₂ξ_Xan[:, 2*N₂], 0) .* Xan[:, 1]
            F₂_Xan[:, 2*N₂+1] = F₂_Xan[:, 1]


            ∂₁F₁_Xan = (F₁_Xan[2:end, :] - F₁_Xan[1:end-1, :]) / Δx
            ∂₂F₂_Xan = (F₂_Xan[:, 2:end] - F₂_Xan[:, 1:end-1]) / Δy


            dXan = -(∂₁F₁_Xan + ∂₂F₂_Xan)


            # Growth dynamics for 1 species

            # LocalMelCount = fetch(LMC)
            LocalXanCount = fetch(LXC)
            # AnnulMelCount = fetch(AMC)
            # AnnulXanCount = fetch(AXC)

            dXan += γ * N_bir * (c₋ .< LocalXanCount) .* (LocalXanCount .< c₊)
      end


      return [fetch(resMel), fetch(resXan)]

end
