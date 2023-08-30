function rhs(t, y, p)
    Mel = y[1]
    Xan = y[2]

    WMM = p[1]
    WMX = p[2]
    WXM = p[3]
    WXX = p[4]
    params = p[5]


    @unpack α, Δt, T, IC, L₁, L₂, N₁, N₂, ν, ϵ = params
    Δx = L₁ / N₁
    Δy = L₂ / N₂

    sMel, ~ = size(Mel)
    k₁, ~ = size(WMM)
    k = Int((k₁ - 1) / 2)

    Mel_extended = extend(Mel, k)
    Xan_extended = extend(Xan, k)

    # second argument, k, is the number of rows
    # we want to extend the density matrix by
    # this is determined by the cut_off_radius
    # in the main file


    # computing convolutions for migration
    Wₘₘ = conv(WMM, Mel_extended)
    Wₘₓ = conv(WMX, Xan_extended)
    Wₓₘ = conv(WXM, Mel_extended)
    Wₓₓ = conv(WXX, Xan_extended)

    Wₘₘ = Δx * Δy * Wₘₘ[2*k+1:end-2*k, 2*k+1:end-2*k]
    Wₘₓ = Δx * Δy * Wₘₓ[2*k+1:end-2*k, 2*k+1:end-2*k]
    Wₓₓ = Δx * Δy * Wₓₓ[2*k+1:end-2*k, 2*k+1:end-2*k]
    Wₓₘ = Δx * Δy * Wₓₘ[2*k+1:end-2*k, 2*k+1:end-2*k]

    ξ_Mel = ϵ * (Mel + Xan) + Wₘₘ + Wₘₓ
    ξ_Xan = ϵ * (Mel + Xan) + Wₓₓ + Wₓₘ

    ∂₁ξ_Mel = zeros(2 * N₁ + 1, 2 * N₂)
    ∂₁ξ_Xan = zeros(2 * N₁ + 1, 2 * N₂)

    ∂₁ξ_Mel[2:2*N₁, :] = -(ξ_Mel[2:end, :] - ξ_Mel[1:end-1, :]) / Δx
    ∂₁ξ_Xan[2:2*N₁, :] = -(ξ_Xan[2:end, :] - ξ_Xan[1:end-1, :]) / Δx

    F₁_Mel = zeros(2 * N₁ + 1, 2 * N₂)
    F₁_Xan = zeros(2 * N₁ + 1, 2 * N₂)


    F₁_Mel[2:2*N₁, :] = max.(∂₁ξ_Mel[2:2*N₁, :], 0) .* Mel[1:end-1, :] .+ min.(∂₁ξ_Mel[2:2*N₁, :], 0) .* Mel[2:end, :]
    F₁_Mel[1, :] = max.(∂₁ξ_Mel[2, :], 0) .* Mel[end, :] .+ min.(∂₁ξ_Mel[2*N₂, :], 0) .* Mel[1, :]
    F₁_Mel[2*N₂+1, :] = F₁_Mel[1, :]

    F₁_Xan[2:2*N₁, :] = max.(∂₁ξ_Xan[2:2*N₁, :], 0) .* Xan[1:end-1, :] .+ min.(∂₁ξ_Xan[2:2*N₁, :], 0) .* Xan[2:end, :]
    F₁_Xan[1, :] = max.(∂₁ξ_Xan[2, :], 0) .* Xan[end, :] .+ min.(∂₁ξ_Xan[2*N₂, :], 0) .* Xan[1, :]
    F₁_Xan[2*N₂+1, :] = F₁_Xan[1, :]


    ∂₂ξ_Mel = zeros(2 * N₁, 2 * N₂ + 1)
    ∂₂ξ_Xan = zeros(2 * N₁, 2 * N₂ + 1)

    ∂₂ξ_Mel[:, 2:2*N₂] = -(ξ_Mel[:, 2:end] - ξ_Mel[:, 1:end-1]) / Δy

    ∂₂ξ_Xan[:, 2:2*N₂] = -(ξ_Xan[:, 2:end] - ξ_Xan[:, 1:end-1]) / Δy

    F₂_Mel = zeros(2 * N₁, 2 * N₂ + 1)
    F₂_Xan = zeros(2 * N₁, 2 * N₂ + 1)


    F₂_Mel[:, 2:2*N₂] = max.(∂₂ξ_Mel[:, 2:2*N₂], 0) .* Mel[:, 1:end-1] + min.(∂₂ξ_Mel[:, 2:2*N₂], 0) .* Mel[:, 2:end]
    F₂_Mel[:, 1] = max.(∂₂ξ_Mel[:, 2], 0) .* Mel[:, end] + min.(∂₂ξ_Mel[:, 2*N₂], 0) .* Mel[:, 1]
    F₂_Mel[:, 2*N₂+1] = F₂_Mel[:, 1]

    F₂_Xan[:, 2:2*N₂] = max.(∂₂ξ_Xan[:, 2:2*N₂], 0) .* Xan[:, 1:end-1] + min.(∂₂ξ_Xan[:, 2:2*N₂], 0) .* Xan[:, 2:end]
    F₂_Xan[:, 1] = max.(∂₂ξ_Xan[:, 2], 0) .* Xan[:, end] + min.(∂₂ξ_Xan[:, 2*N₂], 0) .* Xan[:, 1]
    F₂_Xan[:, 2*N₂+1] = F₂_Xan[:, 1]

    ∂₁F₁_Mel = (F₁_Mel[2:end, :] - F₁_Mel[1:end-1, :]) / Δx
    ∂₂F₂_Mel = (F₂_Mel[:, 2:end] - F₂_Mel[:, 1:end-1]) / Δy

    ∂₁F₁_Xan = (F₁_Xan[2:end, :] - F₁_Xan[1:end-1, :]) / Δx
    ∂₂F₂_Xan = (F₂_Xan[:, 2:end] - F₂_Xan[:, 1:end-1]) / Δy

    dMel = -(∂₁F₁_Mel + ∂₂F₂_Mel)
    dXan = -(∂₁F₁_Xan + ∂₂F₂_Xan)

    return [dMel, dXan]
end
