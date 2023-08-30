using PyCall
using Printf
using DSP
using DelimitedFiles
using BSON

include("../src/InitialData.jl");
include("../src/Parameters.jl");
include("../src/Functions.jl");
include("../src/Rhs_zebrafish.jl");

function run_zebrafish_with_output(params) # returns a dictionary with simulation output
      @unpack α, Δt, T, IC, L₁, L₂, N₁, N₂, ν, ϵ = params;

      # set up domain
      XMID, YMID, ~ = set_up_domain(L₁, L₂, N₁, N₂); # cell centres
      ~, ~, RINT = domain_interfaces(L₁, L₂, N₁, N₂) # cell interfaces

      # initialise time horizon
      t̄ = 0;

      # initialise the cell distributions
      MelIni, XanIni = initialise(XMID, YMID, type = IC)

      ȳ = [MelIni, XanIni];

      WMMb = α[1] * ComputePotential(RINT, "MelMel")
      WMXb = α[2] * ComputePotential(RINT, "MelXan")
      WXMb = α[3] * ComputePotential(RINT, "MelXan")
      WXXb = α[4] * ComputePotential(RINT, "XanXan")
      println(maximum(abs.(WXXb)))

      # truncation of the potentials
      CutOffIdx = minimum([20, N₁, N₂]);
      n₁, ~ = size(WMMb) # all potentials will have the same size...
      m₁ = Int((n₁ - 2 * CutOffIdx + 1) / 2)
      m₂ = Int((n₁ + 2 * CutOffIdx + 1) / 2)

      WMM  = WMMb[m₁:m₂, m₁:m₂];
      WMX  = WMXb[m₁:m₂, m₁:m₂];
      WXM  = WXMb[m₁:m₂, m₁:m₂];
      WXX  = WXXb[m₁:m₂, m₁:m₂];

      # parameters passed to the solver
      p = [WMM, WMX, WXM, WXX, params];

      NumTimeSteps = Int(ceil((T / Δt))) + 1; # The number of time steps is consistent with the ABM
      NumOutputs   = Int(ceil((T / Δt) / ν)) + 1; # The number of time steps is consistent with the ABM

      t   = Array{Float32}(undef, NumOutputs);
      Mel = Array{Float32}(undef, NumOutputs, 2*N₁, 2*N₂);
      Xan = Array{Float32}(undef, NumOutputs, 2*N₁, 2*N₂);

      j = 1; # output counter
      for i = 1:NumTimeSteps # ODE solver
            if mod(i-1, ν) == 0 # save output to array at times t = 0, νΔt, 2νΔt, ...
                  t[j] = t̄;
                  Mel[j, :, :] = ȳ[1];
                  Xan[j, :, :] = ȳ[2];
                  j += 1;
            end

            # update pigment cells
            t̄, ȳ = EEuStep(t̄, ȳ, p, Δt, rhs)
      end

      # save everything:
      data     = @dict t Mel Xan;

      return data;
end




function run_zebrafish(params)
      @unpack α, Δt, T, IC, L₁, L₂, N₁, N₂, ν, ϵ = params;

      # set up domain
      XMID, YMID, ~ = set_up_domain(L₁, L₂, N₁, N₂); # cell centres
      ~, ~, RINT = domain_interfaces(L₁, L₂, N₁, N₂) # cell interfaces

      # initialise time horizon
      t̄ = 0;

      # initialise the cell distributions
      MelIni, XanIni = initialise(XMID, YMID, type = IC)
      ȳ = [MelIni, XanIni];

      WMMb = α[1] * ComputePotential(RINT, "MelMel")
      WMXb = α[2] * ComputePotential(RINT, "MelXan")
      WXMb = α[3] * ComputePotential(RINT, "MelXan")
      WXXb = α[4] * ComputePotential(RINT, "XanXan")

      println(maximum(abs.(WXXb)))

      # truncation of the potentials
      CutOffIdx = minimum([20, N₁, N₂]);
      n₁, ~ = size(WMMb) # all potentials will have the same size...
      m₁ = Int((n₁ - 2 * CutOffIdx + 1) / 2)
      m₂ = Int((n₁ + 2 * CutOffIdx + 1) / 2)

      WMM = WMMb[m₁:m₂, m₁:m₂];
      WMX = WMXb[m₁:m₂, m₁:m₂];
      WXM = WXMb[m₁:m₂, m₁:m₂];
      WXX = WXXb[m₁:m₂, m₁:m₂];

      # parameters passed to the solver
      p = [WMM, WMX, WXM, WXX, params];

      NumTimeSteps = Int(ceil((T / Δt))) + 1; # The number of time steps is consistent with the ABM
      NumOutputs   = Int(ceil((T / Δt) / ν)) + 1; # The number of time steps is consistent with the ABM

      t   = Array{Float32}(undef, NumOutputs);
      Mel = Array{Float32}(undef, NumOutputs, 2*N₁, 2*N₂);
      Xan = Array{Float32}(undef, NumOutputs, 2*N₁, 2*N₂);

      j = 1; # output counter
      for i = 1:NumTimeSteps # ODE solver
            if mod(i-1, ν) == 0 # save output to array at times t = 0, νΔt, 2νΔt, ...
                  t[j] = t̄;
                  Mel[j, :, :] = ȳ[1];
                  Xan[j, :, :] = ȳ[2];
                  j += 1;
            end

            # update pigment cells
            t̄, ȳ = EEuStep(t̄, ȳ, p, Δt, rhs)
      end

      # save everything:
      data     = @dict t Mel Xan;
      filename = savename(params, "bson");
      save(datadir("sims", IC, filename), data)
end
