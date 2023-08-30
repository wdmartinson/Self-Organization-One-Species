using DrWatson
@quickactivate "PDE"

using Plots

include(srcdir("InitialData.jl"))
include(srcdir("Functions.jl"))

function plot_timeseries(params)
      @unpack IC, L₁, L₂, N₁, N₂ = params;
      filename = savename(params, "bson");
      savetag  = savename(params)
      output   = load(datadir("sims", IC, filename));

      XMID, YMID, ~ = set_up_domain(L₁, L₂, N₁, N₂);
      XINT, YINT, ~ = domain_interfaces(L₁, L₂, N₁, N₂);
      x = XMID[1,:];
      y = YMID[:,1];

      t   = output[:t];
      Mel = output[:Mel];
      Xan = output[:Xan];

      for i = 1:1:length(t)
            plot()
            pMel  = contourf(x, y, Mel[i, :, :], c=:grayC, title="melano")
            pXan  = contourf(x, y, Xan[i, :, :], c=cgrad([:white, :orange], [0.1, 0.3, 0.8]), title="xantho")
            plot(pMel, pXan, layout = (1,2), size=(1200,400))
            savefig(datadir("plots", "Timeseries", string(savetag, string("_", string(i)),".png")))
      end
end

function make_melano_movie(params)
      @unpack IC, L₁, L₂, N₁, N₂ = params;
      filename = savename(params, "bson");
      output   = load(datadir("sims", IC, filename));


      XMID, YMID, ~ = set_up_domain(L₁, L₂, N₁, N₂);
      XINT, YINT, ~ = domain_interfaces(L₁, L₂, N₁, N₂);
      x = XMID[1,:];
      y = YMID[:,1];

      t   = output[:t];
      Mel = output[:Mel];
      Xan = output[:Xan];

      anim = @animate for i = 1:1:length(t)
            plot()
            pMel  = contourf(x, y, Mel[i, :, :], c=:grayC, title="t = $(round(t[i], digits = 4))")

      end
      gifname = savename(params, "gif");
      gif(anim, datadir("movies", gifname), fps = 5)
end

function make_melano_movie_or_run(params)
      filename = savename(params, "bson");
      @unpack IC = params;

      if isfile(datadir("sims", IC, filename)) == true # if exists
            make_melano_movie(params);
      else
            run_zebrafish(params);
      end
end

function plot_initial_data(params)
      @unpack IC, L₁, L₂, N₁, N₂ = params

      XMID, YMID, ~ = set_up_domain(L₁, L₂, N₁, N₂)
      XINT, YINT, ~ = domain_interfaces(L₁, L₂, N₁, N₂)
      x = XMID[1,:]
      y = YMID[:,1]

      Mel, Xan = initialise(XMID, YMID; type = IC)

      plot()
      pMel  = contourf(x, y, Mel, c=:grayC, title="melano")
      pXan  = contourf(x, y, Xan, c=cgrad([:white, :orange], [0.1, 0.3, 0.8]), title="xantho")
      plot(pMel, pXan, layout = (1,2), size=(1200,400))

      savefig(datadir("plots", savename(params, "png")))
end

function plot_potential(params, type)
      @unpack L₁, N₁ = params;
      r  = range(0, L₁, length = 10 * N₁)
      W  = ComputePotential(r, type)
      dr = r[2] - r[1]
      rm = (r[2:end] + r[1:end-1]) / 2
      dW = diff(W)

      plot()
      plot!(title = "type: $type", xlabel = "radius")
      plot!(r, W, label = "pot", c = :black)
      plot!(rm, dW / dr, label = "force", c = :red)

      savefig(datadir("plots", "Potentials", string(type, ".png")))
end
