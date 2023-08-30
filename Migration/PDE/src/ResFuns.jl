# Residual Functions for one parameter (Mel/Xan),

include(srcdir("Parameters.jl"));
include(srcdir("Run_zebrafish.jl"));
include(srcdir("Visualisations.jl"));


# Define a function that calculates the residual between the ABM and PDE results
function res_fun_1P_Mel(α, ABM_Data_Mel, ABM_Data_Xan)
      # The input α is vector-valued, so set
      αₘₘ = α[1];
      αₘₓ = 0;
      αₓₘ = 0;
      αₓₓ = 0;
      αs = [αₘₘ, αₘₓ, αₓₘ, αₓₓ];
      println("αₘₘ = $(αₘₘ), αₘₓ = $(αₘₓ), αₓₘ = $(αₓₘ), αₓₓ = $(αₓₓ)")

      # Set the value of α for the PDE params variable:
      params = LoadParametersWithAlpha(αs);
      @unpack N₁, N₂, L₁, L₂, Δt = params;
      Δx = L₁/N₁;
      Δy = L₂/N₂;


      # Run the PDE system with this value of alpha. Output the data dictionary. Load the files.
      data = run_zebrafish_with_output(params); # data dictionary t, Mel, Xan
      @unpack Mel, Xan = data;

      # Metric 1 (L2 space-time norm):
      # Compute the residuals for the Melanophore and Xanthophore populations. Concatenate them so that you have the l2 norm over the whole probability density.
      # ABM is size N₁ x N₂ x Nₜ, PDE is Nₜ x N₁ x N₂, and we use the loop to resolve these different data structures
      Mel_res = zeros(prod(size(Mel)));
      Xan_res = zeros(prod(size(Xan)));
      for i in 1:size(Mel, 1)
          Mel_res[4*N₁*N₂*(i-1)+1:4*N₁*N₂*i] = sqrt(Δx*Δy).*reduce(vcat, ABM_Data_Mel[:, :, i]' .- Mel[i, :, :]);
          Xan_res[4*N₁*N₂*(i-1)+1:4*N₁*N₂*i] = sqrt(Δx*Δy).*reduce(vcat, ABM_Data_Xan[:, :, i]' .- Xan[i, :, :]);
      end

      ## Uncomment the following lines if you would rather use an L∞ norm
      # Metric 2 (L∞ norm for support, smoothed using tanh functions):
      # Mel_res = zeros(prod(size(Mel)));
      # Xan_res = zeros(prod(size(Xan)));
      # ϵ = 0.01;
      # for i in 1:size(Mel, 1)
      #     Mel_res[4*N₁*N₂*(i-1)+1:4*N₁*N₂*i] = reduce(vcat, tanh.(ABM_Data_Mel[:, :, i]' ./ 400 ./ ϵ) .- tanh.(Mel[i, :, :] ./400 ./ ϵ  ) );
      #     Xan_res[4*N₁*N₂*(i-1)+1:4*N₁*N₂*i] = reduce(vcat, tanh.(ABM_Data_Xan[:, :, i]' ./ 772 ./ ϵ) .- tanh.(Xan[i, :, :] ./772 ./ ϵ  ) );
      # end

      return [Mel_res; Xan_res];

end

  # Define a function that calculates the residual between the ABM and PDE results
  # function res_fun(alpha, ABM_DATA, NumBins)
  function res_fun_1P_Xan(α, ABM_Data_Mel, ABM_Data_Xan)
        # The input α is vector-valued, so set
        αₘₘ = 0;
        αₘₓ = 0;
        αₓₘ = 0;
        αₓₓ = α[1];
        αs = [αₘₘ, αₘₓ, αₓₘ, αₓₓ];
        println("αₘₘ = $(αₘₘ), αₘₓ = $(αₘₓ), αₓₘ = $(αₓₘ), αₓₓ = $(αₓₓ)")
        # Set the value of α for the PDE params variable:
        params = LoadParametersWithAlpha(αs);
        @unpack N₁, N₂, L₁, L₂, Δt, α = params;
        Δx = L₁/N₁;
        Δy = L₂/N₂;
        # Run the PDE system with this value of alpha. Output the data dictionary. Load the files.
        data = run_zebrafish_with_output(params); # data dictionary t, Mel, Xan
        @unpack Mel, Xan = data;
        # Metric 1:
        # Compute the residuals for the Melanophore and Xanthophore populations. Concatenate them so that you have the l2 norm over the whole probability density.
        # ABM is size N₁ x N₂ x Nₜ, PDE is Nₜ x N₁ x N₂, which is why we had to use the following for loop (not sure if there's a faster way)
        Mel_res = zeros(prod(size(Mel)));
        Xan_res = zeros(prod(size(Xan)));
        for i in 1:size(Mel, 1)
            Mel_res[4*N₁*N₂*(i-1)+1:4*N₁*N₂*i] = sqrt(Δx*Δy).*reduce(vcat, ABM_Data_Mel[:, :, i]' .- Mel[i, :, :]);
            Xan_res[4*N₁*N₂*(i-1)+1:4*N₁*N₂*i] = sqrt(Δx*Δy).*reduce(vcat, ABM_Data_Xan[:, :, i]' .- Xan[i, :, :]);
        end


        ## Uncomment the following lines if you would rather use an L∞ norm
        # Metric 2 (L∞ norm for support, smoothed using tanh functions):
        # Mel_res = zeros(prod(size(Mel)));
        # Xan_res = zeros(prod(size(Xan)));
        # ϵ = 0.01;
        # for i in 1:size(Mel, 1)
        #     Mel_res[4*N₁*N₂*(i-1)+1:4*N₁*N₂*i] = reduce(vcat, tanh.(ABM_Data_Mel[:, :, i]' ./ 400 ./ ϵ) .- tanh.(Mel[i, :, :] ./400 ./ ϵ  ) );
        #     Xan_res[4*N₁*N₂*(i-1)+1:4*N₁*N₂*i] = reduce(vcat, tanh.(ABM_Data_Xan[:, :, i]' ./ 772 ./ ϵ) .- tanh.(Xan[i, :, :] ./772 ./ ϵ  ) );
        # end


        return [Mel_res; Xan_res];
    end
