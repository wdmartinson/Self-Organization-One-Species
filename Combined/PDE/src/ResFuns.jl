# Residual Functions for one parameter (Mel/Xan),
# two parameters (Mel/Xan cross-terms), and
# four parameters for the full system

include(srcdir("Parameters.jl"));
include(srcdir("Run_zebrafish.jl"));
include(srcdir("Visualisations.jl"));


# Define a function that calculates the residual between the ABM and PDE results
# function res_fun(alpha, ABM_DATA, NumBins)
function res_fun_1P_Mel(α, ABM_Data_Mel, ABM_Data_Xan)
      # The input α is vector-valued, so set
      αₘₘ = α[1];
      αₘₓ = 0;
      αₓₘ = 0;
      αₓₓ = 0;
      # αs = [α[1], 0, 0, 0]; # MM, MX, XM, XX
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

      # Metric 1:
      # Compute the residuals for the Melanophore and Xanthophore populations. Concatenate them so that you have the l2 norm over the whole probability density.
      # ABM is size N₁ x N₂ x Nₜ, PDE is Nₜ x N₁ x N₂, which is why we had to use the following for loop (not sure if there's a faster way)
      Mel_res = zeros(prod(size(Mel)));
      Xan_res = zeros(prod(size(Xan)));
      for i in 1:size(Mel, 1)
          Mel_res[4*N₁*N₂*(i-1)+1:4*N₁*N₂*i] = sqrt(Δx*Δy).*reduce(vcat, ABM_Data_Mel[:, :, i]' .- Mel[i, :, :]);
          Xan_res[4*N₁*N₂*(i-1)+1:4*N₁*N₂*i] = sqrt(Δx*Δy).*reduce(vcat, ABM_Data_Xan[:, :, i]' .- Xan[i, :, :]);
      end


      # # Metric 2:
      # # Compute L2-norm for each time point, and sum over the time points:
      # Mel_res = zeros(size(Mel, 1));
      # Xan_res = zeros(size(Xan, 1));
      # for i in 1:size(Mel,1)
      #     Mel_res[i] = sqrt(Δx*Δy*sum((ABM_Data_Mel[:, :, i]' .- Mel[i, :, :]).^2));
      #     Xan_res[i] = sqrt(Δx*Δy*sum((ABM_Data_Xan[:, :, i]' .- Xan[i, :, :]).^2));
      # end


      # Metric 3:
      # Compute the residuals for the Melanophore and Xanthophore populations. Concatenate them so that you have the l2 norm over the whole probability density.
      # ABM is size N₁ x N₂ x Nₜ, PDE is Nₜ x N₁ x N₂, which is why we had to use the following for loop (not sure if there's a faster way)
      # Mel_res = zeros(prod(size(Mel)));
      # Xan_res = zeros(prod(size(Xan)));
      # for i in 1:size(Mel, 1)
      #     Mel_res[4*N₁*N₂*(i-1)+1:4*N₁*N₂*i] = sqrt(Δx*Δy).*reduce(vcat, (ABM_Data_Mel[:, :, i]' .> (3/240).^2) .- (Mel[i, :, :] .> (3/240).^2) );
      #     Xan_res[4*N₁*N₂*(i-1)+1:4*N₁*N₂*i] = sqrt(Δx*Δy).*reduce(vcat, (ABM_Data_Xan[:, :, i]' .> (3/240).^2) .- (Xan[i, :, :] .> (3/240).^2) );
      # end

      # Metric 4: Area-covered metric (smoothed using tanh functions)
      # Compute the residuals for the Melanophore and Xanthophore populations. Concatenate them so that you have the l2 norm over the whole probability density.
      # ABM is size N₁ x N₂ x Nₜ, PDE is Nₜ x N₁ x N₂, which is why we had to use the following for loop (not sure if there's a faster way)
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
      # αs = [α[1], 0, 0, 0]; # MM, MX, XM, XX
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


      # # Metric 2:
      # # Compute L2-norm for each time point, and sum over the time points:
      # Mel_res = zeros(size(Mel, 1));
      # Xan_res = zeros(size(Xan, 1));
      # for i in 1:size(Mel,1)
      #     Mel_res[i] = sqrt(Δx*Δy*sum((ABM_Data_Mel[:, :, i]' .- Mel[i, :, :]).^2));
      #     Xan_res[i] = sqrt(Δx*Δy*sum((ABM_Data_Xan[:, :, i]' .- Xan[i, :, :]).^2));
      # end

      # Metric 4: Area-covered metric (smoothed using tanh functions)
      # Compute the residuals for the Melanophore and Xanthophore populations. Concatenate them so that you have the l2 norm over the whole probability density.
      # ABM is size N₁ x N₂ x Nₜ, PDE is Nₜ x N₁ x N₂, which is why we had to use the following for loop (not sure if there's a faster way)
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
function res_fun_2P(α, ABM_Data_Mel, ABM_Data_Xan)
      println(α);
      # αₘₘ = 0.4062; # output of LSQ for single species (MEL), 60 x 60 mesh
      # αₘₘ = 0.5296; # output of LSQ for single species (MEL), 100 x 100 mesh
      # αₘₘ = 0.5801; # output of LSQ for single species (MEL), 120 x 120 mesh
      αₘₘ = 0.6276; # output of LSQ for single species (MEL), 240 x 240 mesh
      αₘₓ = α[1];
      αₓₘ = α[2];
      # αₓₓ = 0.1400; # output of LSQ for single species (XAN), 60 x 60 mesh
      # αₓₓ = 0.2928; # output of LSQ for single species (XAN), 100 x 100 mesh
      # αₓₓ = 0.3568; # output of LSQ for single species (XAN), 120 x 120 mesh
      αₓₓ = 0.4809; # output of LSQ for single species (XAN), 240 x 240 mesh

      # αₘₘ = 0.0183; # output of LSQ for single species (MEL), 60 x 60 mesh
      # αₘₘ = 0.0236; # output of LSQ for single species (MEL), 120 x 120 mesh
      # αₘₘ = 0.0245; # output of LSQ for single species (MEL), 240 x 240 mesh
      # αₘₓ = α[1];
      # αₓₘ = α[2];
      # αₓₓ = 0.1361; # output of LSQ for single species (XAN), 60 x 60 mesh
      # αₓₓ = 0.3502; # output of LSQ for single species (XAN), 120 x 120 mesh
      # αₓₓ = 0.4766; # output of LSQ for single species (XAN), 240 x 240 mesh

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

      # Metric 1:
      # Compute the residuals for the Melanophore and Xanthophore populations. Concatenate them so that you have the l2 norm over the whole probability density.
      # ABM is size N₁ x N₂ x Nₜ, PDE is Nₜ x N₁ x N₂, which is why we had to use the following for loop (not sure if there's a faster way)
      Mel_res = zeros(prod(size(Mel)));
      Xan_res = zeros(prod(size(Xan)));
      for i in 1:size(Mel, 1)
          Mel_res[4*N₁*N₂*(i-1)+1:4*N₁*N₂*i] = sqrt(Δx*Δy).*reduce(vcat, ABM_Data_Mel[:, :, i]' .- Mel[i, :, :]);
          Xan_res[4*N₁*N₂*(i-1)+1:4*N₁*N₂*i] = sqrt(Δx*Δy).*reduce(vcat, ABM_Data_Xan[:, :, i]' .- Xan[i, :, :]);
      end


      # # Metric 2:
      # # Compute L2-norm for each time point, and sum over the time points:
      # Mel_res = zeros(size(Mel, 1));
      # Xan_res = zeros(size(Xan, 1));
      # for i in 1:size(Mel,1)
      #     Mel_res[i] = sqrt(Δx*Δy*sum((ABM_Data_Mel[:, :, i]' .- Mel[i, :, :]).^2));
      #     Xan_res[i] = sqrt(Δx*Δy*sum((ABM_Data_Xan[:, :, i]' .- Xan[i, :, :]).^2));
      # end

      return [Mel_res; Xan_res];
end


function res_fun_Combined_Mel(p, ABM_Data_Mel, ABM_Data_Xan)

      Ndart = 150;

      # The input α is vector-valued, so set
      αₘₘ = p[1];
      αₘₓ = 0;
      αₓₘ = 0;
      αₓₓ = 0;
      γₘ  = p[2];
      c₊  = p[3];
      # αs = [p[1], 0, 0, 0]; # MM, MX, XM, XX
      αs = [αₘₘ, αₘₓ, αₓₘ, αₓₓ];
      println("αₘₘ = $(αₘₘ), γₘ = $(γₘ), c₊ = $(c₊)")
      # Set the value of α for the PDE params variable:
      params = LoadCombinedModelParameters(αs, γₘ, 0.0, c₊, 2.0, Ndart);
      @unpack N₁, N₂, L₁, L₂, Δt = params;


      # Run the PDE system with this value of alpha. Output the data dictionary. Load the files.
      data = run_zebrafish_with_output(params); # data dictionary t, Mel, Xan
      @unpack Mel, Xan, CoarserMel, CoarserXan = data;

      N₁ = 15;
      N₂ = 15; 
      Mel = CoarserMel;
      Xan = CoarserXan;
      
      Δx = L₁/N₁;
      Δy = L₂/N₂;

      # Metric 1:
      # Compute the residuals for the Melanophore and Xanthophore populations. Concatenate them so that you have the l2 norm over the whole probability density.
      # ABM is size N₁ x N₂ x Nₜ, PDE is Nₜ x N₁ x N₂, which is why we had to use the following for loop (not sure if there's a faster way)
      Mel_res = zeros(prod(size(Mel)));
      Xan_res = zeros(prod(size(Xan)));
      for i in 1:size(Mel, 1)
          Mel_res[4*N₁*N₂*(i-1)+1:4*N₁*N₂*i] = sqrt(Δx*Δy).*reduce(vcat, ABM_Data_Mel[:, :, i]' .- Mel[i, :, :]);
          Xan_res[4*N₁*N₂*(i-1)+1:4*N₁*N₂*i] = sqrt(Δx*Δy).*reduce(vcat, ABM_Data_Xan[:, :, i]' .- Xan[i, :, :]);
      end


      # # Metric 2:
      # # Compute L2-norm for each time point, and sum over the time points:
      # Mel_res = zeros(size(Mel, 1));
      # Xan_res = zeros(size(Xan, 1));
      # for i in 1:size(Mel,1)
      #     Mel_res[i] = sqrt(Δx*Δy*sum((ABM_Data_Mel[:, :, i]' .- Mel[i, :, :]).^2));
      #     Xan_res[i] = sqrt(Δx*Δy*sum((ABM_Data_Xan[:, :, i]' .- Xan[i, :, :]).^2));
      # end


      # Metric 3:
      # Compute the residuals for the Melanophore and Xanthophore populations. Concatenate them so that you have the l2 norm over the whole probability density.
      # ABM is size N₁ x N₂ x Nₜ, PDE is Nₜ x N₁ x N₂, which is why we had to use the following for loop (not sure if there's a faster way)
      # Mel_res = zeros(prod(size(Mel)));
      # Xan_res = zeros(prod(size(Xan)));
      # for i in 1:size(Mel, 1)
      #     Mel_res[4*N₁*N₂*(i-1)+1:4*N₁*N₂*i] = sqrt(Δx*Δy).*reduce(vcat, (ABM_Data_Mel[:, :, i]' .> (3/240).^2) .- (Mel[i, :, :] .> (3/240).^2) );
      #     Xan_res[4*N₁*N₂*(i-1)+1:4*N₁*N₂*i] = sqrt(Δx*Δy).*reduce(vcat, (ABM_Data_Xan[:, :, i]' .> (3/240).^2) .- (Xan[i, :, :] .> (3/240).^2) );
      # end

      # Metric 4: Area-covered metric (smoothed using tanh functions)
      # Compute the residuals for the Melanophore and Xanthophore populations. Concatenate them so that you have the l2 norm over the whole probability density.
      # ABM is size N₁ x N₂ x Nₜ, PDE is Nₜ x N₁ x N₂, which is why we had to use the following for loop (not sure if there's a faster way)
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
function res_fun_4P(α, ABM_Data_Mel, ABM_Data_Xan)

      println(α);

      αₘₘ = α[1];
      αₘₓ = α[2];
      αₓₘ = α[3];
      αₓₓ = α[4];

      println("αₘₘ = $(αₘₘ), αₘₓ = $(αₘₓ), αₓₘ = $(αₓₘ), αₓₓ = $(αₓₓ)")

      # Set the value of α for the PDE params variable:
      params = LoadParametersWithAlpha(αₘₘ, αₘₓ, αₓₘ, αₓₓ);
      @unpack N₁, N₂, L₁, L₂, Δt = params;
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


      # # Metric 2:
      # # Compute L2-norm for each time point, and sum over the time points:
      # Mel_res = zeros(size(Mel, 1));
      # Xan_res = zeros(size(Xan, 1));
      # for i in 1:size(Mel,1)
      #     Mel_res[i] = sqrt(Δx*Δy*sum((ABM_Data_Mel[:, :, i]' .- Mel[i, :, :]).^2));
      #     Xan_res[i] = sqrt(Δx*Δy*sum((ABM_Data_Xan[:, :, i]' .- Xan[i, :, :]).^2));
      # end

      return [Mel_res; Xan_res];
end


function res_fun_2P_MultipleDataSets(α, ABM_Data_1_Mel, ABM_Data_1_Xan, ABM_Data_2_Mel, ABM_Data_2_Xan, ABM_Data_3_Mel, ABM_Data_3_Xan)
# function res_fun_2P_MultipleDataSets(α, ABM_Data_1_Mel, ABM_Data_1_Xan, ABM_Data_2_Mel, ABM_Data_2_Xan, ABM_Data_3_Mel, ABM_Data_3_Xan, ABM_Data_4_Mel, ABM_Data_4_Xan)

      println(α);
      αₘₘ = 0.4062; # output of LSQ for single species (MEL), 60 x 60 mesh
      # αₘₘ = 0.5296; # output of LSQ for single species (MEL), 100 x 100 mesh
      # αₘₘ = 0.5801; # output of LSQ for single species (MEL), 120 x 120 mesh
      # αₘₘ = 0.6276; # output of LSQ for single species (MEL), 240 x 240 mesh
      αₘₓ = α[1];
      αₓₘ = α[2];
      αₓₓ = 0.1400; # output of LSQ for single species (XAN), 60 x 60 mesh
      # αₓₓ = 0.2928; # output of LSQ for single species (XAN), 100 x 100 mesh
      # αₓₓ = 0.3568; # output of LSQ for single species (XAN), 120 x 120 mesh
      # αₓₓ = 0.4809; # output of LSQ for single species (XAN), 240 x 240 mesh

      αs = [αₘₘ, αₘₓ, αₓₘ, αₓₓ];
      println("αₘₘ = $(αₘₘ), αₘₓ = $(αₘₓ), αₓₘ = $(αₓₘ), αₓₓ = $(αₓₓ)")

      ICs = ["TwoStripes", "OffsetRectangles", "ReplicationExptPerturbedStripes"];
      # ICs = ["TwoStripes", "OffsetRectangles", "ReplicationExptPerturbedStripes", "Checkerboard"];

      Total_Mel_res = [];
      Total_Xan_res = [];

      for k = 1:length(ICs)
            # Set the value of α for the PDE params variable:
            params = LoadParametersWithAlphaAndIC(αs, ICs[k]);
            @unpack N₁, N₂, L₁, L₂, Δt = params;
            Δx = L₁/N₁;
            Δy = L₂/N₂;

            # Run the PDE system with this value of alpha. Output the data dictionary. Load the files.
            data = run_zebrafish_with_output(params); # data dictionary t, Mel, Xan
            @unpack Mel, Xan = data;

            Mel_res = zeros(prod(size(Mel)));
            Xan_res = zeros(prod(size(Xan)));

            if k == 1
                  ABM_Data_Mel = ABM_Data_1_Mel;
                  ABM_Data_Xan = ABM_Data_1_Xan;
            elseif k == 2
                  ABM_Data_Mel = ABM_Data_2_Mel;
                  ABM_Data_Xan = ABM_Data_2_Xan;
            elseif k == 3
                  ABM_Data_Mel = ABM_Data_3_Mel;
                  ABM_Data_Xan = ABM_Data_3_Xan;
            # else
            #       ABM_Data_Mel = ABM_Data_4_Mel;
            #       ABM_Data_Xan = ABM_Data_4_Xan;
            end

            # Metric 1:
            # Compute the residuals for the Melanophore and Xanthophore populations. Concatenate them so that you have the l2 norm over the whole probability density.
            # ABM is size N₁ x N₂ x Nₜ, PDE is Nₜ x N₁ x N₂, which is why we had to use the following for loop (not sure if there's a faster way)
            for i in 1:size(Mel, 1)
                Mel_res[4*N₁*N₂*(i-1)+1:4*N₁*N₂*i] = sqrt(Δx*Δy).*reduce(vcat, ABM_Data_Mel[:, :, i]' .- Mel[i, :, :]);
                Xan_res[4*N₁*N₂*(i-1)+1:4*N₁*N₂*i] = sqrt(Δx*Δy).*reduce(vcat, ABM_Data_Xan[:, :, i]' .- Xan[i, :, :]);
            end

            # # Metric 2:
            # # Compute L2-norm for each time point, and sum over the time points:
            # Mel_res = zeros(size(Mel, 1));
            # Xan_res = zeros(size(Xan, 1));
            # for i in 1:size(Mel,1)
            #     Mel_res[i] = sqrt(Δx*Δy*sum((ABM_Data_Mel[:, :, i]' .- Mel[i, :, :]).^2));
            #     Xan_res[i] = sqrt(Δx*Δy*sum((ABM_Data_Xan[:, :, i]' .- Xan[i, :, :]).^2));
            # end

            Total_Mel_res = vcat(Total_Mel_res, Mel_res)
            Total_Xan_res = vcat(Total_Xan_res, Xan_res)
      end



      return [Total_Mel_res; Total_Xan_res];
end
