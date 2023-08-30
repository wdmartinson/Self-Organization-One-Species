using DrWatson
@quickactivate "PDE"
include(srcdir("Parameters.jl")); 
include(srcdir("Run_zebrafish.jl"));
include(srcdir("Visualisations.jl"));
include(srcdir("ResFuns.jl"))


using PyCall
sio = pyimport("scipy.io"); # import scipy.io as sio

# Step 1: Load ABM data in julia from ABM MC folder.
ABM_Data     = sio.loadmat("../../ABM/data/sims/MC/ABM_Histograms/ABM_Histogram_MonteCarlo_BirthOption=1_dt=0.01_runs=1000_NumBins=240_T=150.mat");
ABM_Data_Mel = ABM_Data["Mel"];
ABM_Data_Xan = ABM_Data["Xan"];

# Step 2: Loop over different values of αₘ:
alpha_M = Array(range(0.0, 1.0, length = 101));
L2_squared = zeros(length(alpha_M));
for i in 1:length(alpha_MX)
      α = [alpha_M[i]];
      residuals = res_fun_1P(α, ABM_Data_Mel, ABM_Data_Xan);
      L2_squared[i] = sum(residuals.^2);
end

# Step 3: Save the arrays for further post-processing
output = @dict alpha_M L2_squared
sio.savemat("../data/sims/MelBox/L2_Squared_Error_BirthOption=1_dt=0.01_runs=1000_NumBins=240_T=150_IC=MelBox.mat", output);
