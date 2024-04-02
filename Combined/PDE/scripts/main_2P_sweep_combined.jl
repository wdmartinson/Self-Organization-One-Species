using DrWatson
@quickactivate "Macro"

include(srcdir("Parameters.jl"));
include(srcdir("Run_zebrafish.jl"));
include(srcdir("Visualisations.jl"));
include(srcdir("ResFuns.jl"))


using PyCall
np  = pyimport("numpy");
sio = pyimport("scipy.io"); # import scipy.io as sio
so  = pyimport("scipy.optimize"); # from scipy.optimize import least_squares; from scipy.optimize import MATLAB STUFF (abm data)

using Threads

# Step 1: Load ABM data 
ABM_Data     = sio.loadmat("../../ABM/ABM_Histogram_Nbin=30_Ndart=150_BirthOption=1_dt=0.1_T=150.mat");
# ABM_Data     = sio.loadmat("../../ABM/ABM_Histogram_Nbin=30_Ndart=150_BirthOption=12_dt=0.1_T=150.mat");
# ABM_Data_Mel = Convert30x30DataInto240x240Data( ABM_Data["MC_Mel"] );
# ABM_Data_Xan = Convert30x30DataInto240x240Data( ABM_Data["MC_Xan"] ); # can include the xanthophore data from the ABM later
ABM_Data_Mel = ABM_Data["MC_Mel"] ;
ABM_Data_Xan = ABM_Data["MC_Xan"] ; # can include the xanthophore data from the ABM later

# Step 2: Loop over different values of αₘₘ, γ:
alpha_MM = Array(range(0.0, 2.0, length = 51));
gamma    = Array(range(0.0, 1.0, length = 26));
L2_squared = zeros(length(alpha_MM), length(gamma));
for i in eachindex(alpha_MM)
  Threads.@threads for j in eachindex(gamma)
      c₊ = 6.628139504174085;
      p  = [alpha_MM[i], gamma[j], c₊];
      residuals = res_fun_Combined_Mel(p, ABM_Data_Mel, ABM_Data_Xan)
      L2_squared[i,j] = sum(residuals.^2);
  end
end

# Step 3: Save the arrays for further post-processing
output = @dict alpha_MM gamma L2_squared
sio.savemat("../data/sims/MelBox/L2_Squared_Error_Nbin=30_Ndart=150_BirthOption=1_dt=0.1_T=150_IC=MelBox.mat", output);
