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

# Step 1: Load ABM data 
# ABM_Data     = sio.loadmat("../../ABM/ABM_Histogram_Nbin=30_Ndart=150_BirthOption=1_dt=0.1_T=150.mat");
ABM_Data     = sio.loadmat("../../ABM/ABM_Histogram_Nbin=30_Ndart=150_BirthOption=12_dt=0.1_T=150.mat");
# ABM_Data_Mel = Convert30x30DataInto240x240Data( ABM_Data["MC_Mel"] );
# ABM_Data_Xan = Convert30x30DataInto240x240Data( ABM_Data["MC_Xan"] ); # can include the xanthophore data from the ABM later
ABM_Data_Mel = ABM_Data["MC_Mel"] ;
ABM_Data_Xan = ABM_Data["MC_Xan"] ; # can include the xanthophore data from the ABM later

# Step 2: Initial guess for parameters:
# p = [αₘₘ, γₘ, c₊ₘ]

# p₀ = [0.642, 0.1274, 8.564]; # αₘₘ, γₘ, c₊ from fitting MelBox parameters separately
p₀ = [0.642, 0.1274, 7.0686]; # αₘₘ and γₘ from fitting MelBox parameters separately, c₊ from theoretical estimate
# p₀ = [0.7, 0.2, 14]; # αₘₘ, γₘ, c₊ from fitting MelBox parameters separately


# Step 3: Optimize for alpha using a nonlinear least squares method:
optimal_params = so.least_squares(res_fun_Combined_Mel, p₀, bounds = (0.0,np.inf), args = (ABM_Data_Mel, ABM_Data_Xan));

# Print to terminal the optimized value of α and the termination conditions
# of the algorithm
println("Optimal parameter values are \n");
println(optimal_params["x"]);
println(optimal_params["message"]);

# Step 4: Save the optimal parameters, the residuals for their values, and the jacobian:
sio.savemat("../data/sims/OffsetRectanglesMelOnly/optimal_params_CombinedModel_Ndart=150_BirthOption=12_dt=0.01_runs=1000_NumBins=240_T=150_CloseGuess.mat", optimal_params);
