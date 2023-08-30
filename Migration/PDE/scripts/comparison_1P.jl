using DrWatson
@quickactivate "PDE"
include(srcdir("Parameters.jl"));
include(srcdir("Run_zebrafish.jl"));
include(srcdir("Visualisations.jl"));
include(srcdir("ResFuns.jl")); # contains the res_fun_1P_Mel, res_fun_1P_Xan functions

using PyCall
np  = pyimport("numpy");
sio = pyimport("scipy.io"); # import scipy.io as sio
so  = pyimport("scipy.optimize"); # from scipy.optimize import least_squares; from scipy.optimize import MATLAB STUFF (abm data)


### Step 1: Load ABM data

ABM_Data     = sio.loadmat("../../ABM/data/sims/MC/ABM_Histograms/ABM_Histogram_MonteCarlo_BirthOption=1_dt=0.01_runs=10000_NumBins=240_T=150.mat");
ABM_Data_Mel = ABM_Data["Chr"];
ABM_Data_Xan = zeros(size(ABM_Data_Mel)); # can include the xanthophore data from the ABM later

# Initial guess for alpha:
# α₀ = 0.5; # for xanthophores
α₀ = 0.6; # for melanophores

# Optimize for alpha using a nonlinear least squares method:
optimal_α = so.least_squares(res_fun_1P_Mel, α₀, bounds = (0.0,np.inf), args = (ABM_Data_Mel, ABM_Data_Xan));

# Print the optimized value of α and the termination conditions of the algorithm
println("Optimal value of α is \n");
println(optimal_α["x"]);
println(optimal_α["message"]);

# Save the value of alpha, the residuals for that value, and the jacobian:
sio.savemat("../data/sims/MelBox/optimal_alpha_BirthOption=1_NumIter=10000_NumBins=240_T=150.mat", optimal_α);
