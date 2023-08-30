using PyCall

function coarsen_histogram(BirthOption, dt, NBin, NumIter, save_yes_no)
    if mod(240, NBin) ~= 0 || NBin >= 240
        error(ArgumentError("Need NBin to be less than 240, and a multiple of 30."))
    end
    sio = pyimport("scipy.io"); # import scipy.io as sio

    # Load the 240 x 240 histogram
    MC_240 = sio.loadmat("../../ABM/data/sims/MC/ABM_Histograms/ABM_Histogram_MonteCarlo_BirthOption=$(BirthOption)_dt=$(dt)_runs=$(NumIter)_NBins=240_T=150.mat", 'MC');
    MC_240 = MC_240["Chr"];

    # Coarsen the mesh based on the value of NBin
    MC = zeros((NBin, NBin, 151));

    g = round(Int, 1.0/(Nbin/240));

    for k = 1:151
        for i = 1:Nbin
            for j = 1:Nbin
                MC[i,j,k] = 1.0/g.^2.*sum( MC_240[(1+g*(i-1)):g*i, (1+g*(j-1)):g*j, k] );
            end
        end
    end

    # If you want to save the coarser mesh, do so.
    if save_yes_no
        sio.savemat("../../ABM/data/sims/MC/ABM_Histograms/ABM_Histogram_MonteCarlo_BirthOption=$(BirthOption)_dt=$(dt)_runs=$(NumIter)_NBins=$(NBin)_T=150.mat", {"MC": MC});
    end

    return MC

end
