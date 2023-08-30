function new_mesh = coarsen_histogram(BirthOption, dt, NBin, NumIter, save_yes_no)
    if mod(NBin, 30) ~= 0 || NBin >= 240
        error('Need NBin to be less than 240, and a multiple of 30.')
    end
    % Load the 240 x 240 histogram
    load(['../data/sims/MC/ABM_Histograms/ABM_Histogram_MonteCarlo_BirthOption=',num2str(BirthOption),'_runs=',num2str(NumIter),'_dt=',num2str(dt),'_NBins=240_T=150.mat'], 'Chr');

    % Coarsen the mesh based on the value of NBin
    MC = zeros(NBin, NBin, 151);

    g = 1.0./(Nbin/240);

    for k = 1:151
        for i = 1:Nbin
            for j = 1:Nbin
                MC(i,j,k) = 1./g.^2.*sum( Chr((1+g*(i-1)):g*i, (1+g*(j-1)):g*j, k), 'all');
            end
        end
    end
    % If you want to save the coarser mesh, do so.
    if save_yes_no
        save(['../data/sims/MC/ABM_Histograms/ABM_Histogram_MonteCarlo_BirthOption=',num2str(BirthOption),'_runs=',num2str(NumIter),'_dt=',num2str(dt),'_NBins=',num2str(NBin),'_T=150.mat'], 'MC');
    end

end
