%% Chose Ndart value to extract
n = 5; % Nvals = [10, 25, 50, 100, 150];

%% MelBox Data
% load("../../ABM/melBox_BirthMoveFixed.mat");

% load('../data/sims/MelBox/MelBox_Combined_PDE_N=10_c+=Theoretical.mat')

% load('../data/sims/MelBox/MelBox_Combined_PDE_N=25_c+=Theoretical.mat')

% load('../data/sims/MelBox/MelBox_Combined_PDE_N=50_c+=Theoretical.mat')

% load('../data/sims/MelBox/MelBox_Combined_PDE_N=100_c+=Theoretical.mat')

% load('../data/sims/MelBox/MelBox_Combined_PDE_N=150.mat')
% load('../data/sims/MelBox/MelBox_Combined_PDE_N=150_c+=Theoretical.mat')
% load('../data/sims/MelBox/MelBox_Combined_PDE_N=150_c+=Theoretical_r=1DFit.mat')
% load('../data/sims/MelBox/MelBox_Combined_PDE_N=150_c+=Theoretical_r=0.2.mat')
% load('../data/sims/MelBox/MelBox_Combined_PDE_N=150_c+=Theoretical_r=0.16.mat')

%% OffsetRectanglesMelOnly Data
load("../../ABM/melTwoBox_BirthMoveFixed.mat");

% load('../data/sims/OffsetRectanglesMelOnly/OffsetRectanglesMelOnly_Combined_PDE_N=10_c+=Theoretical.mat')

% load('../data/sims/OffsetRectanglesMelOnly/OffsetRectanglesMelOnly_Combined_PDE_N=25_c+=Theoretical.mat')

% load('../data/sims/OffsetRectanglesMelOnly/OffsetRectanglesMelOnly_Combined_PDE_N=50_c+=Theoretical.mat')

% load('../data/sims/OffsetRectanglesMelOnly/OffsetRectanglesMelOnly_Combined_PDE_N=100_c+=Theoretical.mat')

% load('../data/sims/OffsetRectanglesMelOnly/MelBox_Combined_PDE_N=150.mat')
load('../data/sims/OffsetRectanglesMelOnly/OffsetRectanglesMelOnly_Combined_PDE_N=150_c+=Theoretical.mat')
% load('../data/sims/OffsetRectanglesMelOnly/OffsetRectanglesMelOnly_Combined_PDE_N=150_c+=Theoretical_r=1DFit.mat')
% load('../data/sims/OffsetRectanglesMelOnly/OffsetRectanglesMelOnly_Combined_PDE_N=150_c+=Theoretical_r=0.2.mat')
% load('../data/sims/OffsetRectanglesMelOnly/OffsetRectanglesMelOnly_Combined_PDE_N=150_c+=Theoretical_r=0.16.mat')


%% Visualize Solutions
% Get mesh for pde
x = (stepBnd(2:end) + stepBnd(1:end-1))./2; % 0.1 mm spacing on the domain, stepBnd contains voxel edges

% for n = 1:length(Nvals)
% n = 4;

%% Total Melanophore Density Plots
fig1 = figure;
figure(fig1);
for i = 1:size(densityM{n}, 3)
    subplot(1, 2, 1);
    imagesc(x, x, densityM{n}(:, :, i)')
    cbar = colorbar;
    caxis([0 max(max(max(CoarserMel)))])
    set(gca, 'fontsize', 20, 'ticklabelinterpreter', 'latex')
    xlabel("$x$", 'Interpreter', 'latex', 'fontsize', 20)
    ylabel("$y$", 'Interpreter', 'latex', 'fontsize', 20)
    title(['ABM, $N = ', num2str(Nvals(n)), '$, $t = ', num2str(t(i)),'$'], 'interpreter', 'latex');
    axis square
    subplot(1, 2, 2);
    imagesc(x, x, squeeze(CoarserMel(i, :, :)))
    cbar = colorbar;
    caxis([0 max(max(max(CoarserMel)))])
    set(gca, 'fontsize', 20, 'ticklabelinterpreter', 'latex')
    xlabel("$x$", 'Interpreter', 'latex', 'fontsize', 20)
    ylabel("$y$", 'Interpreter', 'latex', 'fontsize', 20)
    title(['PDE, $N = ', num2str(Nvals(n)), '$, $t = ', num2str(t(i)),'$'], 'interpreter', 'latex');
    axis square
    drawnow;
end

%% NumberMelanophoresBorn Plot (?)
% fig2 = figure;
% figure(fig2);
% plot(0.5*(tvals(2:end) + tvals(1:end-1)), numMelBornc{n})
% title(['$N = ', num2str(Nvals(n)),'$'], 'interpreter', 'latex');
% set(gca, 'fontsize', 20, 'ticklabelinterpreter', 'latex')
% xlabel("$t$", 'Interpreter', 'latex', 'fontsize', 20)
% ylabel("No. Mel.s Born (?) ", 'Interpreter', 'latex', 'fontsize', 20)

%% Total Melanophore Count
% fig3 = figure;
% figure(fig3);
% plot(tvals, numMelc{n})
% title(['$N = ', num2str(Nvals(n)),'$'], 'interpreter', 'latex');
% set(gca, 'fontsize', 20, 'ticklabelinterpreter', 'latex')
% xlabel("$t$", 'Interpreter', 'latex', 'fontsize', 20)
% ylabel("Total Mel. Count", 'Interpreter', 'latex', 'fontsize', 20)

%% Mean Number1DM -- All of these are zero matrices!

%% NumMelcStd (?)
% fig4 = figure;
% figure(fig4);
% plot(tvals, numMelcStd{n})
% title(['$N = ', num2str(Nvals(n)),'$'], 'interpreter', 'latex');
% set(gca, 'fontsize', 20, 'ticklabelinterpreter', 'latex')
% xlabel("$t$", 'Interpreter', 'latex', 'fontsize', 20)
% ylabel("StDev Mel Count", 'Interpreter', 'latex', 'fontsize', 20)

%% spreadAbs (???)

% end % for n = 1:length(Nvals)