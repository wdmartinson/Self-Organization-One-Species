%% Chose Ndart value to extract
data = 'MelBox'; % 'MelBox', 'OffsetRectanglesMelOnly'
n = 5; % Nvals = [10, 25, 50, 100, 150];


%% MelBox Data
if strcmp(data, 'MelBox')
    load('../../ABM/melBox_birthMoveN150_ABMexample.mat')
    
    load("../../ABM/melBox_BirthMoveFixed.mat");
end

if n == 1 && strcmp(data, 'MelBox')
    load('../data/sims/MelBox/MelBox_Combined_PDE_N=10_c+=Theoretical.mat')
elseif n == 2 && strcmp(data, 'MelBox')
    load('../data/sims/MelBox/MelBox_Combined_PDE_N=25_c+=Theoretical.mat')
elseif n == 3 && strcmp(data, 'MelBox')
    load('../data/sims/MelBox/MelBox_Combined_PDE_N=50_c+=Theoretical.mat')
elseif n == 4 && strcmp(data, 'MelBox')
    load('../data/sims/MelBox/MelBox_Combined_PDE_N=100_c+=Theoretical.mat')
elseif n == 5 && strcmp(data, 'MelBox')
    load('../data/sims/MelBox/MelBox_Combined_PDE_N=150_c+=Theoretical.mat')
%     load('../data/sims/MelBox/MelBox_Combined_PDE_N=150.mat')
end

%% OffsetRectanglesMelOnly Data
if strcmp(data, 'OffsetRectanglesMelOnly')
    load('../../ABM/melTwoBox_birthMoveN150_ABMexample.mat')
    
    load("../../ABM/melTwoBox_BirthMoveFixed.mat");
end

if n == 1 && strcmp(data, 'OffsetRectanglesMelOnly')
    load('../data/sims/OffsetRectanglesMelOnly/OffsetRectanglesMelOnly_Combined_PDE_N=10_c+=Theoretical.mat')
elseif n == 2 && strcmp(data, 'OffsetRectanglesMelOnly')
    load('../data/sims/OffsetRectanglesMelOnly/OffsetRectanglesMelOnly_Combined_PDE_N=25_c+=Theoretical.mat')
elseif n == 3 && strcmp(data, 'OffsetRectanglesMelOnly')
    load('../data/sims/OffsetRectanglesMelOnly/OffsetRectanglesMelOnly_Combined_PDE_N=50_c+=Theoretical.mat')
elseif n == 4 && strcmp(data, 'OffsetRectanglesMelOnly')
    load('../data/sims/OffsetRectanglesMelOnly/OffsetRectanglesMelOnly_Combined_PDE_N=100_c+=Theoretical.mat')
elseif n == 5 && strcmp(data, 'OffsetRectanglesMelOnly')
    load('../data/sims/OffsetRectanglesMelOnly/OffsetRectanglesMelOnly_Combined_PDE_N=150_c+=Theoretical.mat')
end

%% Visualize Solutions
% Get mesh for pde
x = (stepBnd(2:end) + stepBnd(1:end-1))./2; % 0.1 mm spacing on the domain, stepBnd contains voxel edges

%% Total Melanophore Density Plots
% t_view = [11, 41, 71, 101];
% for i = 1:length(t_view)
%     fig1 = figure;
%     figure(fig1);
%     imagesc(x, x, densityM{n}(:, :, t_view(i))')
%     cbar = colorbar; colormap('parula');
%     caxis([0 max(max(max(CoarserMel)))])
%     hold on;
%     scatter(cellsM_example{t_view(i)}(:,1), cellsM_example{t_view(i)}(:,2), 10, 'k', 'filled')
%     set(gca, 'fontsize', 20, 'ticklabelinterpreter', 'latex', 'YDir', 'normal')
%     set(cbar, 'fontsize', 20, 'ticklabelinterpreter', 'latex')
%     xlabel("$x$", 'Interpreter', 'latex', 'fontsize', 20)
%     ylabel("$y$", 'Interpreter', 'latex', 'fontsize', 20)
%     title(['ABM, $N = ', num2str(Nvals(n)), '$, $t = ', num2str(t(t_view(i))),'$'], 'interpreter', 'latex');
%     axis square
%     saveas(fig1, ['../plots/EnsembleAverageABM-MelBox-Fittedc+-Ndarts=',num2str(Nvals(n)),'-t=',num2str(t(t_view(i)))], 'fig')
% 
%     fig2 = figure;
%     figure(fig2);
%     imagesc(x, x, squeeze(CoarserMel(t_view(i), :, :)))
%     cbar = colorbar; colormap('parula');
%     caxis([0 max(max(max(CoarserMel)))])
%     hold on;
%     scatter(cellsM_example{t_view(i)}(:,1), cellsM_example{t_view(i)}(:,2), 10, 'k', 'filled')
%     set(cbar, 'fontsize', 20, 'ticklabelinterpreter', 'latex')
%     set(gca, 'fontsize', 20, 'ticklabelinterpreter', 'latex', 'YDir', 'normal')
%     xlabel("$x$", 'Interpreter', 'latex', 'fontsize', 20)
%     ylabel("$y$", 'Interpreter', 'latex', 'fontsize', 20)
%     title(['PDE, $N = ', num2str(Nvals(n)), '$, $t = ', num2str(t(t_view(i))),'$'], 'interpreter', 'latex');
%     axis square
%     saveas(fig2, ['../plots/PDE-MelBox-Fittedc+-Ndarts=',num2str(Nvals(n)),'-t=',num2str(t(t_view(i)))], 'fig')
% 
%     fig3 = figure;
%     figure(fig3);
%     imagesc(x, x, densityM{n}(:, :, t_view(i))' - squeeze(CoarserMel(t_view(i), :, :)))
%     cbar = colorbar; colormap('cool');
%     caxis([-100 250]) % For MelBox and OffsetRectanglesMelOnly cases
%     hold on;
%     scatter(cellsM_example{t_view(i)}(:,1), cellsM_example{t_view(i)}(:,2), 10, 'k', 'filled')
%     set(gca, 'fontsize', 20, 'ticklabelinterpreter', 'latex', 'YDir', 'normal')
%     set(cbar, 'fontsize', 20, 'ticklabelinterpreter', 'latex')
%     xlabel("$x$", 'Interpreter', 'latex', 'fontsize', 20)
%     ylabel("$y$", 'Interpreter', 'latex', 'fontsize', 20)
%     title(['ABM-PDE, $N = ', num2str(Nvals(n)), '$, $t = ', num2str(t(t_view(i))),'$'], 'interpreter', 'latex');
%     axis square
%     saveas(fig3, ['../plots/Error-MelBox-Fittedc+-Ndarts=',num2str(Nvals(n)),'-t=',num2str(t(t_view(i)))], 'fig')
% 
% end


%% Total Number of Melanophores over time
% fig4 = figure;
% colors_dark  = [0.0000 0.4470 0.7410 ; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880];
% colors_light = [0.6588 0.9098 1.0000 ; 0.9804 0.7843 0.5569; 1.0000 0.9725 0.6000; 0.8588 0.6745 0.9804; 0.7255 1.0000 0.5412];
% for n = 1:length(Nvals)
% 
%     if n == 1 && strcmp(data, 'MelBox')
%     load('../data/sims/MelBox/MelBox_Combined_PDE_N=10_c+=Theoretical.mat')
%     elseif n == 2 && strcmp(data, 'MelBox')
%         load('../data/sims/MelBox/MelBox_Combined_PDE_N=25_c+=Theoretical.mat')
%     elseif n == 3 && strcmp(data, 'MelBox')
%         load('../data/sims/MelBox/MelBox_Combined_PDE_N=50_c+=Theoretical.mat')
%     elseif n == 4 && strcmp(data, 'MelBox')
%         load('../data/sims/MelBox/MelBox_Combined_PDE_N=100_c+=Theoretical.mat')
%     elseif n == 5 && strcmp(data, 'MelBox')
%         load('../data/sims/MelBox/MelBox_Combined_PDE_N=150_c+=Theoretical.mat')
% %         load('../data/sims/MelBox/MelBox_Combined_PDE_N=150.mat')
%     end
% 
%     if n == 1 && strcmp(data, 'OffsetRectanglesMelOnly')
%     load('../data/sims/OffsetRectanglesMelOnly/OffsetRectanglesMelOnly_Combined_PDE_N=10_c+=Theoretical.mat')
%     elseif n == 2 && strcmp(data, 'OffsetRectanglesMelOnly')
%         load('../data/sims/OffsetRectanglesMelOnly/OffsetRectanglesMelOnly_Combined_PDE_N=25_c+=Theoretical.mat')
%     elseif n == 3 && strcmp(data, 'OffsetRectanglesMelOnly')
%         load('../data/sims/OffsetRectanglesMelOnly/OffsetRectanglesMelOnly_Combined_PDE_N=50_c+=Theoretical.mat')
%     elseif n == 4 && strcmp(data, 'OffsetRectanglesMelOnly')
%         load('../data/sims/OffsetRectanglesMelOnly/OffsetRectanglesMelOnly_Combined_PDE_N=100_c+=Theoretical.mat')
%     elseif n == 5 && strcmp(data, 'OffsetRectanglesMelOnly')
%         load('../data/sims/OffsetRectanglesMelOnly/OffsetRectanglesMelOnly_Combined_PDE_N=150_c+=Theoretical.mat')
%     end
% 
%     MassMelPDE = zeros(length(t), 1);
%     for i = 1:length(t)
%         MassMelPDE(i) = (x(2)-x(1))^2*sum(squeeze(CoarserMel(i, :, :)), 'all');
%     end
%     y = numMelc{n};
%     std_dev = numMelcStd{n};
%     curve1 = y + std_dev;
%     curve2 = y - std_dev;
%     t2 = [t, fliplr(t)];
%     inBetween = [curve1, fliplr(curve2)];
%     
% %     fill(t2, inBetween, [0.6118    0.6118    0.6118], 'HandleVisibility', 'off');
% %     hold on;
% %     plot(t, numMelc{n}, 'k', 'linewidth', 1, 'displayname', ['ABM, $N = ', num2str(Nvals(n)), '$']);
% %     plot(t, MassMelPDE, '--r', 'linewidth', 1, 'displayname', 'PDE');
%     fill(t2, inBetween, colors_light(n, :), 'HandleVisibility', 'off');
%     hold on;
%     plot(t, numMelc{n}, 'color', colors_dark(n, :), 'linewidth', 1, 'displayname', ['ABM, $N = ', num2str(Nvals(n)), '$']);
%     plot(t, MassMelPDE, 'color', colors_dark(n, :), 'linestyle', '--', 'linewidth', 1, 'displayname', ['PDE, $N = ', num2str(Nvals(n)), '$']);
% end % for n = 1:length(Nvals)
% legend('location', 'northwest', 'interpreter', 'latex');
% xlabel("$t$", 'Interpreter', 'latex', 'fontsize', 20)
% ylabel("No. Melanophores", 'Interpreter', 'latex', 'fontsize', 20)
% set(gca, 'fontsize', 20, 'ticklabelinterpreter', 'latex', 'YDir', 'normal')
% axis square;
% box on;
% grid on;
% % saveas(fig4, ['../plots/TotalMass-MelBox-Ndarts=',num2str(Nvals(n))], 'fig');
% saveas(fig4, '../plots/TotalMass-MelBox-All-Ndarts', 'fig');

%% Spread of cells over time

fig5 = figure;
colors_dark  = [0.0000 0.4470 0.7410 ; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880];
colors_light = [0.6588 0.9098 1.0000 ; 0.9804 0.7843 0.5569; 1.0000 0.9725 0.6000; 0.8588 0.6745 0.9804; 0.7255 1.0000 0.5412];
for n = 1:length(Nvals)
    
    MeanSpreadABM  = mean(spreadAbs{n}, 1);
    StDevSpreadABM = std(spreadAbs{n}, 0, 1);

    curve1 = MeanSpreadABM + StDevSpreadABM;
    curve2 = MeanSpreadABM - StDevSpreadABM;
    t2 = [t, fliplr(t)];
    inBetween = [curve1, fliplr(curve2)];
    
    
%     fill(t2, inBetween, [0.6118    0.6118    0.6118], 'HandleVisibility', 'off');
%     hold on;
%     plot(t, numMelc{n}, 'k', 'linewidth', 1, 'displayname', ['ABM, $N = ', num2str(Nvals(n)), '$']);
%     plot(t, MassMelPDE, '--r', 'linewidth', 1, 'displayname', 'PDE');
    fill(t2, inBetween, colors_light(n, :), 'HandleVisibility', 'off');
    hold on;
    plot(t, MeanSpreadABM, 'color', colors_dark(n, :), 'linewidth', 1, 'displayname', ['ABM, $N = ', num2str(Nvals(n)), '$']);
end % for n = 1:length(Nvals)
for n = 1:length(Nvals)
        if n == 1 && strcmp(data, 'MelBox')
        load('../data/sims/MelBox/MelBox_Combined_PDE_N=10_c+=Theoretical.mat')
    elseif n == 2 && strcmp(data, 'MelBox')
        load('../data/sims/MelBox/MelBox_Combined_PDE_N=25_c+=Theoretical.mat')
    elseif n == 3 && strcmp(data, 'MelBox')
        load('../data/sims/MelBox/MelBox_Combined_PDE_N=50_c+=Theoretical.mat')
    elseif n == 4 && strcmp(data, 'MelBox')
        load('../data/sims/MelBox/MelBox_Combined_PDE_N=100_c+=Theoretical.mat')
    elseif n == 5 && strcmp(data, 'MelBox')
        load('../data/sims/MelBox/MelBox_Combined_PDE_N=150_c+=Theoretical.mat')
%         load('../data/sims/MelBox/MelBox_Combined_PDE_N=150.mat')
    end

    if n == 1 && strcmp(data, 'OffsetRectanglesMelOnly')
    load('../data/sims/OffsetRectanglesMelOnly/OffsetRectanglesMelOnly_Combined_PDE_N=10_c+=Theoretical.mat')
    elseif n == 2 && strcmp(data, 'OffsetRectanglesMelOnly')
        load('../data/sims/OffsetRectanglesMelOnly/OffsetRectanglesMelOnly_Combined_PDE_N=25_c+=Theoretical.mat')
    elseif n == 3 && strcmp(data, 'OffsetRectanglesMelOnly')
        load('../data/sims/OffsetRectanglesMelOnly/OffsetRectanglesMelOnly_Combined_PDE_N=50_c+=Theoretical.mat')
    elseif n == 4 && strcmp(data, 'OffsetRectanglesMelOnly')
        load('../data/sims/OffsetRectanglesMelOnly/OffsetRectanglesMelOnly_Combined_PDE_N=100_c+=Theoretical.mat')
    elseif n == 5 && strcmp(data, 'OffsetRectanglesMelOnly')
        load('../data/sims/OffsetRectanglesMelOnly/OffsetRectanglesMelOnly_Combined_PDE_N=150_c+=Theoretical.mat')
    end
    X = linspace(-1.5, 1.5, 241);
%     X = (X(2:end)+X(1:end-1))./2;
    SpreadPDE = zeros(length(t), 1);
    for i = 1:length(t)
        tmp = squeeze(Mel(i, (length(X)-1)./2, :));
        tmp(tmp < 1e-7) = 0; % Uncomment if you don't want to correct for using single precision numbers.
        SpreadPDE(i) = abs( X( find( tmp, 1, 'first' ) ) );
    end
    plot(t, SpreadPDE, 'color', colors_dark(n, :), 'linestyle', '--', 'linewidth', 1, 'displayname', ['PDE, $N = ', num2str(Nvals(n)), '$']);
end % for n = 1:length(Nvals)

% legend('location', 'northwest', 'interpreter', 'latex');
xlabel("$t$", 'Interpreter', 'latex', 'fontsize', 20)
ylabel("Distance from Origin (mm)", 'Interpreter', 'latex', 'fontsize', 20)
set(gca, 'fontsize', 20, 'ticklabelinterpreter', 'latex', 'YDir', 'normal')
axis square;
box on;
grid on;
% saveas(fig5, '../plots/RadiusOfSupport-MelBox-All-Ndarts', 'fig');


