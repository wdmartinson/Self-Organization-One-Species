str = computer;
if str(1) == 'M'
    setenv('PATH', '/usr/local/bin:/usr/local/sbin:/usr/bin:/bin:/usr/sbin:/sbin');
end
close all;

%% Chose Ndart value to extract
n = 5; % Nvals = [10, 25, 50, 100, 150];

%% MelBox Data
% load('../../ABM/melBox_birthMoveN150_ABMexample.mat')

% load("../../ABM/melBox_BirthMoveFixed.mat");

% load('../data/sims/MelBox/MelBox_Combined_PDE_N=10_c+=Theoretical.mat')

% load('../data/sims/MelBox/MelBox_Combined_PDE_N=25_c+=Theoretical.mat')

% load('../data/sims/MelBox/MelBox_Combined_PDE_N=50_c+=Theoretical.mat')

% load('../data/sims/MelBox/MelBox_Combined_PDE_N=100_c+=Theoretical.mat')

% load('../data/sims/MelBox/MelBox_Combined_PDE_N=150_c+=Theoretical.mat')


%% OffsetRectanglesMelOnly Data
load('../../ABM/melTwoBox_birthMoveN150_ABMexample.mat')

load("../../ABM/melTwoBox_BirthMoveFixed.mat");

% load('../data/sims/OffsetRectanglesMelOnly/OffsetRectanglesMelOnly_Combined_PDE_N=10_c+=Theoretical.mat')

% load('../data/sims/OffsetRectanglesMelOnly/OffsetRectanglesMelOnly_Combined_PDE_N=25_c+=Theoretical.mat')

% load('../data/sims/OffsetRectanglesMelOnly/OffsetRectanglesMelOnly_Combined_PDE_N=50_c+=Theoretical.mat')

% load('../data/sims/OffsetRectanglesMelOnly/OffsetRectanglesMelOnly_Combined_PDE_N=100_c+=Theoretical.mat')

load('../data/sims/OffsetRectanglesMelOnly/OffsetRectanglesMelOnly_Combined_PDE_N=150_c+=Theoretical.mat')


%% Visualize Solutions
% Get mesh for pde
x = (stepBnd(2:end) + stepBnd(1:end-1))./2; % 0.1 mm spacing on the domain, stepBnd contains voxel edges

%% Total Melanophore Density Video
folder_location = '../plots/';
name = ['OffsetRectanglesMelOnly-N=',num2str(Nvals(n)),'-'];

a_name = [name, 'ABM'];
a = VideoWriter(a_name);
a.Quality = 100;
a.FrameRate = 5;

b_name = [name, 'PDE'];
b = VideoWriter(b_name);
b.Quality = 100;
b.FrameRate = 5;

c_name = [name, 'Error'];
c = VideoWriter(c_name);
c.Quality = 100;
c.FrameRate = 5;

fig1 = figure;
fig2 = figure;
fig3 = figure;
open(a);
open(b);
open(c);
for i = 1:length(t)
    figure(fig1);
    imagesc(x, x, densityM{n}(:, :, i)')
    cbar = colorbar; colormap('parula');
    caxis([0 max(max(max(CoarserMel)))])
    hold on;
    scatter(cellsM_example{i}(:,1), cellsM_example{i}(:,2), 10, 'k', 'filled')
    set(cbar, 'fontsize', 20, 'ticklabelinterpreter', 'latex')
    xlabel("$x$", 'Interpreter', 'latex', 'fontsize', 20)
    ylabel("$y$", 'Interpreter', 'latex', 'fontsize', 20)
    set(get(gca, 'title'), 'string', ['ABM, $N = ', num2str(Nvals(n)), '$, $t = ', num2str(t(i)),'$'], 'interpreter', 'latex');
    set(gca, 'fontsize', 20, 'ticklabelinterpreter', 'latex', 'YDir', 'normal')
    axis square; grid on; hold off;
    drawnow;


    figure(fig2);
    imagesc(x, x, squeeze(CoarserMel(i, :, :)))
    cbar = colorbar; colormap('parula');
    caxis([0 max(max(max(CoarserMel)))])
    hold on;
    scatter(cellsM_example{i}(:,1), cellsM_example{i}(:,2), 10, 'k', 'filled')
    set(cbar, 'fontsize', 20, 'ticklabelinterpreter', 'latex')
    xlabel("$x$", 'Interpreter', 'latex', 'fontsize', 20)
    ylabel("$y$", 'Interpreter', 'latex', 'fontsize', 20)
    set(get(gca, 'title'), 'string', ['PDE, $N = ', num2str(Nvals(n)), '$, $t = ', num2str(t(i)),'$'], 'interpreter', 'latex');
    set(gca, 'fontsize', 20, 'ticklabelinterpreter', 'latex', 'YDir', 'normal')
    axis square; grid on; hold off;
    drawnow;


    figure(fig3);
    imagesc(x, x, densityM{n}(:, :, i)' - squeeze(CoarserMel(i, :, :)))
    cbar = colorbar; colormap('cool');
    caxis([-100 250]) % For MelBox and OffsetRectanglesMelOnly cases
    hold on;
    scatter(cellsM_example{i}(:,1), cellsM_example{i}(:,2), 10, 'k', 'filled')
    set(cbar, 'fontsize', 20, 'ticklabelinterpreter', 'latex')
    xlabel("$x$", 'Interpreter', 'latex', 'fontsize', 20)
    ylabel("$y$", 'Interpreter', 'latex', 'fontsize', 20)
    set(get(gca, 'title'), 'string', ['ABM-PDE, $N = ', num2str(Nvals(n)), '$, $t = ', num2str(t(i)),'$'], 'interpreter', 'latex');
    set(gca, 'fontsize', 20, 'ticklabelinterpreter', 'latex', 'YDir', 'normal')
    axis square; grid on; hold off;
    drawnow;

     % Take the snapshots
    frame1 = getframe(fig1);
    frame2 = getframe(fig2);
    frame3 = getframe(fig3);

    writeVideo(a, frame1);
    writeVideo(b, frame2);
    writeVideo(c, frame3);

end
close(a);
close(b);
close(c);

% Convert the AVI file to an MP4:
system(['ffmpeg -n -i ',a_name,'.avi ',a_name,'.mp4']);
system(['ffmpeg -n -i ',b_name,'.avi ',b_name,'.mp4']);
system(['ffmpeg -n -i ',c_name,'.avi ',c_name,'.mp4']);


% Delete the AVI file:
system(['rm -f ', a_name,'.avi']);
system(['rm -f ', b_name,'.avi']);
system(['rm -f ', c_name,'.avi']);


% Move the MP4 to the storage folder:
system(['mv ', a_name, '.mp4 ', folder_location, a_name, '.mp4']);
system(['mv ', b_name, '.mp4 ', folder_location, b_name, '.mp4']);
system(['mv ', c_name, '.mp4 ', folder_location, c_name, '.mp4']);