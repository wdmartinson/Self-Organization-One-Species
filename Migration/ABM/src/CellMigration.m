%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Purpose: Simplifed agent-based code to simulate two-cell agent-based
%%% model (black and yellow cells) with only cell migration.
%%% First created: 8 March 2015
%%% Simplified/adpated: February 2022
%%% Written by: Alexandria Volkening
%%% Parameters:
%%% time = time (in days) to run simulation
%%% randomseed = seed to initiate random number generator with
%%% path = folder to save results in
%%%
%%% Sample simulation command (to run, just copy this into Matlab):
% CellMigration(10, 1, 'folderName', 'resultsName', 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cellsM, cellsX, path] = CellMigration(time, randomseed, path, name, birthOption, dt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Setting parameters
naturalMdensity = 400;  % in melanophores/mm^2
naturalXdensity = 772;  % in xanthophores/mm^2 (set to 0 for no xanthophores)
% dt = 0.1;              % time step in days for cell movement
%%% Al: dt = 0.1 in my simulations usually (time is always measured in
%%% days, so dt = 0.1 is 0.1 days)
numTimeSteps = 1/dt;    % number of time steps per day
xIC = 3;                % initial length of domain (in mm)
yIC = 3;                % initial height of domain (in mm)
center = 0;             % center of domain

%%% Setting parameters that turn on/off portions of the model
saveon = 0;             % 1 = save data and figures, 0 = do not save

%%% Setting migration parameters
rmm = 0.02;                     % Length scales for migration forces for mel-mel (in mm)
rxx = 0.011;                    % " xan-on-xan (in mm)
rxm = 0.02;                     % " xan-on-mel (in mm)
rmx = 0.02;                     % " mel-on-xan repulsion (in mm)
amx = 0.012;                    % " mel-on-xan attraction (in mm)
Rmm = 0.00124;                  % Repulsion from mel to mel (mm^2/day)
Rxx = 0.00055;                  % Repulsion from xan to xan (mm^2/day)
Rxm = 0.00274;                  % Repulsion force on mel due to xan (mm^2/day)
Rmx = 0.00226;                  % Repulsion force on xan due to mel (mm^2/day)
Amx = 0.001956;                 % Attraction force on xan due to mel (mm^2/day)
cutoff = 0.2;                   % we set forces to 0 for lengthscales >= 0.2 mm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Uncomment to plot migration forces
% figure
% d = 0:.001:0.3;
% plot(d, Rmm*exp(-d./rmm).*(d < 0.2),'LineWidth',3)
% hold on
% plot(d, -(Rmm/rmm)*exp(-d./rmm).*(d < 0.2),'LineWidth',3)
% set(gca,'FontSize',16, 'LineWidth', 2)
% legend('potential (mm^2/day)', 'velocity (mm/day)')
% box('off')
% xlabel('distance (mm)')
% pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Starting time at time_step = 2 (time = 1 is the initial condition).
%%% We measure time in days, and we save our results every save_step time
%%% steps
time_step = 2;
save_step = 1;       % results are saved every save_step time steps

%%% Making a new directory if we are saving the results
if saveon == 1
    mkdir(path)
end

%%% Seeding our random number generator
rng(randomseed);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Defining the initial cell populations; N_M and N_X house the
%%% number of melanophores and xanthophores (which is constant in these
%%% simplified simulations)
[cellsX, cellsM, N_X, N_M] = defineICPop(birthOption,xIC,yIC,time,center,save_step,naturalMdensity, naturalXdensity);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Evolution of pattern -- all in terms of time steps, rather than days,
%%% where days =
migrationPriorM = cellsM(1:N_M,:,1);
migrationPriorX = cellsX(1:N_X,:,1);
save_counter = 2;
while time_step <= time + 1

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Cell migration
    for i = 1:numTimeSteps
        [migrationTempM, migrationTempX] = EulerForward(migrationPriorM, ...
            migrationPriorX, N_M,N_X,dt,Rmm,Rxx,Rxm,Rmx,Amx,rmm,rxx,rxm,rmx,amx,xIC,yIC,cutoff);
%         [migrationTempM, migrationTempX] = StiffODESolver(migrationPriorM,...
%             migrationPriorX, N_M,N_X,dt,Rmm,Rxx,Rxm,Rmx,Amx,rmm,rxx,rxm,rmx,amx,xIC,yIC,cutoff);
        migrationPriorM = migrationTempM;
        migrationPriorX = migrationTempX;
    end
    %%% Saving every save_step days
    if mod(time_step, save_step) == 0
        cellsM(1:N_M,:,save_counter) = migrationTempM;
        cellsX(1:N_X,:,save_counter) = migrationTempX;
        save_counter = save_counter + 1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%% Finally, we step time forward by one day
    time_step = time_step + 1;

end     % end of for loop in time




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Uncomment/comment to plot/not plot cell positions
% h=figure;
% for timeOfInterest = 1:1:time+1
%
%     mel = cellsM(1:N_M,:,timeOfInterest);
%     xan = cellsX(1:N_X,:,timeOfInterest);
%
%     %%% Plotting melanophores
%     plot(mel(:,1), mel(:,2), 'o',...
%         'MarkerEdgeColor','k',...
%         'MarkerFaceColor','k',...
%         'MarkerSize', 2)
%     hold on
%     %%% Plotting xanthophores
%     plot(xan(:,1), xan(:,2), 'o',...
%         'MarkerFaceColor',[1 165/255 0],...
%         'MarkerEdgeColor',[1  165/255 0],...
%         'MarkerSize', 2)
%
%
%     title(['Time: ' num2str((timeOfInterest-1)) ' days']);
%     xlabel('fish length x (mm)')
%     ylabel('fish height y (mm)')
%     rectangle('Position',[-xIC/2 -yIC/2  xIC yIC])
%     set(gca,'XTick', -xIC/2:0.5:xIC/2, 'FontSize', 16, 'FontName', 'Helvetica', 'LineWidth', 2);
%     axis equal
%     box('off')
%     hold off
%
%     F(timeOfInterest) = getframe(h);
%
% end
% % if saveon == 1
% % save([path '/Movie_' name '_' num2str(randomseed) '.mat'], 'F');
% % end
% implay(F,2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Uncomment/comment to create/not create a heatmap with cell densities
% step = 3/34;  % in mm
% numSteps = ceil(xIC/step);
% stepBnd = -xIC/2:step:xIC/2;
% densityM = zeros(ceil(xIC/step), ceil(xIC/step));
% figure
% for timeOfInterest = 1:40:(time/dt +1)
%     for i = 1:numSteps
%         for j = 1:numSteps
%             densityM(i,j) = sum((cellsM(1:N_M,1,timeOfInterest) > stepBnd(i)).*...
%                 (cellsM(1:N_M,1,timeOfInterest) <= stepBnd(i+1)).*...
%             (cellsM(1:N_M,2,timeOfInterest) > stepBnd(j)).*...
%                 (cellsM(1:N_M,2,timeOfInterest) <= stepBnd(j+1)));
%         end
%     end
%    % densityM = densityM/((3/34)^2);
%     heatmap(-xIC/2+step/2:step:xIC/2-step/2,-xIC/2+step/2:step:xIC/2-step/2,densityM)
%     title(['Density (# cells/mm^2) at ' num2str((timeOfInterest-1)*40) ' days (bin width: ' num2str(step) ' mm)']);
%     colorbar
%     caxis([0 6])
%     F2(timeOfInterest) = getframe(gcf);
% end
% % if saveon == 1
% % save([path '/MovieHeatMap_' name '_' num2str(randomseed) '.mat'], 'F2');
% % end
% implay(F2,20)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Saving final results if saveon = 1
if saveon == 1
     save([path '/Mel_' name '_' num2str(randomseed) '.mat'], 'cellsM', 'N_M', 'N_X', 'cellsX','xIC','yIC','save_step','dt', ...
         'Rmm','rmm');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end     % end of main function



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This function defines the initial condition -- the only setting that we
%%% are using in the simplifed case is for birthOption = 100
function [cellsX, cellsM,N_X,N_M] = defineICPop(birthOption,xIC,yIC,time,center,save_step,naturalMdensity,naturalXdensity)

%%% Initializations
counterM = 1;
cellsM = -40000*ones(400,2,time/save_step+1);
cellsX = -40000*ones(400,2,time/save_step+1);
sep = 50/1000;              % cells are separated by 50/1000 mm in grid under birthOption == 0

%%% If-else for different initial conditions
if birthOption == 1   %% random distribution of cells in center 1mm by 1mm square
    N_M = round(naturalMdensity);
    N_X = 0;
    cellsM(1:N_M,:,1) = [-500/1000 + rand(N_M,1), -500/1000 + rand(N_M,1)];
    cellsX(1:N_X,:,1) = [-500/1000 + rand(N_X,1), -500/1000 + rand(N_X,1)];
elseif birthOption == 2   %% one stripe of Melanophores
    N_M = round(naturalMdensity*0.5*3);
    N_X = 0;
    cellsM(1:N_M,:,1) = [-xIC/2 + xIC*rand(N_M,1), -250/1000 + 0.5*rand(N_M,1)];
elseif birthOption == 4   %% random distribution of Xanthophore cells in center 1mm by 1mm square
    N_X = round(naturalXdensity);
    N_M = 0;
    cellsM(1:N_M,:,1) = [-500/1000 + rand(N_M,1), -500/1000 + rand(N_M,1)];
    cellsX(1:N_X,:,1) = [-500/1000 + rand(N_X,1), -500/1000 + rand(N_X,1)];
elseif birthOption == 5   %% one stripe of Xantophores
    N_X = round(naturalXdensity*0.5*3);
    N_M = 0;
    cellsX(1:N_X,:,1) = [-xIC/2 + xIC*rand(N_X,1), -250/1000 + 0.5*rand(N_X,1)];
elseif birthOption == 11 % Two stripes of melanophores

    N_M1 = round(naturalMdensity*0.5*3);
    N_M2 = round(naturalMdensity*0.5*3);
    N_X = 0;
    cellsM(1:N_M1,:,1) = [-xIC/2 + xIC*rand(N_M1,1), 0.5*rand(N_M1,1)];
    cellsM(N_M1+1:N_M1+N_M2, :, 1) = [-xIC/2 + xIC*rand(N_M2,1), -.5820 + 0.5*rand(N_M2,1)];
    N_M = N_M1 + N_M2;

elseif birthOption == 12 % Two offset rectangles of melanophores

    N_M1 = round(naturalMdensity*1*0.5);
    N_M2 = round(naturalMdensity*1*0.5);
    N_X = 0;
    cellsM(1:N_M1,:,1) = [-0.6 + 0.5*rand(N_M1,1), -0.75 + 1.0*rand(N_M1,1)];
    cellsM(N_M1+1:N_M1+N_M2,:,1) = [+0.1 + 0.5*rand(N_M2,1), -0.25 + 1.0*rand(N_M2,1)];
    N_M = N_M1 + N_M2;

elseif birthOption == 13 % Two stripes of xanthophores

    N_M = 0;
    N_X1 = round(naturalXdensity*0.36*3);
    N_X2 = round(naturalXdensity*0.36*3);
    cellsX(1:N_X1,:,1) = [-xIC/2 + xIC*rand(N_X1,1), 0.36*rand(N_X1,1)];
    cellsX(N_X1+1:N_X1+N_X2, :, 1) = [-xIC/2 + xIC*rand(N_X2,1), -.4420 + 0.36*rand(N_X2,1)];
    N_X = N_X1 + N_X2;

elseif birthOption == 14 % Two offset rectangles of xanthophores

    N_M = 0;
    N_X1 = round(naturalXdensity*1*0.5);
    N_X2 = round(naturalXdensity*1*0.5);
    cellsX(1:N_X1,:,1) = [-0.6 + 0.5*rand(N_X1,1), -0.75 + 1.0*rand(N_X1,1)];
    cellsX(N_X1+1:N_X1+N_X2,:,1) = [+0.1 + 0.5*rand(N_X2,1), -0.25 + 1.0*rand(N_X2,1)];
    N_X = N_X1 + N_X2;

else
    error('Birth Option not specified!')
end     % end of if-else on different initial conditions

% Remove the -40000 entries if there are fewer than 400 cells in the
% domain:
cellsM(N_M+1:400, :, :) = [];
cellsX(N_X+1:400, :, :) = [];

end     % end of helper function for initial conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function to implement forward Euler
function [updatedM, updatedX] = EulerForward(M0,X0,NM,NX,dt,Rmm,Rxx,Rxm,Rmx,Amx,rmm,rxx,rxm,rmx,amx,bnd_X,bnd_Y,cutoff)

%%% One Forward Euler step
[dMdt,dXdt] = XanMel_migration(M0, X0,NM,NX,Rmm,Rxx,Rxm,Rmx,Amx,rmm,rxx,rxm,rmx,amx,bnd_X,bnd_Y,cutoff);
updatedM = M0 + dt*dMdt;
updatedX = X0 + dt*dXdt;

%%% Periodic boundary conditions
updatedM(:,1) = mod(bnd_X/2+updatedM(:,1),bnd_X)-bnd_X/2;
updatedX(:,1) = mod(bnd_X/2+updatedX(:,1),bnd_X)-bnd_X/2;
updatedM(:,2) = mod(bnd_Y/2+updatedM(:,2),bnd_Y)-bnd_Y/2;
updatedX(:,2) = mod(bnd_Y/2+updatedX(:,2),bnd_Y)-bnd_Y/2;


end     % end of helper function for forward Euler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Helper function for implementing righthand side of differential
%%% equations for cell migration, called by forward Euler routine
function [rhsM, rhsX] = XanMel_migration(M, X,NM,NX,Rmm,Rxx,Rxm,Rmx,Amx,rmm,rxx,rxm,rmx,amx,lenX,lenY,cutoff)

%%% Initializations
rhsM = zeros(NM,2);
rhsX = zeros(NX,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Defining norms for periodic boundary conditions
myNorm = @(x,y)(sqrt(min((x(:,2)-y(:,2)).^2, (lenY-abs(x(:,2)-y(:,2))).^2) + ...
    min((x(:,1)-y(:,1)).^2, (lenX-abs(x(:,1)-y(:,1))).^2)));
mysignI = @(x,y)((abs(x-y)<=lenX/2).*(x-y)-(abs(x-y)>lenX/2).*(lenX-abs(x-y)).*sign(x-y));
mysignJ = @(x,y)((abs(x-y)<=lenY/2).*(x-y)-(abs(x-y)>lenY/2).*(lenY-abs(x-y)).*sign(x-y));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% These matrices hold the pairwise distances between our cells
distanceMM = pdist2(M,M,myNorm);        % distances between mel cells
distanceMX = pdist2(X,M, myNorm);       % distances between xan and mel (from xan perspective)
distanceXM = pdist2(M,X, myNorm);       % distances between xan and mel (from mel perspective)
distanceXX = pdist2(X, X, myNorm);      % distances between xan cells

%%% We begin by calculating the force strengths; each kernel is a matrix
kernelMM = q_mm(distanceMM,Rmm,rmm,cutoff);
kernelXX = q_xx(distanceXX,Rxx,rxx,cutoff);
kernelXM = q_xm(distanceXM,Rxm,rxm,cutoff);
kernelMX = q_mx(distanceMX,Rmx,rmx,Amx,amx,cutoff);

%%%% Now we calculate the direction of the force in the x and y directions
IvectorMM = pdist2(M(:,1), M(:,1), mysignI);
JvectorMM = pdist2(M(:,2), M(:,2), mysignJ);
IvectorXM = pdist2(M(:,1), X(:,1), mysignI);
JvectorXM = pdist2(M(:,2), X(:,2), mysignJ);
IvectorXX = pdist2(X(:,1),X(:,1), mysignI);
JvectorXX = pdist2(X(:,2),X(:,2), mysignJ);
IvectorMX = pdist2(X(:,1), M(:,1), mysignI);
JvectorMX = pdist2(X(:,2), M(:,2), mysignJ);

%%%% Now we calculate the norms of therse vectors
NormvectorMM = sqrt(IvectorMM.^2 + JvectorMM.^2);
NormvectorXM = sqrt(IvectorXM.^2 + JvectorXM.^2);
NormvectorXX = sqrt(IvectorXX.^2 + JvectorXX.^2);
NormvectorMX = sqrt(IvectorMX.^2 + JvectorMX.^2);

%%%% Dividing by the norms, we obtain unit vectors in the direction of the
%%%% force; Ivector contains the x coord of this vector; Jvector contains
%%%% the y coord of this vector. Note that there will be at least one place
%%%% where dividing by Normvector is division by zero (along the diagonal,
%%%% since we are calculating the distance between a cell and itself).
IvectorMM = IvectorMM./NormvectorMM;
JvectorMM = JvectorMM./NormvectorMM;
IvectorXM = IvectorXM./NormvectorXM;
JvectorXM = JvectorXM./NormvectorXM;
IvectorXX = IvectorXX./NormvectorXX;
JvectorXX = JvectorXX./NormvectorXX;
IvectorMX = IvectorMX./NormvectorMX;
JvectorMX = JvectorMX./NormvectorMX;

%%%% To account for any places we divided by zero (which would give NaN in
%%%% Matlab), we check for any NaN entries and set them to zero.
IvectorMM(isnan(IvectorMM)) = 0;
JvectorMM(isnan(JvectorMM)) = 0;
IvectorXM(isnan(IvectorXM)) = 0;
JvectorXM(isnan(JvectorXM)) = 0;
IvectorXX(isnan(IvectorXX)) = 0;
JvectorXX(isnan(JvectorXX)) = 0;
IvectorMX(isnan(IvectorMX)) = 0;
JvectorMX(isnan(JvectorMX)) = 0;

%%% Computing righthand side of our differential equations
rhsM(:,1) = sum(kernelMM.*IvectorMM,2) + sum(kernelXM.*IvectorXM,2);
rhsM(:,2) = sum(kernelMM.*JvectorMM,2) + sum(kernelXM.*JvectorXM,2);
rhsX(:,1) = sum(kernelXX.*IvectorXX,2) + sum(kernelMX.*IvectorMX,2);
rhsX(:,2) = sum(kernelXX.*JvectorXX,2) + sum(kernelMX.*JvectorMX,2);

end     % end of helper function for cell migration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Series of four helper functions called in XanMel_migration to implement
%%% individual forces between cells

%%% Short range repulsion/long range attraction of melanophores
function forcesMM = q_mm(d,R,r,cutoff)
forcesMM = -(R/r)*exp(-d./r).*(d < cutoff);
end

%%% Short range repulsion of melanophores by xanthophores
function forcesXM = q_xm(d,R,r,cutoff)
forcesXM = -(R/r)*exp(-d./r).*(d < cutoff);
end

%%% Contact attraction & short range repulsion of xanthophores by melanophores
function forcesMX = q_mx(d,R,r,A,a,cutoff)
forcesMX = -(R/r)*exp(-d./r).*(d < cutoff) + (A/a)*exp(-d./a).*(d < cutoff);
end

%%% Short range repulsion of xanthophores
function forcesXX = q_xx(d,R,r,cutoff)
forcesXX = -(R/r)*exp(-d./r).*(d < cutoff);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Series of functions for incorporating a stiff ODE solver into the
%%% agent-based framework.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [updatedM, updatedX] = StiffODESolver(M0,X0,NM,NX,dt,Rmm,Rxx,Rxm,Rmx,Amx,rmm,rxx,rxm,rmx,amx,bnd_X,bnd_Y,cutoff)
    %%% Solve the ODEs for cell migration using a stiff ODE solver
    [~, updatedW] = ode15s(@(z,W) MigrationODEs(z, W, NM,NX,Rmm,Rxx,Rxm,Rmx,Amx,rmm,rxx,rxm,rmx,amx,bnd_X,bnd_Y,cutoff), [0, dt], [M0; X0]);
    updatedM = [updatedW(end, 1:NM)', updatedW(end, NM+NX+1:2*NM+NX)'];
    updatedX = [updatedW(end, NM+1:NM+NX)', updatedW(end, 2*NM+NX+1:end)'];

    %%% Apply periodic boundary conditions
    updatedM(:,1) = mod(bnd_X/2+updatedM(:,1),bnd_X)-bnd_X/2;
    updatedX(:,1) = mod(bnd_X/2+updatedX(:,1),bnd_X)-bnd_X/2;
    updatedM(:,2) = mod(bnd_Y/2+updatedM(:,2),bnd_Y)-bnd_Y/2;
    updatedX(:,2) = mod(bnd_Y/2+updatedX(:,2),bnd_Y)-bnd_Y/2;
end

function dwdt = MigrationODEs(~, W, NM,NX,Rmm,Rxx,Rxm,Rmx,Amx,rmm,rxx,rxm,rmx,amx,bnd_X,bnd_Y,cutoff)
    %%% Extract Melanophore and Xanthophore positions. Note that ode15s
    %%% function makes W a vector, not an array, so we have to reshape the
    %%% matrices to make it compatible with Alexandria's functions:
    M = [W(1:NM), W(NM+NX+1:2*NM+NX)]; X = [W(NM+1:NM+NX), W(2*NM+NX+1:end)];

    %%% Compute the RHS of the ODEs using the function given by Alexandria
    [rhsM, rhsX] = XanMel_migration(M, X,NM,NX,Rmm,Rxx,Rxm,Rmx,Amx,rmm,rxx,rxm,rmx,amx,bnd_X,bnd_Y,cutoff);

    %%% Recombine the Melanophore and Xanthophore positions so that the
    %%% function is compatible with the syntax of ode15s. Note that you
    %%% also have to convert it to a column vector to make its shape
    %%% compatible with the original W supplied by ode15s.
    dwdt = [rhsM; rhsX];
    dwdt = dwdt(:);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%
