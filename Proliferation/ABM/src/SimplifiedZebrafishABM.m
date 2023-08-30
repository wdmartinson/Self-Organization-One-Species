%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Purpose: Comprehensive program to simulate an ABM as in "Modelling 
%%% stripe formation in zebrafish: an agent-based approach" by Alexandria
%%% Volkening and Bjorn Sandstede in Journal of the Royal Society Interface
%%% 12(112):20150812, 2015. This program was adapted and made more flexible
%%% in association with "Linking agent-based and continuum models of cell
%%% birth and migration" by Jose Carrillo, W. Duncan Martinson, Markus 
%%% Schmidtchen, Chandrasekhar Venkataraman, and Alexandria Volkening, 
%%% submitted 2023. This program simulates pattern formation in one or two
%%% dimensions, with options to include cell birth, death, and/or
%%% migration, as well as one or two types of cells (melanophores and dense 
%%% xanthophores). It includes the option to select short- and/or long-
%%% range signals determined by model parameters. We note that a select
%%% portion of these options are considered in the manuscript "Linking 
%%% agent-based and continuum models of cell birth and migration".
%%% 
%%% If you use this code, please cite the references above and email
%%% alexandria.volkening@gmail.com.
%%%
%%% Written by Alexandria Volkening (2023)
%%%
%%% Parameters:
%%%% name = name to give the saved data
%%%% time = time in days to run a simulation
%%%% nrand = number of simulations to run
%%%% initFlag = initial condition (IC) specified
%%%%%%% 'singleMel' = one mel at the center of the domain, no xan
%%%%%%% 'singleXan') = one xan at the center of the domain, no mel
%%%%%%% 'twoMel' = one mel at [0.6,0] and one at [-0.8,0], no xan
%%%%%%% 'twoXan' = one xan at [0.6,0] and one at [-0.8,0], no mel
%%%%%%% 'twoBoxMel' = offset rectangles IC of mel in paper, no xan
%%%%%%% 'twoBoxXan' = offset rectangles IC of xan in paper, no mel
%%%%%%% 'twoStripeMel' = two stripe IC of mel in paper, no xan
%%%%%%% 'twoStripeXan' = two stripe IC of xan in paper, no mel
%%%%%%% 'melXanPair' = one mel and one xan at domain center (0.05 mm apart)
%%%%%%% 'rand' = 10 mel and 10 xan randomly placed (same IC for all sim)
%%%
%%% Suggested simulation (fast):
%%% resultsPath = SimplifiedZebrafishABM('example_result', 100, 1, 'singleMel')
%%%
%%% Suggested simulation (longer):
%%% resultsPath = SimplifiedZebrafishABM('example_result', 150, 1000, 'twoBoxMel')
%%%
%%%
%%% Copyright 2023 Jose Carrillo, W. Duncan Martinson, Markus Schmidtchen, 
%%% Chandrasekhar Venkataraman, and Alexandria Volkening
%%% Licensed under the Apache License, Version 2.0 (the "License"); you may
%%% not use this file except in compliance with the License. You may obtain
%%% a copy of the License at http://www.apache.org/licenses/LICENSE-2.0.
%%% Unless required by applicable law or agreed to in writing, software
%%% distributed under the License is distributed on an "AS IS" BASIS,
%%% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
%%% implied. See the License for the specific language governing
%%% permissions and limitations under the License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [path] = SimplifiedZebrafishABM(name, time, nrand, initFlag)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Flags for selecting options that implement different models
popSwitch = 1;              % 1 if only a population of melanophores, 2 if only xanthophores, 3 if both cell types
overSwitch = 0;             % 0 if overcrowding condition applied, 1 otherwise
longSwitch = 1;             % 0 if long-range signals included, 1 otherwise
randSwitch = 0;             % 1 if random birth is included for all cell types active (depends on popSwitch), 0 otherwise
deathSwitch = 0;            % 1 if cell death is included for all cell types active (depends on popSwitch), 0 otherwise
birthSwitch = 1;            % 1 if cell birth is included for all cell types active, 0 otherwise
moveSwitch = 1;             % 1 if cell movement is included for all cell types active, 0 otherwise
%%% If deathSwitch = 1 but longSwitch = 1, only local death signals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Setting model parameters: Cell migration
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
dt = 0.1;                       % time step for ODE (birth and death time step is always 1 day)
numTimeSteps = 1/dt;            % number of ODE time steps per day 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Setting model parameters: Length and time scales
xBnd = 3;                   % initial length of domain (in mm)
yBnd = 3;                   % initial height of domain (in mm); 0 if a line
domainRad = xBnd/2;         % end points of the domain are -domainRad and domainRad
omegaLoc = 0.075;           % length scale (in mm) over which cells interact (radius of ball capturing cell interactions)
omegaLong = 0.318;          % length scale (in mm) over which cells interact long range (inner radius of annulus capturing cell interactions)
omegaWidth = 0.025;         % width (in mm) of annulus for long-range cell interactions (related to omegaLong)
omegaRand = 0.1;            % length scale (in mm) for random birth of cells at low density
N = 100;                    % number of locations selected (for each cell type) uniformly at random in the domain per day (to be evaluated for cell birth)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Setting model parameters: Melanophore birth
alpha = 1;                          % short-range signals for melanophore birth
beta = 3.5;                         % long-range signals for melanophore birth
eta = 6;                            % cap on max number of melanophores that can be in omegaLoc before it becomes too overcrowed for cell birth
eta = eta + 100000*overSwitch;      % effectively turns off max-density constraint if overSwitch = 1
melRand = 0.03;                     % chance of random melanophore birth when no cells in omegaRand
melRand = melRand*randSwitch;       % if randSwitch = 0, turns off random birth
melRand = melRand*(popSwitch ~= 2); % if no melanophores included, turns off random melanophore birth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Setting model parameters: Xanthophore birth
phi = 1.3;                          % short-range signals for xanthophore birth
psi = 0.7;                          % long-range signals for xanthophore birth
kappa = 10;                         % cap on max number of xanthophores that can be in omegaLoc before overcrowding
kappa = kappa + 100000*overSwitch;  % effectively turns off max-density constraint if overSwitch = 1
xanRand = 0.005;                    % chance of random xanthophore birth when no cells in omegaRand
xanRand = xanRand*randSwitch;       % if randSwitch = 0, turns off random birth
xanRand = xanRand*(popSwitch ~= 1); % if no xanthophores included, turn off random xanthophore birth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Setting model parameters: Cell death
xi = 1.2;                                               % long-range signals for melanophore death
melProb = 0.03;                                         % chance of melanophore death due to long-range interactions
melProb = melProb*deathSwitch*(longSwitch == 0);        % if deathSwitch = 0 or longSwitch = 1, sets probability of melanophore death due to long-range interactions to 0

%%% Note: There are also short-range interactions for cell death (black and
%%% gold cells compete locally). The parameters are all 1 in these
%%% settings, and we turn this off when deathSwitch = 0 later on in the
%%% code below.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Setting other house-keeping parameters
saveon = 1;             % 1 = save data and figures, 0 = do not save
path = 'Results';       % folder to save results in if saveon = 1
if saveon == 1
    %%% If saving data, we create a folder to save in
    mkdir(path);
end

%%% Starting time at time_step = 2. Time = 1 is the initial condition. Time
%%% is measured in days (one time step is one day). periodicBC implements 
%%% periodic boundary conditions.
time_step = 2;
periodicBC = @(x,y)(sqrt(min((x(:,1)-y(:,1)).^2, (xBnd - abs(x(:,1) - y(:,1))).^2) + ...
    min((x(:,2)-y(:,2)).^2, (yBnd - abs(x(:,2) - y(:,2))).^2)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Defining the initial population of cells and initializing matrices
%%%
%%% cellsM is a cell array; cellsM{i}{j}(:,:) corresponds to the
%%% cells present in random simulation i at time j.
%%% numMel(a,b) refers to the number of cells present at time a for
%%% simulation b. cellsX and numXan are structured in the same way.
[cellsM, cellsX, numMel, numXan] = initialCondition(nrand, time, initFlag, xBnd, yBnd);      
numMelBorn = zeros(time,nrand);     % tracking number of cells born per time step
numXanBorn = zeros(time,nrand);
numMelDead = zeros(time,nrand);     % tracking number of cells dying per time step
numXanDead = zeros(time,nrand);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Main loop implementing nrand simulations of our simple model
for randii = 1:nrand

    %%% Looping through time for one model simulations
    while time_step <= time + 1
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% For ease, we use temporary matrices to refer to the cell
        %%% positions and numbers at the current time.
        N_M = numMel(time_step-1, randii);              % number of mel in the domain
        N_X = numXan(time_step-1, randii);              % number of xan in the domain
        tempM = cellsM{randii}{time_step-1};            % current positions of melanophores
        tempX = cellsX{randii}{time_step-1};            % current positions of xanthophores
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Implementing cell migration
        migrationPriorM = tempM;
        migrationPriorX = tempX;
        for i = 1:numTimeSteps
            [migrationTempM, migrationTempX] = EulerForward(migrationPriorM, ...
                migrationPriorX, N_M,N_X,dt,Rmm,Rxx,Rxm,Rmx,Amx,rmm,rxx,rxm,rmx,amx,xBnd,yBnd,cutoff,moveSwitch);
            migrationPriorM = migrationTempM;
            migrationPriorX = migrationTempX;
        end
        tempM = migrationTempM;
        tempX = migrationTempX;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Calculating cell-cell distances and selecting cells in
        %%% short- and long-range neighborhoods for later use in cell birth
        %%%
        %%% Choosing N locations (uniformly at random in the domain) to
        %%% evaluate for possible melanophore birth and another N for
        %%% possible xanthophore birth
        locationsPos = [-domainRad+xBnd*rand(2*N,1),-yBnd/2 + yBnd*rand(2*N,1)];
        
        %%% Calculating the distance between possible birth locations and
        %%% the current cell population using periodic boundary conditions
        distanceFromM = pdist2(locationsPos, tempM, periodicBC);
        distanceFromX = pdist2(locationsPos, tempX, periodicBC);
        
        %%% Counting the number of cells in a local ball
        sumBallM = sum((distanceFromM <= omegaLoc),2);
        sumBallX = sum((distanceFromX <= omegaLoc),2);
        
        %%% Counting the number of cells in a slightly larger local ball
        %%% for possible random cell birth
        sumRandM = sum((distanceFromM <= omegaRand),2);
        sumRandX = sum((distanceFromX <= omegaRand),2);
        
        %%% Counting the number of cells a in a long-range disk
        sumAnnulusM = sum((distanceFromM >= omegaLong).*(distanceFromM <= (omegaLong+omegaWidth)),2);
        sumAnnulusX = sum((distanceFromX >= omegaLong).*(distanceFromX <= (omegaLong+omegaWidth)),2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Birth of melanophores (black cells)
        %%%
        %%% Positions are approved for birth depending on short-range
        %%% signals, overcrowding conditions, and long-range signals.
        %%% Switch variables have already removed any signals indicated
        %%% (e.g., we set eta = a large number if overcrowding is off), and
        %%% longSwitch is 1 if long-range signals off, 0 otherwise
        %%%
        %%% Mborn(i) = 1 if the ith location will have birth, 0 otherwise
        Mborn = logical((sumBallM(1:N) > alpha*sumBallX(1:N)));                             % short-range condition 
        Mborn = logical(Mborn.*((sumBallM(1:N) + sumBallX(1:N)) < eta));                    % overcrowding condition
        Mborn = logical(Mborn.*(sumAnnulusX(1:N) + 100000*longSwitch > beta*sumAnnulusM(1:N)));    % long-range condition
        randomMborn = double((sumRandM(1:N) + sumRandX(1:N) == 0).*(binornd(1,melRand, [N,1])).*(sumAnnulusX(1:N) + 100000*longSwitch > beta*sumAnnulusM(1:N)));    
        %%% above is random birth constrained on long-range dynamics
        Mborn = logical(birthSwitch*(double(Mborn) + randomMborn));     % combining all of our birth
        numMborn = sum(Mborn);              % number of new melanophores born
        newM = locationsPos(1:N,:); 	    % gathering possible melanophore locations
        newM = newM(Mborn,:);               % gathering locations selected for melanophore birth
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Birth of xanthophores (gold cells)
        %%%
        %%% See above for implementation details (similar to melanophores)
        Xborn = logical((sumBallX(N+1:end) > phi*sumBallM(N+1:end)));                       % short-range condition
        Xborn = logical(Xborn.*((sumBallX(N+1:end) + sumBallM(N+1:end)) < kappa));          % overcrowding condition 
        Xborn = logical(Xborn.*(sumAnnulusM(N+1:end) + 100000*longSwitch > psi*sumAnnulusX(N+1:end)));      % long-range condition 
        randomXborn = double((sumRandM(N+1:end) + sumRandX(N+1:end) == 0).*(binornd(1,xanRand, [N,1])).*(sumAnnulusM(N+1:end) + 100000*longSwitch > psi*sumAnnulusX(N+1:end)));   
        %%% above is random birth constrained on long-range dynamics
        Xborn = logical(birthSwitch*(double(Xborn) + randomXborn));     % combining all of our birth
        numXborn = sum(Xborn);              % number of new xanthophores born
        newX = locationsPos(N+1:end,:);     % gathering possible xanthophore locations
        newX = newX(Xborn,:);               % gathering locations selected for xanthophore birth
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Calculating cell-cell distances and selecting cells in
        %%% short- and long-range neighborhoods for later use in cell death
        distanceMX = pdist2(tempM, tempX, periodicBC);      % distance between mel and xan
        distanceXX = pdist2(tempX, tempX, periodicBC);      % distance between xan
        distanceMM = pdist2(tempM, tempM, periodicBC);      % distance between mel
        distanceMXforM = distanceMX;        % MX distance from melanophores' perspective
        distanceMXforX = distanceMX';       % MX distance from xanthophores' perspective
        
        %%% Counting the number of cells a in a local ball
        closeXtoM = sum((distanceMXforM <= omegaLoc),2);       % = 1 when the measured mel-xan distance implies the xan is in the disc centered at the mel
        closeMtoX = sum((distanceMXforX <= omegaLoc),2);       % = 1 when the measured mel-xan distance implies the mel is in the disc centered at the xan
        closeXtoX = sum((distanceXX <= omegaLoc),2);
        closeMtoM = sum((distanceMM <= omegaLoc),2);
        
        %%% Counting the number of cells a in a long-range disk
        farMX = (distanceMXforM >= omegaLong).*(distanceMXforM <= omegaLong+omegaWidth);      % = 1 for each mel-xan pair that is in the annulus (0 otherwise)
        farMM = (distanceMM >= omegaLong).*(distanceMM <= omegaLong+omegaWidth);              % = 1 for each mel-mel pair that is in the disc (0 otherwise)
        farXtoM = sum(farMX,2);         % counting number of xanthophores that are far way from melanophores
        farMtoM = sum(farMM,2);         % counting number of mel far away from mel
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Death of melanophores and xanthophores
        deadX = deathSwitch*logical((closeMtoX > closeXtoX));
        bernoulli = binornd(1,melProb, [length(N_M),1]);
        deadM = deathSwitch*logical(double((closeXtoM > closeMtoM)) + double((farMtoM > xi*farXtoM)).*bernoulli);  

        %%% For cases with only one cell type, we have place a single cell
        %%% of the opposite type at [-40000,0] just to make our matrices
        %%% all work and be nonempty; we keep this cell from dying by
        %%% requiring that cells be in the domain in order to die
        inDomainM = (tempM(:,1) >= -domainRad);    % cells can only die if they are in the domain
        inDomainX = (tempX(:,1) >= -domainRad);    % cells can only die if they are in the domain
        deadM = deadM.*inDomainM;
        deadX = deadX.*inDomainX;
        numXdead = sum(deadX);          % number of xanthophores dying
        numMdead = sum(deadM);          % number of melanophores dying
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Updating our cells due to birth and death for the next time 
        cellsM{randii}{time_step} = [tempM(~deadM,:); newM];
        cellsX{randii}{time_step} = [tempX(~deadX,:); newX];
        numMel(time_step, randii) = N_M + numMborn - numMdead;
        numXan(time_step, randii) = N_X + numXborn - numXdead;
        numMelBorn(time_step-1,randii) = numMborn;
        numXanBorn(time_step-1,randii) = numXborn;
        numMelDead(time_step-1,randii) = numMdead;
        numXanDead(time_step-1,randii) = numXdead;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Finally, we step time forward by one day for this simulation
        time_step = time_step + 1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end     % end of loop over time for a given simulation
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Resetting time for our next simulation
    time_step = 2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end     % end of loop over number of random simulations



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plotting mean number of cells in time. Adjust this code based on cell
%%% type.
meanNumMel = mean(numMel,2);
stdNumMel = std(numMel,[],2);
figure
plot(0:time, meanNumMel, 'k', 'LineWidth', 3)
hold on
plot(0:time, meanNumMel + stdNumMel, 'b', 'LineWidth', 2)
plot(0:time, meanNumMel - stdNumMel, 'b', 'LineWidth', 2)
legend({'mean number of melanophores', '+/- standard deviation'})
xlabel('time steps (days)')
ylabel('number of melanophores')
box('off')
set(gca,'FontSize', 16, 'FontName', 'Helvetica', 'LineWidth', 2)
title(['Ball radius = ' num2str(omegaLoc) ' mm and ' num2str(N) ' positions evaluated per day']);
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Creating a heatmap with mean cell densities in time. Adjust this code
%%% based on cell type.
step = 0.1;                                         % grid length in mm
numSteps = ceil(xBnd/step);                         % number of squares in x grid
numStepsY = ceil(yBnd/step);                        % number of squares in y grid
stepBnd = -domainRad:step:domainRad;                % left endpoints of grid
stepBndY = -yBnd/2:step:yBnd/2;                     % y grid
timeStep = 10;                                      % calculated each timeStep steps
numTsteps = ceil(time/timeStep);                    % total number of times evaluated
times = 1:timeStep:time;                            % vector of times to be evaluated at
densityM = zeros(numSteps, max(1,numStepsY), numTsteps);   % matrix to hold mean cell density
for j = 1:numTsteps
    timeOfInterest = times(j);
    for randii = 1:nrand
        for i = 1:numSteps
            for k = 1:numStepsY
                %%% Counting up the cells in each grid square across each
                %%% simulation and dividing by the grid size to express in
                %%% cells/mm
                densityM(i,k,j) = densityM(i,k,j) + (sum((cellsM{randii}{timeOfInterest}(:,1) > stepBnd(i)).*...
                    (cellsM{randii}{timeOfInterest}(:,1) <= stepBnd(i+1)).*...
                    (cellsM{randii}{timeOfInterest}(:,2) > stepBndY(k)).*...
                    (cellsM{randii}{timeOfInterest}(:,2) <= stepBndY(k+1))));
            end
            %%% Handling the one-dimensional case
            if numStepsY == 0
                densityM(i,1,j) = densityM(i,1,j) + (sum((cellsM{randii}{timeOfInterest}(:,1) > stepBnd(i)).*...
                    (cellsM{randii}{timeOfInterest}(:,1) <= stepBnd(i+1))));
            end
        end
    end     % loop across all simulations
end
if yBnd > 0
    densityM = densityM/(step*step);
else
    densityM = densityM/step;
end

%%% Taking the mean across all simulations
densityM = densityM/nrand;

%%% Plotting the mean cell density in time
figure
if yBnd > 0
    heatmap(1:timeStep:time,stepBnd(1:end-1),squeeze(densityM(:,floor(numStepsY/2),:)));
else
    heatmap(1:timeStep:time,stepBnd(1:end-1),squeeze(densityM(:,1,:)));
end
colorbar
xlabel('time step (days)')
ylabel('left bin edges (mm)')
set(gca,'FontSize', 12, 'FontName', 'Helvetica')
title(['Mean #cells/mm at x = 0; ball radius = ' num2str(omegaLoc) ' mm and ' num2str(N) ' positions evaluated per day']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Saving results if saveon = 1
if saveon == 1
    save([path '/' name '.mat'], 'cellsM', 'cellsX','time', 'dt', 'numMel', 'alpha', 'nrand', 'xBnd', 'yBnd',...
        'cellsX','numMelBorn','numXanBorn','initFlag','numMelDead','numXanDead','numXan','beta','eta','phi','psi','kappa','xi','N',...
        'popSwitch','deathSwitch','longSwitch','randSwitch','-v7.3');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end     % end of main function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Helper function to set initial conditions: nrand specifies the number
%%% of random simulations, time is the number of time steps, and initFlag
%%% is a flag for different initial conditions. Written by W. Duncan
%%% Martinson, Markus Schmidtchen, and Alexandria Volkening.
function [cellsM,cellsX,numMel,numXan] = initialCondition(nrand, time, initFlag, xBnd, yBnd)
cellsM = cell(nrand,time+1);
cellsX = cell(nrand,time+1);

%%% Defining our estimates of natural cell density
naturalMdensity = 400;      % in melanophores/mm^2
naturalXdensity = 772;      % in xanthophores/mm^2

%%% Also definining cell numbers
numMel = zeros(time+1,nrand);        
numXan = zeros(time+1,nrand); 

%%% Note: We place one xanthophore at [-40000,0] when there are only
%%% melanophores, and one melanophore at [-40000,0] when there are only
%%% xanthophores, just to make things work smoothly without empty matrices

for randii = 1:nrand    
    if strcmp(initFlag, 'singleMel') == 1
        %%% Single melanophore at the center of the domain, no xanthophores
        cellsM{randii}{1} = [0,0];
        cellsX{randii}{1} =[-40000,0];      % outside of the domain of interest
        numMel(1,randii) = 1;
    elseif strcmp(initFlag, 'singleXan') == 1
        %%% Single xanthophore at the center of the domain, no melanophores
        cellsM{randii}{1} = [-40000,0];     % outside of the domain of interest
        cellsX{randii}{1} =[0,0];
        numXan(1,randii) = 1;
    elseif strcmp(initFlag, 'twoMel') == 1
        %%% One melanophore at [0.6,0] and one at [-0.8,0]
        cellsM{randii}{1} = [0.6,0;-0.8,0];
        cellsX{randii}{1} =[-40000,0];
        numMel(1,randii) = 2;
    elseif strcmp(initFlag, 'twoXan') == 1
        %%% One xanthophore at [0.6,0] and one at [-0.8,0]
        cellsM{randii}{1} = [-40000,0];
        cellsX{randii}{1} =[0.6,0;-0.8,0];
        numXan(1,randii) = 2;
    elseif strcmp(initFlag, 'melXanPair') == 1
        %%% A melanohpore and an xanthophore at the center of the domain
        cellsM{randii}{1} = [0.025,0];
        cellsX{randii}{1} =[-0.025,0];
        numXan(1,randii) = 1;
        numMel(1,randii) = 1;
    elseif strcmp(initFlag, 'boxMel') == 1
        %%% A 1 mm^2 box of melanophores at their natural density at the
        %%% center of the domain
        cellsM{randii}{1} = [-500/1000 + rand(naturalMdensity,1), -500/1000 + rand(naturalMdensity,1)];
        cellsX{randii}{1} = [-40000,0];
        numMel(1,randii) = naturalMdensity;
    elseif strcmp(initFlag, 'boxXan') == 1
        %%% A 1 mm^2 box of xanthophores at their natural density at the
        %%% center of the domain
        cellsM{randii}{1} = [-40000,0];
        cellsX{randii}{1} = [-500/1000 + rand(naturalXdensity,1), -500/1000 + rand(naturalXdensity,1)];
        numXan(1,randii) = naturalXdensity;
    elseif strcmp(initFlag, 'stripeMel') == 1
        N_M = round(naturalMdensity*3*0.5);
        cellsM{randii}{1} = [-xIC/2 + xIC*rand(N_M,1), -250/1000 + 0.5*rand(N_M,1)];
        cellsX{randii}{1} = [-40000,0];
    elseif strcmp(initFlag, 'stripeXan') == 1
        N_X = round(naturalXdensity*3*0.5);
        cellsM{randii}{1} = [-40000,0]; 
        cellsX{randii}{1} = [-xIC/2 + xIC*rand(N_X,1), -250/1000 + 0.5*rand(N_X,1)];
    elseif strcmp(initFlag, 'twoBoxMel') == 1
        %%% Two boxes of melanophores at their natural density
        N_M1 = round(naturalMdensity*1*0.5);
        N_M2 = round(naturalMdensity*1*0.5);
        cellsM{randii}{1} = [-0.6 + 0.5*rand(N_M1,1), -0.75 + 1.0*rand(N_M1,1); ...
            0.1 + 0.5*rand(N_M2,1), -0.25 + 1.0*rand(N_M2,1)];
        cellsX{randii}{1} = [-40000,0];
        numMel(1,randii) = N_M1+N_M2;
    elseif strcmp(initFlag, 'twoStripeMel') == 1
        %%% Two stripes of melanophores at their natural density
        N_M1 = round(naturalMdensity*0.5*3);
        N_M2 = round(naturalMdensity*0.5*3);
        cellsM{randii}{1} = [-xIC/2 + xIC*rand(N_M1,1), 0.5*rand(N_M1,1); ...
            -xIC/2 + xIC*rand(N_M2,1), -.5820 + 0.5*rand(N_M2,1)];
        cellsX{randii}{1} = [-40000,0];
        numMel(1,randii) = N_M1+N_M2;
    elseif strcmp(initFlag, 'twoBoxXan') == 1
        %%% Two offset rectangles of xanthophores at their natural density
        N_X1 = round(naturalXdensity*1*0.5);
        N_X2 = round(naturalXdensity*1*0.5);
        cellsX{randii}{1} = [-0.6 + 0.5*rand(N_X1,1), -0.75 + 1.0*rand(N_X1,1); ...
            0.1 + 0.5*rand(N_X2,1), -0.25 + 1.0*rand(N_X2,1)];
        cellsM{randii}{1} = [-40000,0];
        numXan(1,randii) = N_X1+N_X2;
    elseif strcmp(initFlag, 'twoStripeXan') == 1
        %%% Two stripes of xanthophores at their natural density
        N_X1 = round(naturalXdensity*0.36*3);
        N_X2 = round(naturalXdensity*0.36*3);
        cellsM{randii}{1} = [-40000,0];
        cellsX{randii}{1} = [-xIC/2 + xIC*rand(N_X1,1), 0.36*rand(N_X1,1); ...
            -xIC/2 + xIC*rand(N_X2,1), -.4420 + 0.36*rand(N_X2,1)];
        numXan(1,randii) = N_X1+N_X2;
    end
end

if strcmp(initFlag, 'rand') == 1
    %%% 10 melanophores and xanthophores placed uniformly at random (same
    %%% random initial condition for all simulations
    randLocations = [-xBnd/2+xBnd*rand(20,1),-yBnd/2 + yBnd*rand(20,1)];
    for randii = 1:nrand
        cellsM{randii}{1} = randLocations(1:10,:);
        cellsX{randii}{1} = randLocations(11:20,:);
        numXan(1,randii) = 10;
        numMel(1,randii) = 10;
    end
end
end     % end of helper function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function to implement forward Euler
function [updatedM, updatedX] = EulerForward(M0,X0,NM,NX,dt,Rmm,Rxx,Rxm,Rmx,Amx,rmm,rxx,rxm,rmx,amx,bnd_X,bnd_Y,cutoff,moveSwitch)

%%% One Forward Euler step
[dMdt,dXdt] = XanMel_migration(M0, X0,NM,NX,Rmm,Rxx,Rxm,Rmx,Amx,rmm,rxx,rxm,rmx,amx,bnd_X,bnd_Y,cutoff);
updatedM = M0 + dt*moveSwitch*dMdt;
updatedX = X0 + dt*moveSwitch*dXdt;

%%% Periodic boundary conditions
updatedM(:,1) = mod(bnd_X/2+updatedM(:,1),bnd_X)-bnd_X/2;
updatedX(:,1) = mod(bnd_X/2+updatedX(:,1),bnd_X)-bnd_X/2;
updatedM(:,2) = mod(bnd_Y/2+updatedM(:,2),bnd_Y)-bnd_Y/2;
updatedX(:,2) = mod(bnd_Y/2+updatedX(:,2),bnd_Y)-bnd_Y/2;
end
% end of helper function for forward Euler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Helper function for implementing righthand side of our differential
%%% equations for cell migration, called by forward Euler routine
function [rhsM, rhsX] = XanMel_migration(M, X,NM,NX,Rmm,Rxx,Rxm,Rmx,Amx,rmm,rxx,rxm,rmx,amx,lenX,lenY,cutoff)

%%% Initializations
rhsM = zeros(NM,2);
rhsX = zeros(NX,2);

%%% Defining norms for periodic boundary conditions
myNorm = @(x,y)(sqrt(min((x(:,2)-y(:,2)).^2, (lenY-abs(x(:,2)-y(:,2))).^2) + ...
    min((x(:,1)-y(:,1)).^2, (lenX-abs(x(:,1)-y(:,1))).^2)));
mysignI = @(x,y)((abs(x-y)<=lenX/2).*(x-y)-(abs(x-y)>lenX/2).*(lenX-abs(x-y)).*sign(x-y));
mysignJ = @(x,y)((abs(x-y)<=lenY/2).*(x-y)-(abs(x-y)>lenY/2).*(lenY-abs(x-y)).*sign(x-y));

%%% These matrices hold the pairwise distances between our cells
distanceMM = pdist2(M,M,myNorm);    % distances between mel cells
distanceMX = pdist2(X,M, myNorm);   % distances between xan and mel (from xan perspective)
distanceXM = pdist2(M,X, myNorm);   % distances between xan and mel (from mel perspective)
distanceXX = pdist2(X, X, myNorm);  % distances between xan cells

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
end
% end of helper function for cell migration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Series of four helper functions called in XanMel_migration to implement
%%% individual forces between cells
%%% Short-range repulsion and long-range attraction of melanophores
function forcesMM = q_mm(d,R,r,cutoff)
forcesMM = -(R/r)*exp(-d./r).*(d < cutoff);
end

%%% Short-range repulsion of melanophores by xanthophores
function forcesXM = q_xm(d,R,r,cutoff)
forcesXM = -(R/r)*exp(-d./r).*(d < cutoff);
end

%%% Contact attraction and short-range repulsion of xanthophores by 
%%%% melanophores
function forcesMX = q_mx(d,R,r,A,a,cutoff)
forcesMX = -(R/r)*exp(-d./r).*(d < cutoff) + (A/a)*exp(-d./a).*(d < cutoff);
end

%%% Short-range repulsion of xanthophores
function forcesXX = q_xx(d,R,r,cutoff)
forcesXX = -(R/r)*exp(-d./r).*(d < cutoff);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



