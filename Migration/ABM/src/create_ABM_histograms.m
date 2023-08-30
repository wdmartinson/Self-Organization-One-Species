function create_ABM_histograms(BirthOption, dt, Tf)

% runs = [400, 300, 200, 100];
NBin = [30, 60, 100, 120, 240];

runs = [500];
NBin = [240];

Len = 1.5;
NumIter = sum(runs);
NaturalMDensity = 400;
NaturalXDensity = 772;
if BirthOption == 1 % MelBox
    NumCells = round(NaturalMDensity);
    NumOutputs = 1;
elseif BirthOption == 2 % MelStripe
    NumCells = round(NaturalMDensity*0.5*3);
    NumOutputs = 1;
elseif BirthOption == 4 % XanBox
    NumCells = round(NaturalXDensity);
    NumOutputs = 1;
elseif BirthOption == 5 % XanStripe
    NumCells = round(NaturalXDensity*0.5*3);
    NumOutputs = 1;
elseif BirthOption == 11
    N_M1 = round(NaturalMDensity*0.5*3);
    N_M2 = round(NaturalMDensity*0.5*3);

    NumMelCells = N_M1 + N_M2;
    NumXanCells = 0;
    NumOutputs = 1;
elseif BirthOption == 12
    N_M1 = round(NaturalMDensity*0.5*1);
    N_M2 = round(NaturalMDensity*0.5*1);

    NumMelCells = N_M1 + N_M2;
    NumXanCells = 0;
    NumOutputs = 1;
elseif BirthOption == 13
    NumMelCells = 0;

    N_X1 = round(NaturalXDensity*0.36*3);
    N_X2 = round(NaturalXDensity*0.36*3);
    NumXanCells = N_X1 + N_X2;

    NumOutputs = 1;

elseif BirthOption == 14
    NumMelCells = 0;

    N_X1 = round(NaturalXDensity*0.5*1);
    N_X2 = round(NaturalXDensity*0.5*1);
    NumXanCells = N_X1 + N_X2;

    NumOutputs = 1;
else
    error('BirthOption not specified / accounted for.')
end

for i = 1:length(NBin)
    if isfile(['../data/sims/MC/ABM_Histograms/ABM_Histogram_MonteCarlo_BirthOption=',num2str(BirthOption),'_dt=',num2str(dt),'_runs=',num2str(NumIter),'_NumBins=',num2str(NBin(i)),'_T=',num2str(Tf),'.mat'])
        NBin(i) = NaN;
    end
end
NBin(isnan(NBin)) = [];

temp_Chr = cell(length(NBin), 1);
if  NumOutputs ~= 1
    temp_Chr_2 = cell(length(NBin), 1);
end

for k = 1:length(runs)
    if NumOutputs == 1
        load(['../data/sims/MC/MonteCarlo_BirthOption=',num2str(BirthOption),'_NumIter=',num2str(runs(k)),'_dt=',num2str(dt),'_T=',num2str(Tf),'_Cell_Locations.mat'], 'MC_Mel');
        MC = MC_Mel;
        for i = 1:length(NBin)
            if k == 1
                temp_Chr{i} = zeros(NBin(i), NBin(i), Tf + 1);
            end
            S  = linspace(-Len, Len, NBin(i) + 1);
            for j = 1:size(MC, 3)
                [temp, ~ , ~] = histcounts2(MC(:,1,j), MC(:,2,j), S, S);
                temp_Chr{i}(:,:,j) = temp_Chr{i}(:,:,j) + temp(:,:);
            end
            if k == length(runs)
                Chr = temp_Chr{i};
                % Compute the average number density
                Chr = Chr./NumIter./(2.*Len./NBin(i)).^2;
                save(['../data/sims/MC/ABM_Histograms/ABM_Histogram_MonteCarlo_BirthOption=',num2str(BirthOption),'_dt=',num2str(dt),'_runs=',num2str(NumIter),'_NumBins=',num2str(NBin(i)),'_T=',num2str(Tf),'.mat'], 'Chr');
            end
        end
    else % If NumOutputs == 2
        load(['../data/sims/MC/MonteCarlo_BirthOption=',num2str(BirthOption),'_NumIter=',num2str(runs(k)),'_dt=',num2str(dt),'_T=',num2str(Tf),'_Cell_Locations.mat'], 'MC_Mel', 'MC_Xan');
        for i = 1:length(NBin)
            if k == 1
               temp_Chr{i} = zeros(NBin(i), NBin(i), Tf + 1);
               temp_Chr_2{i} = zeros(NBin(i), NBin(i), Tf + 1);
            end % if k == 1
            S  = linspace(-Len, Len, NBin(i) + 1);
            for j = 1:size(MC_Mel, 3)
                [temp_Mel, ~ , ~] = histcounts2(MC_Mel(:,1,j), MC_Mel(:,2,j), S, S);
                [temp_Xan, ~ , ~] = histcounts2(MC_Xan(:,1,j), MC_Xan(:,2,j), S, S);
                temp_Chr{i}(:,:,j) = temp_Chr{i}(:,:,j) + temp_Mel(:,:);
                temp_Chr_2{i}(:,:,j) = temp_Chr_2{i}(:,:,j) + temp_Xan(:,:);
            end % for j
            if k == length(runs)
                Mel = temp_Chr{i};
                Xan = temp_Chr_2{i};
                % Compute the average number density
                Mel = Mel./NumIter./(2.*Len./NBin(i)).^2;
                Xan = Xan./NumIter./(2.*Len./NBin(i)).^2;
                save(['../data/sims/MC/ABM_Histograms/ABM_Histogram_MonteCarlo_BirthOption=',num2str(BirthOption),'_dt=',num2str(dt),'_runs=',num2str(NumIter),'_NumBins=',num2str(NBin(i)),'_T=',num2str(Tf),'.mat'], 'Mel', 'Xan');
            end % if
        end % for i = 1:length(NBin)
    end % if NumOutputs == 1
end
end
