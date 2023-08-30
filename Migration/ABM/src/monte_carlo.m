function [MC_Mel, MC_Xan] = monte_carlo(NumIter, Len, FinalTime, dt, BirthOption)

NumOutputs = FinalTime./1 + 1; % Divide FinalTime by 1 because you save in time increments of 1 day in CellMigration.m


MC_Mel = cell(NumIter, 1);
MC_Xan = cell(NumIter, 1);

parfor (i = 1:NumIter, 6)
    disp(i)
    % Generate a Realization
    [cellsM, cellsX, ~] = CellMigration(FinalTime, i, 'folderName', 'resultsName', BirthOption, dt);

    % Extract cell locations at each time step.
    if BirthOption == 1 || BirthOption == 2 || BirthOption == 11 || BirthOption == 12
        MC_Mel{i} = cellsM(:,:,:); % Extract only melanophore data if you only have melanophores in the simulation
    elseif BirthOption == 4 || BirthOption == 5 || BirthOption == 13 || BirthOption == 14
        MC_Xan{i} = cellsX(:,:,:); % Extract only xanthophore data if you only have xanthophores in the simulation
    elseif BirthOption == 3 || BirthOption == 6 || BirthOption == 7 || BirthOption == 8 || BirthOption == 9 || BirthOption == 10
        MC_Mel{i} = cellsM(:,:,:);
        MC_Xan{i} = cellsX(:,:,:);
    else
        MC_Mel{i} = cellsM(:,:,:);
        MC_Xan{i} = cellsX(:,:,:);
        fprintf('Code has not been written yet to extract both melanophore and xanthophore data for this case!');
    end
end % parfor i


MC_Mel      = cell2mat(MC_Mel);
MC_Xan      = cell2mat(MC_Xan);

% save the array
save(['../data/sims/MC/MonteCarlo_BirthOption=',num2str(BirthOption),'_NumIter=',num2str(NumIter),'_dt=',num2str(dt),'_T=',num2str(FinalTime),'_Cell_Locations.mat'], 'MC_Mel', 'MC_Xan', '-v7.3');


end % function
