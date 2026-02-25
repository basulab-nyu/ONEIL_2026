%% Processing STP Data from MASS files (not individual mice)
% Does not utilize individual files 
% Last Edit: 10/20/2025

clear all; close all; 
addpath '/Volumes/KO_Portable/MATLAB/EphysAnalysis_2025'

path = ['/Volumes/KO_Portable/MEC_LEC_CA3/analysis/VC/STP/csv']; cd(path); 
plt = 0; %Do you want to plt (0 or 1)

files = dir('*_chronos_only.csv');
files = files(~startsWith(extractfield(files, 'name'), '.'));
kinetics = ["sum", "amp","indvAuc"];

files_kinetics = files(contains(extractfield(files, 'name'), kinetics));
files_mec = files_kinetics(contains(extractfield(files_kinetics,'name'), 'mec'));
files_lec = files_kinetics(contains(extractfield(files_kinetics,'name'), 'lec'));

for i = 1:length(files_mec)
    stp_temp = readtable(files_mec(i).name);
    stp_sweep = string(table2array(stp_temp(:,1)));

    stp_temp = stp_temp(~contains(stp_sweep, 'average','IgnoreCase',true),:); 
    stp_sweep = stp_sweep(~contains(stp_sweep, 'average','IgnoreCase',true),:);
    stp_temp = stp_temp(~contains(stp_sweep, 'drug','IgnoreCase',true),:); 
    stp_sweep = stp_sweep(~contains(stp_sweep, 'drug','IgnoreCase',true),:);
    %stp_temp = stp_temp(~contains(stp_sweep, '_10v','IgnoreCase',true),:); 
    %stp_sweep = stp_sweep(~contains(stp_sweep, '_10v','IgnoreCase',true),:);

    stp_temp = table2array(stp_temp(:,2:end));

    if contains(files_mec(i).name, 'hz10') && ~contains(files_mec(i).name, 'indvAuc') 
        stp_temp(stp_temp < 0) = NaN; 
    elseif contains(files_mec(i).name, 'hz80') && ~contains(files_mec(i).name, 'indvAuc') 
        stp_temp(stp_temp > 0) = NaN; 
    elseif contains(files_mec(i).name, 'indvAuc') 
        stp_temp(stp_temp < 0) = NaN; 
    end

    temp_name = extractBefore(files_mec(i).name, '_mec');
    stpMec.(temp_name).cellsAll = stp_sweep;
    stpMec.(temp_name).dataAll = stp_temp;
    stpMec.(temp_name).pprAll = stp_temp./stp_temp(:,1);  
    stpMec.(temp_name).pprAll(stpMec.(temp_name).pprAll >25) = NaN;
    stpMec.(temp_name).pprEnd = mean(stpMec.(temp_name).pprAll(:, end-2:end),2, 'omitnan');

    stpMec.(temp_name).cellsAll(:,2) = replace(string(stp_sweep), '_7p', '_1v');
    stpMec.(temp_name).cellsAll(:,2) = replace(string(stpMec.(temp_name).cellsAll(:,2)), '_5p', '_1v');
    stpMec.(temp_name).cellsAll(:,2) = replace(string(stpMec.(temp_name).cellsAll(:,2)), '_3p', '_1v');
    stpMec.(temp_name).cellsAll(:,2) = replace(string(stpMec.(temp_name).cellsAll(:,2)), '_2p', '_1v');

    [u, ~, g] = unique(extractAfter(stpMec.(temp_name).cellsAll(:,2), '20p'), 'stable'); %Find unique cell by intensity values
    stpMec.(temp_name).dataCell = splitapply (@(x) mean(x, 1, 'omitnan'), stpMec.(temp_name).dataAll, g); %Average for each cell maintain intensity (ie average across all sweeps)
    stpMec.(temp_name).cvCell = splitapply (@(x) 100*(std(x, 0, 1, 'omitnan') ./ mean(x, 1, 'omitnan')), stpMec.(temp_name).dataAll, g); %Average for each cell maintain intensity (ie average across all sweeps)
    stpMec.(temp_name).dataCellId = u; %Cell Ids
    stpMec.(temp_name).pprCell = splitapply (@(x) mean(x, 1, 'omitnan'), stpMec.(temp_name).pprAll, g); %Average ppr for each cell maintain intensity
    stpMec.(temp_name).cvPprCell = splitapply (@(x) 100*(std(x, 0,1, 'omitnan') ./ mean(x, 1, 'omitnan')), stpMec.(temp_name).pprAll, g); %Average for each cell maintain intensity (ie average across all sweeps)
    stpMec.(temp_name).pprEndCell =  splitapply (@(x) mean(x, 1, 'omitnan'), stpMec.(temp_name).pprEnd, g);

    [u, ~, g] = unique(extractBefore(u, '202'), 'stable'); % Find Unique Intensity Values
    stpMec.(temp_name).dataAvg = splitapply (@(x) mean(x, 1, 'omitnan'), stpMec.(temp_name).dataCell, g); %Average across cells for each intensity
    stpMec.(temp_name).stdAvg = splitapply (@(x) std(x,0, 1, 'omitnan'), stpMec.(temp_name).dataCell, g);
    stpMec.(temp_name).nAvg = histc(g, unique(g,'stable'));
    stpMec.(temp_name).semAvg = stpMec.(temp_name).stdAvg./sqrt(stpMec.(temp_name).nAvg);

    stpMec.(temp_name).pprAvg = splitapply (@(x) mean(x, 1, 'omitnan'), stpMec.(temp_name).pprCell, g); %Average across cells for each intensity
    stpMec.(temp_name).pprStdAvg = splitapply (@(x) std(x,0, 1, 'omitnan'), stpMec.(temp_name).pprCell, g);
    stpMec.(temp_name).pprNAvg = histc(g, unique(g,'stable'));
    stpMec.(temp_name).pprSemAvg = stpMec.(temp_name).pprStdAvg./sqrt(stpMec.(temp_name).pprNAvg);
    
    stpMec.(temp_name).pprEndAvg = splitapply (@(x) mean(x, 1, 'omitnan'), stpMec.(temp_name).pprEndCell, g); %Average across cells for each intensity
    stpMec.(temp_name).pprEndStdAvg = splitapply (@(x) std(x,0, 1, 'omitnan'), stpMec.(temp_name).pprEndCell, g);
    stpMec.(temp_name).pprEndNAvg = histc(g, unique(g, 'stable'));
    stpMec.(temp_name).pprEndSemAvg = stpMec.(temp_name).pprEndStdAvg./sqrt(stpMec.(temp_name).pprEndNAvg);

    stpMec.(temp_name).cvAvg = splitapply (@(x) mean(x, 1, 'omitnan'), stpMec.(temp_name).cvCell, g); %Average across cells for each intensity
    stpMec.(temp_name).cvStdAvg = splitapply (@(x) std(x,0, 1, 'omitnan'), stpMec.(temp_name).cvCell, g);
    stpMec.(temp_name).cvNAvg = histc(g, unique(g, 'stable'));
    stpMec.(temp_name).cvSemAvg = stpMec.(temp_name).cvStdAvg./sqrt(stpMec.(temp_name).cvNAvg);

    stpMec.(temp_name).cvPprAvg = splitapply (@(x) mean(x, 1, 'omitnan'), stpMec.(temp_name).cvPprCell, g); %Average across cells for each intensity
    stpMec.(temp_name).cvPprStdAvg = splitapply (@(x) std(x,0, 1, 'omitnan'), stpMec.(temp_name).cvPprCell, g);
    stpMec.(temp_name).cvPprNAvg = histc(g, unique(g, 'stable'));
    stpMec.(temp_name).cvPprSemAvg = stpMec.(temp_name).cvPprStdAvg./sqrt(stpMec.(temp_name).cvPprNAvg);
    
    stpMec.(temp_name).avgId = u;
end

for i = 1:length(files_lec)
    stp_temp = readtable(files_lec(i).name);
    stp_sweep = string(table2array(stp_temp(:,1)));

    stp_temp = stp_temp(~contains(stp_sweep, 'average','IgnoreCase',true),:); 
    stp_sweep = stp_sweep(~contains(stp_sweep, 'average','IgnoreCase',true),:);
    stp_temp = stp_temp(~contains(stp_sweep, 'drug','IgnoreCase',true),:); 
    stp_sweep = stp_sweep(~contains(stp_sweep, 'drug','IgnoreCase',true),:);

    stp_temp = table2array(stp_temp(:,2:end));

    if contains(files_lec(i).name, 'hz10') && ~contains(files_lec(i).name, 'indvAuc') 
        stp_temp(stp_temp < 0) = NaN; 
    elseif contains(files_lec(i).name, 'hz80') && ~contains(files_lec(i).name, 'indvAuc') 
        stp_temp(stp_temp > 0) = NaN; 
    elseif contains(files_lec(i).name, 'indvAuc') 
        stp_temp(stp_temp < 0) = NaN; 
    end

    temp_name = extractBefore(files_lec(i).name, '_lec');
    stpLec.(temp_name).cellsAll = stp_sweep;
    stpLec.(temp_name).dataAll = stp_temp;
    stpLec.(temp_name).pprAll = stp_temp./stp_temp(:,1); 
    stpLec.(temp_name).pprAll(stpLec.(temp_name).pprAll > 25) = NaN;
    stpLec.(temp_name).pprEnd = mean(stpLec.(temp_name).pprAll(:, end-2:end),2, 'omitnan');
    
    %stpLec.(temp_name).cellsAll(:,2) = replace(string(stp_sweep), '_2p', '_7p');
    %stpLec.(temp_name).cellsAll(:,2) = replace(string(stpLec.(temp_name).cellsAll(:,2)), '_3p', '_1v');
    %stpLec.(temp_name).cellsAll(:,2) = replace(string(stpLec.(temp_name).cellsAll(:,2)), '_4p', '_1v');
    %stpLec.(temp_name).cellsAll(:,2) = replace(string(stpLec.(temp_name).cellsAll(:,2)), '_5p', '_1v');
    stpLec.(temp_name).cellsAll(:,2) = replace(string(stp_sweep), '_7p', '_1v');
    stpLec.(temp_name).cellsAll(:,2) = replace(string(stpLec.(temp_name).cellsAll(:,2)), '_5p', '_1v');
    stpLec.(temp_name).cellsAll(:,2) = replace(string(stpLec.(temp_name).cellsAll(:,2)), '_3p', '_1v');
    stpLec.(temp_name).cellsAll(:,2) = replace(string(stpLec.(temp_name).cellsAll(:,2)), '_4p', '_1v');
    stpLec.(temp_name).cellsAll(:,2) = replace(string(stpLec.(temp_name).cellsAll(:,2)), '_2p', '_1v');

    [u, ~, g] = unique(extractAfter(stpLec.(temp_name).cellsAll(:,2), '20p'), 'stable'); %Find unique cell by intensity values
    stpLec.(temp_name).dataCell = splitapply (@(x) mean(x, 1, 'omitnan'), stpLec.(temp_name).dataAll, g) %Average for each cell maintain intensity (ie average across all sweeps)
    stpLec.(temp_name).cvCell = splitapply (@(x) 100*(std(x, 0,1, 'omitnan') ./ mean(x, 1, 'omitnan')), stpLec.(temp_name).dataAll, g); %Average for each cell maintain intensity (ie average across all sweeps)
    stpLec.(temp_name).dataCellId = u; %Cell Ids
    stpLec.(temp_name).pprCell = splitapply (@(x) mean(x, 1, 'omitnan'), stpLec.(temp_name).pprAll, g); %Average ppr for each cell maintain intensity
    stpLec.(temp_name).cvPprCell = splitapply (@(x) 100*(std(x, 0,1, 'omitnan') ./ mean(x, 1, 'omitnan')), stpLec.(temp_name).pprAll, g); %Average for each cell maintain intensity (ie average across all sweeps)
    stpLec.(temp_name).pprEndCell =  splitapply (@(x) mean(x, 1, 'omitnan'), stpLec.(temp_name).pprEnd, g);

    [u, ~, g] = unique(extractBefore(u, '202'), 'stable'); % Find Unique Intensity Values
    stpLec.(temp_name).dataAvg = splitapply (@(x) mean(x, 1, 'omitnan'), stpLec.(temp_name).dataCell, g); %Average across cells for each intensity
    stpLec.(temp_name).stdAvg = splitapply (@(x) std(x,0, 1, 'omitnan'), stpLec.(temp_name).dataCell, g);
    stpLec.(temp_name).nAvg = histc(g, unique(g, 'stable'));
    stpLec.(temp_name).semAvg = stpLec.(temp_name).stdAvg./sqrt(stpLec.(temp_name).nAvg);

    stpLec.(temp_name).pprAvg = splitapply (@(x) mean(x, 1, 'omitnan'), stpLec.(temp_name).pprCell, g); %Average across cells for each intensity
    stpLec.(temp_name).pprStdAvg = splitapply (@(x) std(x,0, 1, 'omitnan'), stpLec.(temp_name).pprCell, g);
    stpLec.(temp_name).pprNAvg = histc(g, unique(g, 'stable'));
    stpLec.(temp_name).pprSemAvg = stpLec.(temp_name).pprStdAvg./sqrt(stpLec.(temp_name).pprNAvg);

    stpLec.(temp_name).pprEndAvg = splitapply (@(x) mean(x, 1, 'omitnan'), stpLec.(temp_name).pprEndCell, g); %Average across cells for each intensity
    stpLec.(temp_name).pprEndStdAvg = splitapply (@(x) std(x,0, 1, 'omitnan'), stpLec.(temp_name).pprEndCell, g);
    stpLec.(temp_name).pprEndNAvg = histc(g, unique(g, 'stable'));
    stpLec.(temp_name).pprEndSemAvg = stpLec.(temp_name).pprEndStdAvg./sqrt(stpLec.(temp_name).pprEndNAvg);

    stpLec.(temp_name).cvAvg = splitapply (@(x) mean(x, 1, 'omitnan'), stpLec.(temp_name).cvCell, g); %Average across cells for each intensity
    stpLec.(temp_name).cvStdAvg = splitapply (@(x) std(x,0, 1, 'omitnan'), stpLec.(temp_name).cvCell, g);
    stpLec.(temp_name).cvNAvg = histc(g, unique(g, 'stable'));
    stpLec.(temp_name).cvSemAvg = stpLec.(temp_name).cvStdAvg./sqrt(stpLec.(temp_name).cvNAvg);

    stpLec.(temp_name).cvPprAvg = splitapply (@(x) mean(x, 1, 'omitnan'), stpLec.(temp_name).cvPprCell, g); %Average across cells for each intensity
    stpLec.(temp_name).cvPprStdAvg = splitapply (@(x) std(x,0, 1, 'omitnan'), stpLec.(temp_name).cvPprCell, g);
    stpLec.(temp_name).cvPprNAvg = histc(g, unique(g, 'stable'));
    stpLec.(temp_name).cvPprSemAvg = stpLec.(temp_name).cvPprStdAvg./sqrt(stpLec.(temp_name).cvPprNAvg);
    
    stpLec.(temp_name).avgId = u;
end

datasets = {'stpLec', 'stpMec'};

for h = 1:length(datasets)
    s = evalin('base', datasets{h});
    f = fieldnames(s);
    for i = 1:length(f)
        % Step 1: Preprocess cell identifiers
        s.(f{i}).scellsAll(:,2) = replace(s.(f{i}).cellsAll(:,1), '_5p', '_7p');   % Replace "_5p" with "_7p"
        s.(f{i}).scellsAll(:,2) = replace(s.(f{i}).cellsAll(:,2), '_3p', '_7p');   % Replace "_3p" with "_7p"
        s.(f{i}).scellsAll(:,2) = replace(s.(f{i}).cellsAll(:,2), '_4p', '_7p');   % Replace "_2p" with "_7p"

        % Add two new columns:
        s.(f{i}).scellsAll(:,end:end+1) = [extractAfter(s.(f{i}).cellsAll(:,2), '_1_'), ... % Extract after "_1_"
                                          extractBetween(s.(f{i}).cellsAll(:,2), '_1_','202')]; % Extract between "_1_" and "202"

        s.(f{i}).ssweepAll = str2double(extractBefore(s.(f{i}).cellsAll(:,end), '_2'));    % Extract sweep numbers
        s.(f{i}).spprAll = s.(f{i}).dataAll ./ s.(f{i}).dataAll(:,1);                      % Compute PPR normalized to first column


        % Step 2: Compute averages grouped by unique cell IDs
        [u, ~, g] = unique(s.(f{i}).scellsAll(:,2), 'stable');                            % Unique cell IDs (keep order)

        s.(f{i}).sdataCellAvg = splitapply(@(x) mean(x,1,'omitnan'), s.(f{i}).dataAll, g); % Average data per cell per sweep
        s.(f{i}).scellAvgStim = u;                                                         % Store cell ID list
        s.(f{i}).spprCellAvg = splitapply(@(x) mean(x,1,'omitnan'), s.(f{i}).pprAll, g);   % Average PPR per cell


        % Step 3: Compute averages grouped by sweep (ignore "202..." suffix)
        [u, ~, g] = unique(extractBefore(u, '202'), 'stable');                            % Collapse IDs by removing "202..."

        s.(f{i}).sdataSweepAvg = splitapply(@(x) mean(x,1,'omitnan'), s.(f{i}).sdataCellAvg, g); % Sweep-averaged data
        s.(f{i}).sdataSweepStd = splitapply(@(x) std(x,0,1,'omitnan'), s.(f{i}).sdataCellAvg, g);    % Sweep-wise std
        s.(f{i}).sdataSweepN = histc(g, unique(g,'stable'));                                            % Count cells per sweep
        s.(f{i}).sdataSweepSem = s.(f{i}).sdataSweepStd ./ sqrt(s.(f{i}).sdataSweepN);           % Sweep-wise SEM


        % Step 4: Compute sweep-level statistics for PPR
        s.(f{i}).spprSweepAvg = splitapply(@(x) mean(x,1,'omitnan'), s.(f{i}).spprCellAvg, g);  % Sweep-averaged PPR
        s.(f{i}).spprSweepStd = splitapply(@(x) std(x,0,1,'omitnan'), s.(f{i}).spprCellAvg, g);   % Sweep-wise std of PPR
        s.(f{i}).spprSweepN = histc(g, unique(g, 'stable'));                                             % Count cells per sweep (PPR)
        s.(f{i}).spprSweepSem = s.(f{i}).spprSweepStd ./ sqrt(s.(f{i}).spprSweepN);              % Sweep-wise SEM of PPR

        s.(f{i}).ssweepAvgStim = u;                                                            % Store sweep identifiers

        % Step 5: Average across cells 
        [u, ~, g] = unique(extractAfter(s.(f{i}).scellsAll(:,2), 'hz'), 'stable');          %Unique ID's (keep order)

        s.(f{i}).sdataIndvAvg = splitapply(@(x) mean(x,1,'omitnan'), s.(f{i}).dataAll, g); % Average data across sweeps per cell by intensity
        s.(f{i}).sindvAvgStim = u;                                                         % Store cell ID list
        s.(f{i}).spprIndvAvg = splitapply(@(x) mean(x,1,'omitnan'), s.(f{i}).pprAll, g);   % Average PPR per cell

        % Step 6: Average across intensities
        [u, ~, g] = unique(extractBefore(u, '202'), 'stable');          %Unique ID's (keep order)

        s.(f{i}).sdataAvg = splitapply(@(x) mean(x,1,'omitnan'),  s.(f{i}).sdataIndvAvg, g); % Average data across sweeps per cell by intensity
        s.(f{i}).savgStim = u;                                                         % Store cell ID list
        s.(f{i}).spprAvg = splitapply(@(x) mean(x,1,'omitnan'), s.(f{i}).spprIndvAvg, g);   % Average PPR per cell
        s.(f{i}).spprAvgStd = splitapply(@(x) std(x,0,1,'omitnan'), s.(f{i}).spprIndvAvg, g);    % Sweep-wise std
        s.(f{i}).spprAvgN = histc(g, unique(g, 'stable'));        
        s.(f{i}).spprAvgSem = s.(f{i}).spprAvgStd ./ sqrt(s.(f{i}).spprAvgN);  
    end
    assignin('base', datasets{h}, s);
end

for i = 1:5
    ms_axis(i,:) = [0.1:0.1:(0.1*20)]+(15*(i-1));
end

% Add up IndvAuc

for h = 1:length(datasets)
    s = evalin('base', datasets{h});
    f = fieldnames(s);
    f = f(contains(fieldnames(s), 'Auc'));
    for i = 1:length(f)
        s.(f{i}).dataSum = sum(s.(f{i}).dataCell,2);
        [u, ~, g] = unique(extractBefore(s.(f{i}).dataCellId, '202'), 'stable');
        s.(f{i}).dataSumAvg = splitapply(@(x) mean(x, 'omitnan'), s.(f{i}).dataSum, g);
        s.(f{i}).dataSumStd = splitapply(@(x) std(x,'omitnan'), s.(f{i}).dataSum, g);
        s.(f{i}).dataSumN = histc(g, unique(g,'stable'));
        s.(f{i}).dataSumSem = s.(f{i}).dataSumStd ./ sqrt(s.(f{i}).dataSumN);
    end
    assignin('base', datasets{h}, s);
end




%% Plotting 
if plt

figure; %20p10hz at -80 LEC v MEC PPR
    errorbar(1:20, stpLec.amp_10hz80.pprAvg(1,:), stpLec.amp_10hz80.pprSemAvg(1,:), '.-','MarkerSize', 30, 'Color', '#289e1e','LineWidth',2); hold on;
    errorbar(1:20, stpMec.amp_10hz80.pprAvg(1,:), stpMec.amp_10hz80.pprSemAvg(1,:), 'r.-', 'MarkerSize', 30, 'LineWidth',2);
    yline(1, '--', 'LineWidth',3);
    legend('LEC(22)', 'MEC(54)', 'FontSize', 16); ylabel('PPR (Pn/P1)'); xlabel('Pulse Number', 'FontSize', 16);
    set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
    ylim([0.5 2]); title('10hz PPR at -80');

figure; %20p20hz at -80 LEC v MEC PPR
    errorbar(1:20, stpLec.amp_20hz80.pprAvg(1,:), stpLec.amp_20hz80.pprSemAvg(1,:), '.-','MarkerSize', 30, 'Color', '#289e1e','LineWidth',2); hold on;
    errorbar(1:20, stpMec.amp_20hz80.pprAvg(2,:), stpMec.amp_20hz80.pprSemAvg(2,:), 'r.-', 'MarkerSize', 30, 'LineWidth',2);
    yline(1, '--', 'LineWidth',3);
    legend('LEC(19)', 'MEC(49)', 'FontSize', 16); ylabel('PPR (Pn/P1)'); xlabel('Pulse Number', 'FontSize', 16);
    set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
    ylim([0.5 2]); title('20hz PPR at -80');

 figure; %20p10hz at -80 LEC v MEC Sum  PPR
    errorbar(1:20, stpLec.sum_10hz80.avgPpr(1,:), stpLec.sum_10hz80.semPpr(1,:), '.-','MarkerSize', 30, 'Color', '#289e1e','LineWidth',2); hold on;
    errorbar(1:20, stpMec.sum_10hz80.avgPpr(1,:), stpMec.sum_10hz80.semPpr(1,:), 'r.-', 'MarkerSize', 30, 'LineWidth',2);
    yline(1, '--', 'LineWidth',3);
    legend('LEC(23)', 'MEC(54)', 'FontSize', 16); ylabel('Sum PPR (Pn/P1)'); xlabel('Pulse Number', 'FontSize', 16);
    set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
    ylim([0.5 2.5]); title('10hz Sum PPR at -80');

figure; %20p10hz at -80 LEC v MEC Sum PPR
    errorbar(1:20, stpLec.sum_20hz80.avgPpr(1,:), stpLec.sum_20hz80.semPpr(1,:), '.-','MarkerSize', 30, 'Color', '#289e1e','LineWidth',2); hold on;
    errorbar(1:20, stpMec.sum_20hz80.avgPpr(1,:), stpMec.sum_20hz80.semPpr(1,:), 'r.-', 'MarkerSize', 30, 'LineWidth',2);
    yline(1, '--', 'LineWidth',3);
    legend('LEC(21)', 'MEC(49)', 'FontSize', 16); ylabel('Sum PPR (Pn/P1)'); xlabel('Pulse Number', 'FontSize', 16);
    set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
    ylim([0.5 2.5]); title('20hz Sum PPR at -80');  


% Plot +10
figure; %20p10hz at -80 LEC v MEC PPR
    errorbar(1:20, stpLec.amp_10hz10.avgPpr(1,:), stpLec.amp_10hz10.semPpr(1,:), '.-','MarkerSize', 30, 'Color', '#289e1e','LineWidth',2); hold on;
    errorbar(1:20, stpMec.amp_10hz10.avgPpr(1,:), stpMec.amp_10hz10.semPpr(1,:), 'r.-', 'MarkerSize', 30, 'LineWidth',2);
    yline(1, '--', 'LineWidth',3);
    legend('LEC(12)', 'MEC(31)', 'FontSize', 16); ylabel('PPR (Pn/P1)'); xlabel('Pulse Number', 'FontSize', 16);
    set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
    ylim([0.2 1.4]); title('10hz PPR at +10');

figure; %20p20hz at +10 LEC v MEC PPR
    errorbar(1:20, stpLec.amp_20hz10.avgPpr(1,:), stpLec.amp_20hz10.semPpr(1,:), '.-','MarkerSize', 30, 'Color', '#289e1e','LineWidth',2); hold on;
    errorbar(1:20, stpMec.amp_20hz10.avgPpr(1,:), stpMec.amp_20hz10.semPpr(1,:), 'r.-', 'MarkerSize', 30, 'LineWidth',2);
    yline(1, '--', 'LineWidth',3);
    legend('LEC(12)', 'MEC(32)', 'FontSize', 16); ylabel('PPR (Pn/P1)'); xlabel('Pulse Number', 'FontSize', 16);
    set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
    ylim([0 1.4]); title('20hz PPR at +10');

end

cd ../ 
save('STP_amp', 'stpMec', 'stpLec','ms_axis');