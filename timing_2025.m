%% Timing Analysis -- Updated 7/17/2025


%close all; clear all;

path = ['/Volumes/KO_Portable/MEC_LEC_CA3/data_converted/']; cd(path);
path2 = ['/Volumes/KO_Portable/MEC_LEC_CA3/analysis/CC/timing/'];

indv = 0; %Do you want to recollect individual cell data from data converted (0/1)

%Find list of folders in data converted (parent) - all experimental animals
parent = dir(pwd); parent = parent([parent(:).isdir]); 
parent = parent(~ismember({parent(:).name},{'.','..'}));

if indv
    for i = 1:length(parent)
        cd(parent(i).name) %go to indv animal folder
        child = dir(pwd); child = child([child(:).isdir]); 
        child = child(~ismember({child(:).name},{'.','..'})); %get list of cell folders
           
        for ii = 1:length(child)
            cd(child(ii).name); %navigate to an individual cell
            if isdir('timing') %do a bunch of stuff if there is timing data analyzed in the correct folder labelled "timing"
                timingStruct = []; id = child(ii).name; timingStruct.id = string(id);
                if isfile('pathway.txt')
                   pathway = readcell('pathway'); 
                else
                   fprintf('pathway txt file does not exist'); pwd
                   fileID = fopen('pathway', 'w'); ptwy_temp = inputdlg({'blue', 'red'});
                   ptwy_temp = [{['blue_',ptwy_temp{1}]}; {['red_', ptwy_temp{2}]}];
                   fprintf(fileID,  '%s\n', ptwy_temp{:});
                   pathway = readcell('pathway')
                end 
                
                cd('timing'); %navigate to the correct folder
                files = dir(['PSP*csv']); %load amplitude file
                timing = readtable(files.name); %make file a table and save as variable
                timingStruct.dataKey = timing.Properties.VariableNames(2:end);
                timingStruct.data = table2array(timing(:, 2:end));
    
                idx = find(contains(timingStruct.dataKey, 'spike')==1); %change NaNs in spike prob to 0s
                for iii = 1:length(idx)
                    timingStruct.data(isnan(timingStruct.data(:,idx(iii))), idx(iii)) = 0;
                end
                
                idx = find(contains(timingStruct.dataKey, 'amp')==1); %Nan any values that are less than 0 (negative values) 
                for iii = 1:length(idx)
                    timingStruct.data(find(timingStruct.data(:,idx(iii)) <=0), idx(iii)) = NaN;
                end
    
                idx = find(contains(timingStruct.dataKey, 'sum')==1);  %Nan any values that are less than 0 (negative values) 
                for iii = 1:length(idx)
                    timingStruct.data(find(timingStruct.data(:,idx(iii)) <=0), idx(iii)) = NaN;
                end
                
                %Find mean of the sweeps
                temp = [string(cellfun(@(x) extractAfter(x, 'timing'), table2cell(timing(:,1)), 'UniformOutput', false)) ...
                    +  string(cellfun(@(x) extractBetween(x, '_1_', '_timing'),table2cell(timing(:,1))))];
                [u,~,g] = unique(temp,'stable');
                timingStruct.dataAvg = splitapply(@(x) mean(x,1, 'omitnan'), timingStruct.data, g);
                timingStruct.stimAvg = u;
                timingStruct.stim = temp;
                timingStruct.blue = pathway{string(pathway)=='blue',2}; timingStruct.red = pathway{string(pathway)=='red',2}; %store pathway information
                timingStruct.path = strcat(strcat(timingStruct.red, 'r'),strcat(timingStruct.blue, 'b'), temp);
                timingStruct.pathAvg = strcat(strcat(timingStruct.red, 'r'),strcat(timingStruct.blue, 'b'), u);
                timingStruct.pathAvgId = strcat(timingStruct.id, 'x', timingStruct.pathAvg);
                
                temp = [find(strcmp(extractAfter(timingStruct.pathAvg, 'v'), '1')== 1);find(strcmp(extractAfter(timingStruct.pathAvg, 'p'), '1') ==1); ...
                    find(strcmp(extractAfter(timingStruct.pathAvg, 'r1'), '01')== 1)];
                temp = sort(temp);
                timingStruct.ratioAvg = [];
                timingStruct.diffAvg = [];
                % calculate ratio within a cell then go below and change
                % averaging
                for i = 1:length(temp)
                    idx = temp(i);
                    temp_m = mean(timingStruct.dataAvg(idx:(idx+3), :), 'omitnan');
                    temp_ratio = timingStruct.dataAvg((i:(i+10)),:)./temp_m;
                    temp_diff = timingStruct.dataAvg((i:(i+10)),:) - temp_m;
                    timingStruct.ratioAvg = [timingStruct.ratioAvg; temp_ratio];
                    timingStruct.diffAvg = [timingStruct.diffAvg; temp_diff];
                end
                   
                cd ..
                save(['timingStruct_', id], 'timingStruct'); path = pwd;
                cd(path2)
                save(['timingStruct_', id], 'timingStruct')
                cd(path);
    %             if sum(contains(u, 'b1r10')) > 1
    %                 fprintf(id)
    %             end
            end
            cd ../
        end
       cd ../
    end
end

cd(path2); %Go to timing analysis folder

files = dir('timingStruct*.mat'); %find all individual timing files

timeAll.data = [];
timeAll.dataAvg = [];
timeAll.dataKey = [];
timeAll.blue = [];
timeAll.red = [];
timeAll.stim = [];
timeAll.stimAvg = [];
timeAll.path = [];
timeAll.pathAvg = [];
timeAll.id = [];
timeAll.ratioAvg = [];
timeAll.diffAvg = [];


for i = 1:length(files)
    load(files(i).name);
    timeAll.data = [timeAll.data; timingStruct.data];
    timeAll.dataAvg = [timeAll.dataAvg; timingStruct.dataAvg];
    timeAll.dataKey = [timeAll.dataKey; timingStruct.dataKey];
    timeAll.blue = [timeAll.blue; string(timingStruct.blue), length(timingStruct.stim), length(timingStruct.stimAvg)];
    timeAll.red = [timeAll.red; string(timingStruct.red), length(timingStruct.stim), length(timingStruct.stimAvg)];
    timeAll.stim = [timeAll.stim; timingStruct.stim];
    timeAll.stimAvg = [timeAll.stimAvg; timingStruct.stimAvg];
    timeAll.id = [timeAll.id; repelem({timingStruct.id}, length(timingStruct.pathAvg),1)]; %change this to keep track of cell numbers
    timeAll.path =[timeAll.path; timingStruct.path];
    timeAll.pathAvg =[timeAll.pathAvg; timingStruct.pathAvg];   
    timeAll.ratioAvg =[timeAll.ratioAvg; timingStruct.ratioAvg];
    timeAll.diffAvg = [timeAll.diffAvg; timingStruct.diffAvg];
end

%Combine all the low (2v and <) and high intensity (5v and >), combine all
%the no hold and 50 (depolarizing cond) 

timeAll.pathAvg(:,2) = strrep(timeAll.pathAvg, '2v', '1v'); 
timeAll.pathAvg(:,2) = strrep(timeAll.pathAvg(:,2), '4p', '1v');
timeAll.pathAvg(:,2) = strrep(timeAll.pathAvg(:,2), '7p', '1v');
timeAll.pathAvg(:,2) = strrep(timeAll.pathAvg(:,2), '5v', '10v');
timeAll.pathAvg(:,2) = strrep(timeAll.pathAvg(:,2), 'b1r10', '1v');
timeAll.pathAvg(:,2) = strrep(timeAll.pathAvg(:,2), '_nohold_', '50_');
timeAll.pathAvg(:,2) = strrep(timeAll.pathAvg(:,2), 'nohold', '50');

timeAll.pathAvgId = strcat(string(timeAll.id), 'x', timeAll.pathAvg(:,2));

%find averages for high and low intensity by cell
[u,~,g] = unique(timeAll.pathAvgId,'stable');
timeAll.dataCell = splitapply(@(x) mean(x, 1,'omitnan'), timeAll.dataAvg, g);
timeAll.ratioCell = splitapply(@(x) mean(x, 1,'omitnan'), timeAll.ratioAvg, g);
timeAll.diffCell = splitapply(@(x) mean(x, 1,'omitnan'), timeAll.diffAvg, g);
timeAll.pathCell = u;

%Find averages across pathway (mecrlecb)
timeAll.pathCell(:,2) = extractAfter(timeAll.pathCell, 'x');
[u,~,g] = unique(timeAll.pathCell(:,2),'stable');
timeAll.dataPool = splitapply(@(x) mean(x, 1,'omitnan'), timeAll.dataCell, g);
timeAll.ratioPool = splitapply(@(x) mean(x, 1,'omitnan'), timeAll.ratioCell, g);
timeAll.diffPool = splitapply(@(x) mean(x, 1,'omitnan'), timeAll.diffCell, g);
timeAll.pathPool = u;
 
%Combine all the mec data (opsin doesn't matter here)
idxEc = find(contains(timeAll.pathCell(:,2), 'mecrlecb') == 1);
tempPath = timeAll.pathCell(idxEc,:);
tempData = timeAll.dataCell(idxEc,:);
tempRatio = timeAll.ratioCell(idxEc,:);
tempDiff = timeAll.diffCell(idxEc,:);
temp2 = extractAfter(tempPath(:,2), 'v');
idxR = [3:13];
idxB = [14:24];

%Store Sweeps 1:6 Data (mecrlecb)
idxSweep = logical(sum(temp2 == string(1:6),2)); %looking for sweeps 1:6 in mecrlecb data
timeAll.mec = tempData(idxSweep,idxR); %all the 1:6 sweeps for the red led (2nd response)
timeAll.mecRatio = tempRatio(idxSweep,idxR);
timeAll.mecDiff = tempDiff(idxSweep,idxR);
timeAll.mecPath = tempPath(idxSweep,:);
timeAll.lecFirst = tempData(idxSweep,idxB); %Store blue data under first response for lec
timeAll.lecRatioFirst = tempRatio(idxSweep,idxB);
timeAll.lecDiffFirst = tempDiff(idxSweep,idxB);
timeAll.lecPathFirst = tempPath(idxSweep,:)
 
%Store Sweeps 6:11 Data (mecrlecb)
idxSweep = logical(sum(temp2 == string(6:11),2));
timeAll.lec = tempData(idxSweep,idxB); %store 6:11 sweeps for blue (fixed) led, 2nd response
timeAll.lecRatio = tempRatio(idxSweep,idxB);
timeAll.lecDiff = tempDiff(idxSweep,idxB);
timeAll.lecPath = tempPath(idxSweep,:);
timeAll.mecFirst = tempData(idxSweep,idxR); %red respones before blue responses for sweeps 6:11
timeAll.mecRatioFirst = tempRatio(idxSweep,idxR);
timeAll.mecDiffFirst = tempDiff(idxSweep,idxR);
timeAll.mecPathFirst = tempPath(idxSweep,:);

%Pull new data 
idxEc = find(contains(timeAll.pathCell(:,2), 'lecrmecb') == 1);
tempPath = timeAll.pathCell(idxEc,:);
tempData = timeAll.dataCell(idxEc,:);
tempRatio = timeAll.ratioCell(idxEc,:);
tempDiff = timeAll.diffCell(idxEc,:);
temp2 = extractAfter(tempPath(:,2), 'v');

%Store Sweeps 6:11 data (lecrmecb)
idxSweep = logical(sum(temp2 == string(6:11),2));
timeAll.mec =[timeAll.mec; tempData(idxSweep,idxB)]; %2nd response for sweeps 6:11 blue led
timeAll.mecRatio =[timeAll.mecRatio; tempRatio(idxSweep,idxB)];
timeAll.mecDiff =[timeAll.mecDiff; tempDiff(idxSweep,idxB)];
timeAll.mecPath = [timeAll.mecPath; tempPath(idxSweep,:)];
timeAll.lecFirst =[timeAll.lecFirst; tempData(idxSweep,idxR)]; %first response for sweeps 6:11 red led
timeAll.lecRatioFirst =[timeAll.lecRatioFirst; tempRatio(idxSweep,idxR)];
timeAll.lecDiffFirst =[timeAll.lecDiffFirst; tempDiff(idxSweep,idxR)];
timeAll.lecPathFirst = [timeAll.lecPathFirst; tempPath(idxSweep,:)];


%Store Sweeps 1:6 Data (lecrmecb)
idxSweep = logical(sum(temp2 == string(1:6),2));
timeAll.lec =[timeAll.lec; tempData(idxSweep,idxR)]; %2nd response for sweeps 1:6 red led 
timeAll.lecRatio =[timeAll.lecRatio; tempRatio(idxSweep,idxR)];
timeAll.lecDiff =[timeAll.lecDiff; tempDiff(idxSweep,idxR)];
timeAll.lecPath = [timeAll.lecPath; tempPath(idxSweep,:)]; 
timeAll.mecFirst =[timeAll.mecFirst; tempData(idxSweep,idxB)]; %First response for sweeps 1:6, blue led
timeAll.mecRatioFirst =[timeAll.mecRatioFirst; tempRatio(idxSweep,idxB)];
timeAll.mecDiffFirst =[timeAll.mecDiffFirst; tempDiff(idxSweep,idxB)];
timeAll.mecPathFirst = [timeAll.mecPathFirst; tempPath(idxSweep,:)];


%Change all sweeps 7:10 to 5:1
timeAvg.lecPath = replace(timeAll.lecPath, 'v11', 'v1');
timeAvg.lecPath = replace(timeAvg.lecPath, 'v10', 'v2');
timeAvg.lecPath = replace(timeAvg.lecPath, 'v9', 'v3');
timeAvg.lecPath = replace(timeAvg.lecPath, 'v8', 'v4');
timeAvg.lecPath = replace(timeAvg.lecPath, 'v7', 'v5');

timeAvg.mecPath = replace(timeAll.mecPath, 'v11', 'v1');
timeAvg.mecPath = replace(timeAvg.mecPath, 'v10', 'v2');
timeAvg.mecPath = replace(timeAvg.mecPath, 'v9', 'v3');
timeAvg.mecPath = replace(timeAvg.mecPath, 'v8', 'v4');
timeAvg.mecPath = replace(timeAvg.mecPath, 'v7', 'v5');

%Change all sweeps 7:10 to 5:1
timeAvg.lecPathFirst = replace(timeAll.lecPathFirst, 'v11', 'v1');
timeAvg.lecPathFirst = replace(timeAvg.lecPathFirst, 'v10', 'v2');
timeAvg.lecPathFirst = replace(timeAvg.lecPathFirst, 'v9', 'v3');
timeAvg.lecPathFirst = replace(timeAvg.lecPathFirst, 'v8', 'v4');
timeAvg.lecPathFirst = replace(timeAvg.lecPathFirst, 'v7', 'v5');

timeAvg.mecPathFirst = replace(timeAll.mecPathFirst, 'v11', 'v1');
timeAvg.mecPathFirst = replace(timeAvg.mecPathFirst, 'v10', 'v2');
timeAvg.mecPathFirst = replace(timeAvg.mecPathFirst, 'v9', 'v3');
timeAvg.mecPathFirst = replace(timeAvg.mecPathFirst, 'v8', 'v4');
timeAvg.mecPathFirst = replace(timeAvg.mecPathFirst, 'v7', 'v5');


timeAvg.mecPath(:,3) = extractAfter(timeAvg.mecPath(:,2), 'b');
timeAvg.lecPath(:,3) = extractAfter(timeAvg.lecPath(:,2), 'b');
timeAvg.mecPathFirst(:,3) = extractAfter(timeAvg.mecPathFirst(:,2), 'b');
timeAvg.lecPathFirst(:,3) = extractAfter(timeAvg.lecPathFirst(:,2), 'b');

[u,~,g] = unique(timeAvg.lecPath(:,3),'stable');
timeAvg.dataLec = splitapply(@(x) mean(x, 1,'omitnan'), timeAll.lec, g);
timeAvg.ratioLec = splitapply(@(x) mean(x, 1,'omitnan'), timeAll.lecRatio, g);
timeAvg.diffLec = splitapply(@(x) mean(x, 1,'omitnan'), timeAll.lecDiff, g);
timeAvg.stimLec = u;

[u,~,g] = unique(timeAvg.lecPathFirst(:,3),'stable');
timeAvg.dataLecFirst = splitapply(@(x) mean(x, 1,'omitnan'), timeAll.lecFirst, g);
timeAvg.ratioLecFirst = splitapply(@(x) mean(x, 1,'omitnan'), timeAll.lecRatioFirst, g);
timeAvg.diffLecFirst = splitapply(@(x) mean(x, 1,'omitnan'), timeAll.lecDiffFirst, g);
timeAvg.stimLecFirst = u;


[u,~,g] = unique(timeAvg.mecPath(:,3),'stable');
timeAvg.dataMec = splitapply(@(x) mean(x, 1,'omitnan'), timeAll.mec, g);
timeAvg.ratioMec = splitapply(@(x) mean(x, 1,'omitnan'), timeAll.mecRatio, g);
timeAvg.diffMec = splitapply(@(x) mean(x, 1,'omitnan'), timeAll.mecDiff, g);
timeAvg.stimMec = u;

[u,~,g] = unique(timeAvg.mecPathFirst(:,3),'stable');
timeAvg.dataMecFirst = splitapply(@(x) mean(x, 1,'omitnan'), timeAll.mecFirst, g);
timeAvg.ratioMecFirst = splitapply(@(x) mean(x, 1,'omitnan'), timeAll.mecRatioFirst, g);
timeAvg.diffMecFirst = splitapply(@(x) mean(x, 1,'omitnan'), timeAll.mecDiffFirst, g);
timeAvg.stimMecFirst = u;

timeAvg.dataKey = timeAll.dataKey(1,idxR);


%% Hard Coded Plots for Combined All MEC and LEC
%Low Intensity at -70, plot summation, amplitude, peak, spike probability 
figure; subplot(1,4,1); %summation
plot(6:11, timeAvg.dataLec(1:6, 10), 'b-', 'LineWidth',2); hold on;
plot(1:6, timeAvg.dataLecFirst(1:6,10), 'b-', 'LineWidth',2)
plot(6:11,fliplr(timeAvg.dataMec(1:6, 10)'), 'r-', 'LineWidth',2);
plot(1:6, fliplr(timeAvg.dataMecFirst(1:6, 10)'), 'r-', 'LineWidth',2);

ylabel('Summation (mV)'); xlabel('Time between Input Stim (ms)');
legend('LEC', '','MEC');
title('Low Intensity Timing Experiment at -70');
xticks([1:1:11]);
xticklabels({'-50', '-40', '-30', '-20', '-10', '0', ...
    '10','20','30','40','50'})
set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
ylim([1, 5]); xlim([1,11]);

subplot(1,4,2); %amplitude
plot(6:11, timeAvg.dataLec(1:6, 8), 'b-', 'LineWidth',2); hold on;
plot(1:6, timeAvg.dataLecFirst(1:6,8), 'b-', 'LineWidth',2)
plot(6:11,fliplr(timeAvg.dataMec(1:6, 8)'), 'r-', 'LineWidth',2);
plot(1:6, fliplr(timeAvg.dataMecFirst(1:6, 8)'), 'r-', 'LineWidth',2);

ylabel('Amplitude (mV)'); xlabel('Time between Input Stim (ms)');
legend('LEC', '','MEC');
xticks([1:1:11]);
xticklabels({'-50', '-40', '-30', '-20', '-10', '0', ...
    '10','20','30','40','50'})
set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
ylim([0, 5]); xlim([1,11]);

subplot(1,4,3); %peak
plot(6:11, timeAvg.dataLec(1:6, 6), 'b-', 'LineWidth',2); hold on;
plot(1:6, timeAvg.dataLecFirst(1:6,6), 'b-', 'LineWidth',2)
plot(6:11,fliplr(timeAvg.dataMec(1:6, 6)'), 'r-', 'LineWidth',2);
plot(1:6, fliplr(timeAvg.dataMecFirst(1:6, 6)'), 'r-', 'LineWidth',2);

ylabel('Peak (mV)'); xlabel('Time between Input Stim (ms)');
legend('LEC', '','MEC');
xticks([1:1:11]);
xticklabels({'-50', '-40', '-30', '-20', '-10', '0', ...
    '10','20','30','40','50'})
set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
xlim([1,11]);

subplot(1,4,4); %spike probability
plot(6:11, timeAvg.dataLec(1:6, 11), 'b-', 'LineWidth',2); hold on;
plot(1:6, timeAvg.dataLecFirst(1:6,11), 'b-', 'LineWidth',2)
plot(6:11,fliplr(timeAvg.dataMec(1:6, 11)'), 'r--', 'LineWidth',2);
plot(1:6, fliplr(timeAvg.dataMecFirst(1:6, 11)'), 'r--', 'LineWidth',2);

ylabel('Spike Probability'); xlabel('Time between Input Stim (ms)');
legend('LEC', '','MEC');
xticks([1:1:11]);
xticklabels({'-50', '-40', '-30', '-20', '-10', '0', ...
    '10','20','30','40','50'})
set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
xlim([1,11]); 

%% Ratios

figure; subplot(1,3,1); %summation
plot(1:6, timeAvg.ratioLec(1:6, 10), 'b-', 'LineWidth',2); hold on;
plot(1:6,fliplr(timeAvg.ratioMec(1:6, 10)'), 'r-', 'LineWidth',2);

ylabel('Summation Ratio 2nd/1st'); xlabel('Time between Input Stim (ms)');
legend('LEC','MEC');
title('Ratios: Low Intensity Timing Experiment at -70');
xticks([1:6]);
xticklabels({ '0','10','20','30','40','50'})
set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
xlim([1,6]);

subplot(1,3,2); %amplitude
plot(1:6, timeAvg.ratioLec(1:6, 8), 'b-', 'LineWidth',2); hold on;
plot(1:6,fliplr(timeAvg.ratioMec(1:6, 8)'), 'r-', 'LineWidth',2);

ylabel('Amplitude Ratio 2nd/1st'); xlabel('Time between Input Stim (ms)');
legend('LEC','MEC');
xticks([1:6]);
xticklabels({ '0','10','20','30','40','50'})
set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
xlim([1,6]);

subplot(1,3,3); %peak
plot(1:6, timeAvg.ratioLec(1:6, 6), 'b-', 'LineWidth',2); hold on;
plot(1:6,fliplr(timeAvg.ratioMec(1:6, 6)'), 'r-', 'LineWidth',2);

ylabel('Peak Ratio 2nd/1st'); xlabel('Time between Input Stim (ms)');
legend('LEC', 'MEC');
xticks([1:6]);
xticklabels({ '0','10','20','30','40','50'})
set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
xlim([1,6]);;

%% Difference
figure; subplot(1,3,1); %summation
plot(1:6, timeAvg.diffLec(1:6, 10), 'b-', 'LineWidth',2); hold on;
plot(1:6,fliplr(timeAvg.diffMec(1:6, 10)'), 'r-', 'LineWidth',2);

ylabel('Summation Diff 2nd-1st'); xlabel('Time between Input Stim (ms)');
legend('LEC','MEC');
title('Difference: Low Intensity Timing Experiment at -70');
xticks([1:6]);
xticklabels({ '0','10','20','30','40','50'})
set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
xlim([1,6]);

subplot(1,3,2); %amplitude
plot(1:6, timeAvg.diffLec(1:6, 8), 'b-', 'LineWidth',2); hold on;
plot(1:6,fliplr(timeAvg.diffMec(1:6, 8)'), 'r-', 'LineWidth',2);

ylabel('Amplitude Diff 2nd-1st'); xlabel('Time between Input Stim (ms)');
legend('LEC', 'MEC');
xticks([1:6]);
xticklabels({ '0','10','20','30','40','50'})
set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
xlim([1,6]);

subplot(1,3,3); %peak
plot(1:6, timeAvg.diffLec(1:6, 6), 'b-', 'LineWidth',2); hold on;
plot(1:6,fliplr(timeAvg.diffMec(1:6, 6)'), 'r-', 'LineWidth',2);

ylabel('Peak Diff 2nd-1st'); xlabel('Time between Input Stim (ms)');
legend('LEC', 'MEC');
xticks([1:6]);
xticklabels({ '0','10','20','30','40','50'})
set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
xlim([1,6]);

%%
%Low Intensity at -50, plot summation, amplitude, peak, spike probability 
figure; subplot(1,4,1); %summation
plot(6:11, timeAvg.dataLec(7:12, 10), 'b-', 'LineWidth',2); hold on;
plot(1:6, timeAvg.dataLecFirst(7:12,10), 'b-', 'LineWidth',2)
plot(6:11,fliplr(timeAvg.dataMec(7:12, 10)'), 'r-', 'LineWidth',2);
plot(1:6, fliplr(timeAvg.dataMecFirst(7:12, 10)'), 'r-', 'LineWidth',2);

ylabel('Summation (mV)'); xlabel('Time between Input Stim (ms)');
legend('LEC', '','MEC');
title('Low Intensity Timing Experiment at -50');
xticks([1:1:11]);
xticklabels({'-50', '-40', '-30', '-20', '-10', '0', ...
    '10','20','30','40','50'})
set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
ylim([1, 5]); xlim([1,11]);

subplot(1,4,2); %amplitude
plot(6:11, timeAvg.dataLec(7:12, 8), 'b-', 'LineWidth',2); hold on;
plot(1:6, timeAvg.dataLecFirst(7:12,8), 'b-', 'LineWidth',2)
plot(6:11,fliplr(timeAvg.dataMec(7:12, 8)'), 'r-', 'LineWidth',2);
plot(1:6, fliplr(timeAvg.dataMecFirst(7:12, 8)'), 'r-', 'LineWidth',2);

ylabel('Amplitude (mV)'); xlabel('Time between Input Stim (ms)');
legend('LEC', '','MEC');
xticks([1:1:11]);
xticklabels({'-50', '-40', '-30', '-20', '-10', '0', ...
    '10','20','30','40','50'})
set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
ylim([0, 5]); xlim([1,11]);

subplot(1,4,3); %peak
plot(6:11, timeAvg.dataLec(7:12, 6), 'b-', 'LineWidth',2); hold on;
plot(1:6, timeAvg.dataLecFirst(7:12,6), 'b-', 'LineWidth',2)
plot(6:11,fliplr(timeAvg.dataMec(7:12, 6)'), 'r-', 'LineWidth',2);
plot(1:6, fliplr(timeAvg.dataMecFirst(7:12, 6)'), 'r-', 'LineWidth',2);

ylabel('Peak (mV)'); xlabel('Time between Input Stim (ms)');
legend('LEC', '','MEC');
xticks([1:1:11]);
xticklabels({'-50', '-40', '-30', '-20', '-10', '0', ...
    '10','20','30','40','50'})
set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
xlim([1,11]);

subplot(1,4,4); %spike probability
plot(6:11, timeAvg.dataLec(7:12, 11), 'b-', 'LineWidth',2); hold on;
plot(1:6, timeAvg.dataLecFirst(7:12,11), 'b-', 'LineWidth',2)
plot(6:11,fliplr(timeAvg.dataMec(7:12, 11)'), 'r--', 'LineWidth',2);
plot(1:6, fliplr(timeAvg.dataMecFirst(7:12, 11)'), 'r--', 'LineWidth',2);

ylabel('Spike Probability'); xlabel('Time between Input Stim (ms)');
legend('LEC', '','MEC');
xticks([1:1:11]);
xticklabels({'-50', '-40', '-30', '-20', '-10', '0', ...
    '10','20','30','40','50'})
set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
xlim([1,11]); 

%%
%High Intensity at -70, plot summation, amplitude, peak, spike probability 
figure; subplot(1,4,1); %summation
plot(6:11, timeAvg.dataLec(13:18, 10), 'b-', 'LineWidth',2); hold on;
plot(1:6, timeAvg.dataLecFirst(13:18,10), 'b-', 'LineWidth',2)
plot(6:11,fliplr(timeAvg.dataMec(13:18, 10)'), 'r-', 'LineWidth',2);
plot(1:6, fliplr(timeAvg.dataMecFirst(13:18, 10)'), 'r-', 'LineWidth',2);

ylabel('Summation (mV)'); xlabel('Time between Input Stim (ms)');
legend('LEC', '','MEC');
title('High Intensity Timing Experiment at -70');
xticks([1:1:11]);
xticklabels({'-50', '-40', '-30', '-20', '-10', '0', ...
    '10','20','30','40','50'})
set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
 xlim([1,11]);

subplot(1,4,2); %amplitude
plot(6:11, timeAvg.dataLec(13:18, 8), 'b-', 'LineWidth',2); hold on;
plot(1:6, timeAvg.dataLecFirst(13:18,8), 'b-', 'LineWidth',2)
plot(6:11,fliplr(timeAvg.dataMec(13:18, 8)'), 'r-', 'LineWidth',2);
plot(1:6, fliplr(timeAvg.dataMecFirst(13:18, 8)'), 'r-', 'LineWidth',2);

ylabel('Amplitude (mV)'); xlabel('Time between Input Stim (ms)');
legend('LEC', '','MEC');
xticks([1:1:11]);
xticklabels({'-50', '-40', '-30', '-20', '-10', '0', ...
    '10','20','30','40','50'})
set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
 xlim([1,11]);

subplot(1,4,3); %peak
plot(6:11, timeAvg.dataLec(13:18, 6), 'b-', 'LineWidth',2); hold on;
plot(1:6, timeAvg.dataLecFirst(13:18,6), 'b-', 'LineWidth',2)
plot(6:11,fliplr(timeAvg.dataMec(13:18, 6)'), 'r-', 'LineWidth',2);
plot(1:6, fliplr(timeAvg.dataMecFirst(13:18, 6)'), 'r-', 'LineWidth',2);

ylabel('Peak (mV)'); xlabel('Time between Input Stim (ms)');
legend('LEC', '','MEC');
xticks([1:1:11]);
xticklabels({'-50', '-40', '-30', '-20', '-10', '0', ...
    '10','20','30','40','50'})
set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
xlim([1,11]);

subplot(1,4,4); %spike probability
plot(6:11, timeAvg.dataLec(13:18, 11), 'b-', 'LineWidth',2); hold on;
plot(1:6, timeAvg.dataLecFirst(13:18,11), 'b-', 'LineWidth',2)
plot(6:11,fliplr(timeAvg.dataMec(13:18, 11)'), 'r--', 'LineWidth',2);
plot(1:6, fliplr(timeAvg.dataMecFirst(13:18, 11)'), 'r--', 'LineWidth',2);

ylabel('Spike Probability'); xlabel('Time between Input Stim (ms)');
legend('LEC', '','MEC');
xticks([1:1:11]);
xticklabels({'-50', '-40', '-30', '-20', '-10', '0', ...
    '10','20','30','40','50'})
set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
xlim([1,11]); 

%%
%High Intensity at -50, plot summation, amplitude, peak, spike probability 
figure; subplot(1,4,1); %summation
plot(6:11, timeAvg.dataLec(19:24, 10), 'b-', 'LineWidth',2); hold on;
plot(1:6, timeAvg.dataLecFirst(19:24,10), 'b-', 'LineWidth',2)
plot(6:11,fliplr(timeAvg.dataMec(19:24, 10)'), 'r-', 'LineWidth',2);
plot(1:6, fliplr(timeAvg.dataMecFirst(19:24, 10)'), 'r-', 'LineWidth',2);

ylabel('Summation (mV)'); xlabel('Time between Input Stim (ms)');
legend('LEC', '','MEC');
title('High Intensity Timing Experiment at -50');
xticks([1:1:11]);
xticklabels({'-50', '-40', '-30', '-20', '-10', '0', ...
    '10','20','30','40','50'})
set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
 xlim([1,11]);

subplot(1,4,2); %amplitude
plot(6:11, timeAvg.dataLec(19:24, 8), 'b-', 'LineWidth',2); hold on;
plot(1:6, timeAvg.dataLecFirst(19:24,8), 'b-', 'LineWidth',2)
plot(6:11,fliplr(timeAvg.dataMec(19:24, 8)'), 'r-', 'LineWidth',2);
plot(1:6, fliplr(timeAvg.dataMecFirst(19:24, 8)'), 'r-', 'LineWidth',2);

ylabel('Amplitude (mV)'); xlabel('Time between Input Stim (ms)');
legend('LEC', '','MEC');
xticks([1:1:11]);
xticklabels({'-50', '-40', '-30', '-20', '-10', '0', ...
    '10','20','30','40','50'})
set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
 xlim([1,11]);

subplot(1,4,3); %peak
plot(6:11, timeAvg.dataLec(19:24, 6), 'b-', 'LineWidth',2); hold on;
plot(1:6, timeAvg.dataLecFirst(19:24,6), 'b-', 'LineWidth',2)
plot(6:11,fliplr(timeAvg.dataMec(19:24, 6)'), 'r-', 'LineWidth',2);
plot(1:6, fliplr(timeAvg.dataMecFirst(19:24, 6)'), 'r-', 'LineWidth',2);

ylabel('Peak (mV)'); xlabel('Time between Input Stim (ms)');
legend('LEC', '','MEC');
xticks([1:1:11]);
xticklabels({'-50', '-40', '-30', '-20', '-10', '0', ...
    '10','20','30','40','50'})
set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
xlim([1,11]);

subplot(1,4,4); %spike probability
plot(6:11, timeAvg.dataLec(19:24, 11), 'b-', 'LineWidth',2); hold on;
plot(1:6, timeAvg.dataLecFirst(19:24,11), 'b-', 'LineWidth',2)
plot(6:11,fliplr(timeAvg.dataMec(19:24, 11)'), 'r--', 'LineWidth',2);
plot(1:6, fliplr(timeAvg.dataMecFirst(19:24, 11)'), 'r--', 'LineWidth',2);

ylabel('Spike Probability'); xlabel('Time between Input Stim (ms)');
legend('LEC', '','MEC');
xticks([1:1:11]);
xticklabels({'-50', '-40', '-30', '-20', '-10', '0', ...
    '10','20','30','40','50'})
set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
xlim([1,11]); 