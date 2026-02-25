%% Timing Analysis -- Updated 7/17/2025
% Compile timingStructures in Analysis Folder and average

close all; clear all;

path2 = ['/Volumes/KO_Portable/MEC_LEC_CA3/analysis/CC/timing/']; cd(path2);

files = dir('timingStruct*.mat'); %find all individual timing files

load(files(1).name); 
fields = fieldnames(timingStruct);
for i = 1:length(fields);
    tField = fields{i};
    timeCompile.(tField) = [];
end

for i = 1:length(files)
    load(files(i).name); 
    for ii = 1:length(fields)
        tField = fields{ii};
        timeCompile.(tField) = [timeCompile.(tField); timingStruct.(tField)];
    end
end

fields = fieldnames(timeCompile);
fields_idx = find(contains(fields, 'pathAvg','IgnoreCase',true)==1);

for i = 1:length(fields_idx)
    k = fields_idx(i);
    timeCompile.(fields{k})(:,2) = strrep(timeCompile.(fields{k}), '2v', '1v');
    timeCompile.(fields{k})(:,2)  = strrep(timeCompile.(fields{k})(:,2) , '4p', '1v');
    timeCompile.(fields{k})(:,2)  = strrep(timeCompile.(fields{k})(:,2) , '7p', '1v');
    timeCompile.(fields{k})(:,2)  = strrep(timeCompile.(fields{k})(:,2) , '5v', '10v');
    timeCompile.(fields{k})(:,2)  = strrep(timeCompile.(fields{k})(:,2) , 'b1r10', '1v');
    timeCompile.(fields{k})(:,2)  = strrep(timeCompile.(fields{k})(:,2) , '_nohold_', '50_');
    timeCompile.(fields{k})(:,2)  = strrep(timeCompile.(fields{k})(:,2) , 'nohold', '50');
end

timeCompile.sRatioValPathAvg = extractAfter(timeCompile.sRatioValPathAvgId,'x');
%Find averages for high and low intensity by cell (Cell Averages)
[u,~,g] = unique(timeCompile.pathAvgId(:,2),'stable');
timeCompile.dataCell = splitapply(@(x) mean(x, 1,'omitnan'), timeCompile.dataAvg, g);
timeCompile.ratioCell = splitapply(@(x) mean(x, 1,'omitnan'), timeCompile.ratioAvg, g);
timeCompile.sRatioCell = splitapply(@(x) mean(x, 1,'omitnan'), timeCompile.sRatioAvg, g);
timeCompile.diffCell = splitapply(@(x) mean(x, 1,'omitnan'), timeCompile.diffAvg, g);
timeCompile.brRatioCell = [splitapply(@(x) mean(x, 1,'omitnan'), timeCompile.brRatioAvg, g), splitapply(@(x) mean(x, 1,'omitnan'), timeCompile.brRatioSumAvg, g)];
timeCompile.pathCell = u;
[u,~,g] = unique(timeCompile.sRatioValPathAvgId(:,2),'stable');
timeCompile.sRatioValCell = splitapply(@(x) mean(x, 1,'omitnan'), timeCompile.sRatioVal, g);
timeCompile.sRatioValPathCell = u;
timeCompile.sRatioValPathCell(:,2) = extractAfter(timeCompile.sRatioValPathCell, 'x');

%Find averages for high/low intensity (opsins remain separate)
timeCompile.pathCell(:,2) = extractAfter(timeCompile.pathCell, 'x');
[u,~,g] = unique(timeCompile.pathCell(:,2),'stable');
timeCompile.dataPool = splitapply(@(x) mean(x, 1,'omitnan'), timeCompile.dataCell, g);
timeCompile.ratioPool = splitapply(@(x) mean(x, 1,'omitnan'), timeCompile.ratioCell, g);
timeCompile.sRatioPool = splitapply(@(x) mean(x, 1,'omitnan'), timeCompile.sRatioCell, g);
timeCompile.diffPool = splitapply(@(x) mean(x, 1,'omitnan'), timeCompile.diffCell, g);
timeCompile.brRatioPool = splitapply(@(x) mean(x, 1,'omitnan'), timeCompile.brRatioCell, g);
timeCompile.pathPool = u;
[u,~,g] = unique(timeCompile.sRatioValPathCell(:,2),'stable');
timeCompile.sRatioValPool = splitapply(@(x) mean(x, 1,'omitnan'), timeCompile.sRatioValCell, g);
timeCompile.sRatioValpathPool = u;

%Find averages for high/low intensity for indv inputs (ie mec and lec)
%irrespective of opsin
idxEc = contains(timeCompile.pathCell(:,2), 'mecrlecb') == 1;
tempPath = timeCompile.pathCell(idxEc,:);
tempData = timeCompile.dataCell(idxEc,:);
tempRatio = timeCompile.ratioCell(idxEc,:);
tempsRatio = timeCompile.sRatioCell(idxEc,:);
%tempsRatioVal = timeCompile.sRatioValCell(idxEc,:);
tempDiff = timeCompile.diffCell(idxEc,:);
tempbrRatio = timeCompile.brRatioCell(idxEc,:);
temp = extractAfter(tempPath(:,2), 'v');

idxEc2 = find(contains(timeCompile.pathCell(:,2), 'lecrmecb') == 1);
tempPath2 = timeCompile.pathCell(idxEc2,:);
tempData2 = timeCompile.dataCell(idxEc2,:);
tempRatio2 = timeCompile.ratioCell(idxEc2,:);
tempsRatio2 = timeCompile.sRatioCell(idxEc2,:);
%tempsRatioVal2 = timeCompile.sRatioValCell(idxEc2,:);
tempDiff2 = timeCompile.diffCell(idxEc2,:);
tempbrRatio2 = timeCompile.brRatioCell(idxEc2,:);
temp2 = extractAfter(tempPath2(:,2), 'v');

idxR = [3:13];
idxB = [14:24];

%Store Sweeps for indv inputs
idxSweep = logical(sum(temp == string(1:6),2)); %looking for sweeps 1:6 (lecb first, mecr second)
idxSweep2 = logical(sum(temp == string(6:11),2)); %looking for sweeps 6:11 (lecb second, mecr first)
idxSweep3 = logical(sum(temp2 == string(1:6),2)); %mecb first, lecr second
idxSweep4 = logical(sum(temp2 == string(6:11),2)); %mecb second, lecr first

%store all mec second response data
timeCompile.mec = [tempData(idxSweep,idxR); tempData2(idxSweep4, idxB)];
timeCompile.mecRatio = [tempRatio(idxSweep,idxR); tempRatio2(idxSweep4, idxB)];
timeCompile.mecsRatio = [tempsRatio(idxSweep,idxR); tempsRatio2(idxSweep4, idxB)];
timeCompile.mecsRatioVal = [tempsRatioVal(idxSweep,idxR); tempsRatioVal2(idxSweep4, idxB)];
timeCompile.mecDiff = [tempDiff(idxSweep,idxR); tempDiff2(idxSweep4, idxB)];
timeCompile.mecPath = [tempPath(idxSweep,:); tempPath2(idxSweep4, :)];
timeCompile.mecbrRatio = [tempbrRatio(idxSweep,:); tempbrRatio(idxSweep4,:)];


%store all mec first response data 
timeCompile.mecFirst = [tempData(idxSweep2,idxR); tempData2(idxSweep3,idxB)];
timeCompile.mecRatioFirst = [tempRatio(idxSweep2,idxR);tempRatio2(idxSweep3,idxB)];
timeCompile.mecsRatioFirst = [tempsRatio(idxSweep2,idxR);tempsRatio2(idxSweep3,idxB)];
timeCompile.mecsRatioValFirst = [tempsRatioVal(idxSweep2,idxR);tempsRatioVal2(idxSweep3,idxB)];
timeCompile.mecDiffFirst = [tempDiff(idxSweep2,idxR);tempDiff2(idxSweep3,idxB)];
timeCompile.mecPathFirst = [tempPath(idxSweep2,:); tempPath2(idxSweep3,:)];


%store all lec second response data
timeCompile.lec = [tempData(idxSweep2,idxB); tempData2(idxSweep3,idxR)];
timeCompile.lecRatio =[tempRatio(idxSweep2,idxB); tempRatio2(idxSweep3,idxR)];
timeCompile.lecsRatio =[tempsRatio(idxSweep2,idxB); tempsRatio2(idxSweep3,idxR)];
timeCompile.lecsRatioVal =[tempsRatioVal(idxSweep2,idxB); tempsRatioVal2(idxSweep3,idxR)];
timeCompile.lecDiff = [tempDiff(idxSweep2,idxB);tempDiff2(idxSweep3,idxR)];
timeCompile.lecPath = [tempPath(idxSweep2,:);tempPath2(idxSweep3,:)];
timeCompile.lecbrRatio = [tempbrRatio(idxSweep2,:); tempbrRatio(idxSweep3,:)];


%store all lec first response data
timeCompile.lecFirst = [tempData(idxSweep,idxB); tempData2(idxSweep4,idxR)];
timeCompile.lecRatioFirst = [tempRatio(idxSweep,idxB);tempRatio2(idxSweep4,idxR)];
timeCompile.lecsRatioFirst = [tempsRatio(idxSweep,idxB);tempsRatio2(idxSweep4,idxR)];
timeCompile.lecsRatioValFirst = [tempsRatioVal(idxSweep,idxB);tempsRatioVal2(idxSweep4,idxR)];
timeCompile.lecDiffFirst = [tempDiff(idxSweep,idxB); tempDiff2(idxSweep4,idxR)];
timeCompile.lecPathFirst = [tempPath(idxSweep,:);tempPath2(idxSweep4,:)];

%Change all sweeps 7:10 to 5:1
timeAvg.lecPath = [timeCompile.lecPath(:,1:2), replace(timeCompile.lecPath(:,2), 'v11', 'v1')];
timeAvg.lecPath(:,3) = replace(timeAvg.lecPath(:,3), 'v10', 'v2');
timeAvg.lecPath(:,3) = replace(timeAvg.lecPath(:,3), 'v9', 'v3');
timeAvg.lecPath(:,3) = replace(timeAvg.lecPath(:,3), 'v8', 'v4');
timeAvg.lecPath(:,3) = replace(timeAvg.lecPath(:,3), 'v7', 'v5');

timeAvg.mecPath = [timeCompile.mecPath(:,1:2), replace(timeCompile.mecPath(:,2), 'v11', 'v1')];
timeAvg.mecPath(:,3) = replace(timeAvg.mecPath(:,3), 'v10', 'v2');
timeAvg.mecPath(:,3) = replace(timeAvg.mecPath(:,3), 'v9', 'v3');
timeAvg.mecPath(:,3) = replace(timeAvg.mecPath(:,3), 'v8', 'v4');
timeAvg.mecPath(:,3) = replace(timeAvg.mecPath(:,3), 'v7', 'v5');

%Change all sweeps 7:10 to 5:1
timeAvg.lecPathFirst = [timeCompile.lecPathFirst(:,1:2), replace(timeCompile.lecPathFirst(:,2), 'v11', 'v1')];
timeAvg.lecPathFirst(:,3) = replace(timeAvg.lecPathFirst(:,3), 'v10', 'v2');
timeAvg.lecPathFirst(:,3) = replace(timeAvg.lecPathFirst(:,3), 'v9', 'v3');
timeAvg.lecPathFirst(:,3) = replace(timeAvg.lecPathFirst(:,3), 'v8', 'v4');
timeAvg.lecPathFirst(:,3) = replace(timeAvg.lecPathFirst(:,3), 'v7', 'v5');

timeAvg.mecPathFirst = [timeCompile.mecPathFirst(:,1:2), replace(timeCompile.mecPathFirst(:,2), 'v11', 'v1')];
timeAvg.mecPathFirst(:,3) = replace(timeAvg.mecPathFirst(:,3), 'v10', 'v2');
timeAvg.mecPathFirst(:,3) = replace(timeAvg.mecPathFirst(:,3), 'v9', 'v3');
timeAvg.mecPathFirst(:,3) = replace(timeAvg.mecPathFirst(:,3), 'v8', 'v4');
timeAvg.mecPathFirst(:,3) = replace(timeAvg.mecPathFirst(:,3), 'v7', 'v5');

timeAvg.mecPath(:,4) = extractAfter(timeAvg.mecPath(:,3), 'b');
timeAvg.lecPath(:,4) = extractAfter(timeAvg.lecPath(:,3), 'b');
timeAvg.mecPathFirst(:,4) = extractAfter(timeAvg.mecPathFirst(:,3), 'b');
timeAvg.lecPathFirst(:,4) = extractAfter(timeAvg.lecPathFirst(:,3), 'b');

%lec data
[u,~,g] = unique(timeAvg.lecPath(:,4),'stable');
timeAvg.dataLec = splitapply(@(x) mean(x, 1,'omitnan'), timeCompile.lec, g);
timeAvg.dataLecSd = splitapply(@(x) std(x, 1,'omitnan'), timeCompile.lec, g);
timeAvg.lecNum =  histc(g, unique(g));
timeAvg.dataLecSem = timeAvg.dataLecSd ./ sqrt(timeAvg.lecNum); 
timeAvg.dataAvgLecPath = u;

timeAvg.ratioLec = splitapply(@(x) mean(x, 1,'omitnan'), timeCompile.lecRatio, g);
timeAvg.ratioLecSd = splitapply(@(x) std(x, 1,'omitnan'), timeCompile.lecRatio, g);
timeAvg.ratioLecSem = timeAvg.ratioLecSd ./ sqrt(timeAvg.lecNum);

timeAvg.diffLec = splitapply(@(x) mean(x, 1,'omitnan'), timeCompile.lecDiff, g);
timeAvg.diffLecSd = splitapply(@(x) std(x, 1,'omitnan'), timeCompile.lecDiff, g);
timeAvg.diffLecSem = timeAvg.diffLecSd ./ sqrt(timeAvg.lecNum);

timeAvg.brRatioLec = splitapply(@(x) mean(x,1,'omitnan'), timeCompile.lecbrRatio, g);
timeAvg.brRatioLecSd = splitapply(@(x) std(x, 1,'omitnan'), timeCompile.lecbrRatio, g);
timeAvg.brRatioLecSem = timeAvg.brRatioLecSd ./ sqrt(timeAvg.lecNum);

timeAvg.stimLec = u;

%lec first data
[u,~,g] = unique(timeAvg.lecPathFirst(:,4),'stable');
timeAvg.dataLecFirst = splitapply(@(x) mean(x, 1,'omitnan'), timeCompile.lecFirst, g);
timeAvg.dataLecFirstSd = splitapply(@(x) std(x, 1,'omitnan'), timeCompile.lecFirst, g);
timeAvg.lecFirstNum =  histc(g, unique(g));
timeAvg.dataLecFirstSem = timeAvg.dataLecFirstSd ./ sqrt(timeAvg.lecFirstNum); 
timeAvg.dataAvgLecFirstPath = u;

timeAvg.ratioLecFirst = splitapply(@(x) mean(x, 1,'omitnan'), timeCompile.lecRatioFirst, g);
timeAvg.ratioLecFirstSd = splitapply(@(x) std(x, 1,'omitnan'), timeCompile.lecRatioFirst, g);
timeAvg.ratioLecFirstSem = timeAvg.ratioLecFirstSd ./ sqrt(timeAvg.lecFirstNum);

timeAvg.diffLecFirst = splitapply(@(x) mean(x, 1,'omitnan'), timeCompile.lecDiffFirst, g);
timeAvg.diffLecFirstSd = splitapply(@(x) std(x, 1,'omitnan'), timeCompile.lecDiffFirst, g);
timeAvg.diffLecFirstSem = timeAvg.diffLecFirstSd ./ sqrt(timeAvg.lecFirstNum);

timeAvg.stimLecFirst = u;

%mec data
[u,~,g] = unique(timeAvg.mecPath(:,4),'stable');
timeAvg.dataMec = splitapply(@(x) mean(x, 1,'omitnan'), timeCompile.mec, g);
timeAvg.dataMecSd = splitapply(@(x) std(x, 1,'omitnan'), timeCompile.mec, g);
timeAvg.mecNum =  histc(g, unique(g));
timeAvg.dataMecSem = timeAvg.dataMecSd ./ sqrt(timeAvg.mecNum); 
timeAvg.dataAvgMecPath = u;

timeAvg.ratioMec = splitapply(@(x) mean(x, 1,'omitnan'), timeCompile.mecRatio, g);
timeAvg.ratioMecSd = splitapply(@(x) std(x, 1,'omitnan'), timeCompile.mecRatio, g);
timeAvg.ratioMecSem = timeAvg.ratioMecSd ./ sqrt(timeAvg.mecNum);

timeAvg.diffMec = splitapply(@(x) mean(x, 1,'omitnan'), timeCompile.mecDiff, g);
timeAvg.diffMecSd = splitapply(@(x) std(x, 1,'omitnan'), timeCompile.mecDiff, g);
timeAvg.diffMecSem = timeAvg.diffMecSd ./ sqrt(timeAvg.mecNum);

timeAvg.brRatioMec = splitapply(@(x) mean(x,1,'omitnan'), timeCompile.mecbrRatio, g);
timeAvg.brRatioMecSd = splitapply(@(x) std(x, 1,'omitnan'), timeCompile.mecbrRatio, g);
timeAvg.brRatioMecSem = timeAvg.brRatioMecSd ./ sqrt(timeAvg.mecNum);

timeAvg.stimMec = u;

%mec first data
[u,~,g] = unique(timeAvg.mecPathFirst(:,4),'stable');
timeAvg.dataMecFirst = splitapply(@(x) mean(x, 1,'omitnan'), timeCompile.mecFirst, g);
timeAvg.dataMecFirstSd = splitapply(@(x) std(x, 1,'omitnan'), timeCompile.mecFirst, g);
timeAvg.mecFirstNum =  histc(g, unique(g));
timeAvg.dataMecFirstSem = timeAvg.dataMecSd ./ sqrt(timeAvg.mecFirstNum); 
 imeAvg.dataAvgMecFirstPath = u;

timeAvg.ratioMecFirst = splitapply(@(x) mean(x, 1,'omitnan'), timeCompile.mecRatioFirst, g);
timeAvg.ratioMecFirstSd = splitapply(@(x) std(x, 1,'omitnan'), timeCompile.mecRatioFirst, g);
timeAvg.ratioMecFirstSem = timeAvg.ratioMecFirstSd ./ sqrt(timeAvg.mecFirstNum);

timeAvg.diffMecFirst = splitapply(@(x) mean(x, 1,'omitnan'), timeCompile.mecDiffFirst, g);
timeAvg.diffMecFirstSd = splitapply(@(x) std(x, 1,'omitnan'), timeCompile.mecDiffFirst, g);
timeAvg.diffMecFirstSem = timeAvg.diffMecFirstSd ./ sqrt(timeAvg.mecFirstNum);

timeAvg.stimMecFirst = u;

timeAvg.dataKey = timeCompile.dataKey(1,idxR);



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
shadedErrorBar(1:6', timeAvg.ratioLec(1:6, 10), timeAvg.ratioLecSem(1:6, 10),'lineprops', {'b','LineWidth',2}); hold on;
shadedErrorBar(1:6', fliplr(timeAvg.ratioMec(1:6, 10)'), fliplr(timeAvg.ratioMecSem(1:6, 10)'),'lineprops', {'r','LineWidth',2});

ylabel('Summation Ratio 2nd/1st'); xlabel('Time between Input Stim (ms)');
legend('LEC','MEC');
title('Ratios: Low Intensity Timing Experiment at -70');
xticks([1:6]);
xticklabels({ '0','10','20','30','40','50'})
set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
xlim([1,6]);

subplot(1,3,2); %amplitude
shadedErrorBar(1:6, timeAvg.ratioLec(1:6, 8), timeAvg.ratioLecSem(1:6, 8),'lineprops', {'b','LineWidth',2}); hold on;
shadedErrorBar(1:6, fliplr(timeAvg.ratioMec(1:6, 8)'), fliplr(timeAvg.ratioMecSem(1:6, 8)'),'lineprops', {'r','LineWidth',2});

ylabel('Amplitude Ratio 2nd/1st'); xlabel('Time between Input Stim (ms)');
legend('LEC','MEC');
xticks([1:6]);
xticklabels({ '0','10','20','30','40','50'})
set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
xlim([1,6]);

subplot(1,3,3); %peak
shadedErrorBar(1:6, timeAvg.ratioLec(1:6, 6), timeAvg.ratioLecSem(1:6, 6),'lineprops', {'b','LineWidth',2}); hold on;
shadedErrorBar(1:6, fliplr(timeAvg.ratioMec(1:6, 6)'), fliplr(timeAvg.ratioMecSem(1:6, 6)'),'lineprops', {'r','LineWidth',2});

ylabel('Peak Ratio 2nd/1st'); xlabel('Time between Input Stim (ms)');
legend('LEC', 'MEC');
xticks([1:6]);
xticklabels({ '0','10','20','30','40','50'})
set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
xlim([1,6]);;

%% Difference
figure; subplot(1,3,1); %summation
shadedErrorBar(1:6, timeAvg.diffLec(1:6, 10), timeAvg.diffLecSem(1:6, 10),'lineprops', {'b','LineWidth',2}); hold on;
shadedErrorBar(1:6, fliplr(timeAvg.diffMec(1:6, 10)'), fliplr(timeAvg.diffMecSem(1:6, 10)'),'lineprops', {'r','LineWidth',2});

ylabel('Summation Diff 2nd-1st'); xlabel('Time between Input Stim (ms)');
legend('LEC','MEC');
title('Difference: Low Intensity Timing Experiment at -70');
xticks([1:6]);
xticklabels({ '0','10','20','30','40','50'})
set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
xlim([1,6]);

subplot(1,3,2); %amplitude
shadedErrorBar(1:6, timeAvg.diffLec(1:6, 8), timeAvg.diffLecSem(1:6, 8),'lineprops', {'b','LineWidth',2}); hold on;
shadedErrorBar(1:6, fliplr(timeAvg.diffMec(1:6, 8)'), fliplr(timeAvg.diffMecSem(1:6, 8)'),'lineprops', {'r','LineWidth',2});

ylabel('Amplitude Diff 2nd-1st'); xlabel('Time between Input Stim (ms)');
legend('LEC', 'MEC');
xticks([1:6]);
xticklabels({ '0','10','20','30','40','50'})
set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
xlim([1,6]);

subplot(1,3,3); %peak
shadedErrorBar(1:6, timeAvg.diffLec(1:6, 6), timeAvg.diffLecSem(1:6, 6),'lineprops', {'b','LineWidth',2}); hold on;
shadedErrorBar(1:6, fliplr(timeAvg.diffMec(1:6, 6)'), fliplr(timeAvg.diffMecSem(1:6, 6)'),'lineprops', {'r','LineWidth',2});

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

 save(['timeMat'], 'timeAvg', 'timeCompile');