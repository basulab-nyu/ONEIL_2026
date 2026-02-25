
clear all; close all; 

path = ['/Volumes/KO_Portable/MEC_LEC_CA3/analysis/VC/TTX/io']; cd(path); 
path2 =  ['/Volumes/KO_Portable/MEC_LEC_CA3/data_converted'];
path3 = ['/Volumes/KO_Portable/MEC_LEC_CA3/analysis/mW'];

files = dir('*_io80ttx*.csv');
files = files(~startsWith(extractfield(files, 'name'), ['.']));
kinetics = ["amp"];
files_kinetics = files(contains(extractfield(files, 'name'), kinetics));

ttx_tempTbl = readtable(files_kinetics(1).name);
led = extractBetween(ttx_tempTbl.Properties.VariableNames(2:end), 'ttx_', '_');
ttx_temp = ttx_tempTbl{:,2:end};
ttx_temp(find(ttx_temp >=-1)) = NaN;
sweeps = extractAfter(ttx_tempTbl{:,1}, '_1_');

ttx.dataAll = ttx_temp;
ttx.sweepAll = sweeps;
ttx.ledAll = led;

temp_idx = [find(matches(led, 'bb'),1,'first') find(matches(led, 'br'),1,'first')];
temp_idx2 =[find(matches(led, 'rr'),1,'first') find(matches(led, 'rb'),1,'first')];
temp_idx3 = find(matches(led, 'brbr'),1,'first');
ttx_temp2 = [mean(ttx_temp(:,temp_idx), 2, 'omitnan'), mean(ttx_temp(:,temp_idx2), 2, 'omitnan'), ttx_temp(:,temp_idx3)];
ttx.ledCell = ["bb" "rr" "combo"];

cd(path3); 
load('led_power.mat');

[u, ~, g] = unique(sweeps,'stable');
ttx.dataCell = splitapply(@(x) mean(x,1, 'omitnan'), ttx_temp2, g);
ttx.sweepCell = string(u);

temp_sweep = erase(ttx.sweepCell, 'TTX'); temp_sweep = erase(temp_sweep, 'ttx');
temp_sweep = erase(temp_sweep, '80'); temp_sweep = erase(temp_sweep, '_drug');
temp_sweep = erase(temp_sweep, '_srcgp_'); temp_sweep = extractBefore(temp_sweep, '20');
ttx.sweepCell(:,2) = temp_sweep;
ttx.sweepCell(:,3) =  extractAfter(string(u), '_K');

cd(path2);
fold = dir('*202*');
anNum = extractBefore(ttx.sweepCell(:,3), '_S');

for i = 1:length(anNum)
    cd(fold(contains(extractfield(fold,'name'), anNum(i))).name);
    fold_temp = dir('*202*');
    cd(fold_temp(contains(extractfield(fold_temp,'name'), ttx.sweepCell(i,3))).name);
    pathway = readcell('pathway');
    ttx.blueCellPath{i,1} = pathway{string(pathway)=='blue',2};
    ttx.redCellPath{i,1} = pathway{string(pathway)=='red',2};
    cd(path2);
end

% %Find values for binned data
% ttx.binCell(:,1) = extractBefore(ttx.sweepCell(:,2), '_');
% 
% for i = 1:length(ttx.sweepCell(:,2))
%     if contains(ttx.sweepCell(i,2), 'iof') && contains(extractBefore(ttx.sweepCell(i,2), '_'), '10')
%        ttx.binCell(i,2) = 1;
%     elseif contains(ttx.sweepCell(i,2), 'iof')
%        ttx.binCell(i,2) = str2double(strcat('0.',extractBefore(ttx.sweepCell(i,2), '_')));
%     elseif contains(ttx.sweepCell(i,2), 'ioc') 
%        ttx.binCell(i,2) = extractBefore(ttx.sweepCell(i,2), '_');
%     end
%     ttx.binCell(i,3) = blue(find(ttx.binCell(i,2) == string(blue(:,1))),4);
%     ttx.binCell(i,4) = red(find(ttx.binCell(i,2) == string(red(:,1))),4);
% end
% 
% %Bin Data by Cell using ttx.binCell
% 
% ttx.binCellId(:,1) = strcat(string(ttx.binCell(:,3)), 'K', ttx.sweepCell(:,3));
% ttx.binCellId(:,2) = strcat(string(ttx.binCell(:,4)), 'K', ttx.sweepCell(:,3));
% ttx.binCellId(find(ismissing(ttx.binCellId))) = string(0);
% 
% %Bin Data by Cell using ttx.binCell
% [u, ~, g] = unique(ttx.binCellId(:,1),'stable');
% ttx.dataBinCellBlue = splitapply(@(x) mean(x,1, 'omitnan'), ttx.dataCell(:,1), g);
% ttx.sweepBinCellBlue = u; 
% 
% [u, ~, g] = unique(ttx.binCellId(:,2),'stable');
% ttx.dataBinCellRed = splitapply(@(x) mean(x,1, 'omitnan'), ttx.dataCell(:,2), g);
% ttx.sweepBinCellRed = u;
% 
% ttx.dataBinCellBlue(find(ttx.sweepBinCellBlue == '0')) = [];
% ttx.sweepBinCellBlue(find(ttx.sweepBinCellBlue == '0')) = [];
% ttx.dataBinCellRed(find(ttx.sweepBinCellRed == '0')) = [];
% ttx.sweepBinCellRed(find(ttx.sweepBinCellRed == '0')) = [];


ttx.sweepCell(:,4) = append(ttx.sweepCell(:,2), ttx.blueCellPath);
[u, ~, g] = unique(ttx.sweepCell(:,4),'stable');
ttx.dataAvg = splitapply(@(x) mean(x,1, 'omitnan'), ttx.dataCell, g);
ttx.sdAvg = splitapply(@(x) std(x,1, 'omitnan'), ttx.dataCell, g);
ttx.numAvg = histc(g, unique(g));
ttx.semAvg = ttx.sdAvg ./ sqrt(ttx.numAvg);
ttx.sweepAvg = u;

tempIdx = find(contains(ttx.sweepCell(:,4), 'lec'));
ttx.dataEc(tempIdx,1) = ttx.dataCell(tempIdx, 1);
ttx.dataEc(tempIdx,2) = ttx.dataCell(tempIdx, 2);
ttx.sweepEc(tempIdx,1:2) = ttx.sweepCell(tempIdx,[2,4]);

tempIdx = find(contains(ttx.sweepCell(:,4), 'mec'));
ttx.dataEc(tempIdx,1) = ttx.dataCell(tempIdx, 2);
ttx.dataEc(tempIdx,2) = ttx.dataCell(tempIdx, 1);
ttx.sweepEc(tempIdx,1:2) = ttx.sweepCell(tempIdx,[2,4]);

[u, ~, g] = unique(ttx.sweepEc(:,1),'stable');
ttx.dataEcAvg = splitapply(@(x) mean(x,1, 'omitnan'), ttx.dataEc, g);
ttx.dataEcAvg = splitapply(@(x) mean(x,1, 'omitnan'), ttx.dataEc, g);
ttx.sdEcAvg = splitapply(@(x) std(x,1, 'omitnan'), ttx.dataEc, g);
ttx.numEcAvg = histc(g, unique(g));
ttx.semEcAvg = ttx.sdEcAvg ./ sqrt(ttx.numEcAvg);
ttx.sweepEcAvg = u;


cd('/Volumes/KO_Portable/MEC_LEC_CA3/analysis/VC/IO');
load('Compile_Kinetics.mat');

%Get Max Data
maxData = [];
rowIdx = find(contains(ttx.sweepAvg, 'lec') & contains(ttx.sweepAvg, '10_ioc'));
maxData.dataAvgLecbMecr = ttx.dataAvg(rowIdx, :);
maxData.semAvgLecbMecr = ttx.semAvg(rowIdx,:);
rowIdx = find(contains(ttx.sweepCell(:,4), 'lec') & contains(ttx.sweepCell(:,4), '10_ioc'));
maxData.dataLecbMecr = ttx.dataCell(rowIdx,:);

rowIdx = find(contains(ttx.sweepEcAvg, '10_ioc'));
maxData.dataAvgEc = ttx.dataEcAvg(rowIdx,:);
maxData.semAvgEc = ttx.semAvg(rowIdx,:);
rowIdx = find(contains(ttx.sweepEc(:,1), '10_ioc'));
maxData.dataEc = ttx.dataEc(rowIdx,:);

rowIdx = find(io.lecBDualStim{2} == 10);
acsfMaxData = io.lecBDual{2}(rowIdx, 1); 
rowIdx = find(io.mecRDualStim{2} == 10);
acsfMaxData(:,2) = io.mecRDual{2}(rowIdx, 1);
rowIdx = find(io.comboStim{2} == 10);
acsfMaxDataCombo = io.combo{2}(rowIdx, 1); 

%Get Min Data
minData = [];
rowIdx = find(contains(ttx.sweepAvg, 'lec') & contains(ttx.sweepAvg, '1_iof'));
minData.dataAvgLecbMecr = ttx.dataAvg(rowIdx, :);
minData.semAvgLecbMecr = ttx.semAvg(rowIdx,:);
rowIdx = find(contains(ttx.sweepCell(:,4), 'lec') & contains(ttx.sweepCell(:,4), '1_iof'));
minData.dataLecbMecr = ttx.dataCell(rowIdx,:);

rowIdx = find(contains(ttx.sweepEcAvg, '1_iof'));
minData.dataAvgEc = ttx.dataEcAvg(rowIdx,:);
minData.semAvgEc = ttx.semAvg(rowIdx,:);
rowIdx = find(contains(ttx.sweepEc(:,1), '1_iof'));
minData.dataEc = ttx.dataEc(rowIdx,:);

 cd(path); save(['ttx'], 'ttx', 'minData','maxData');

figure;
plot(1, maxData.dataLecbMecr(:,1), 'ks','LineWidth',1,'MarkerFaceColor', 'k'); hold on;
plot(2, maxData.dataLecbMecr(:,2), 'ko', 'LineWidth',1,'MarkerFaceColor', 'k');
plot(3, maxData.dataLecbMecr(:,3), 'k^','LineWidth',1,'MarkerFaceColor', 'k');
plot(4, acsfMaxData(:,1), 'ks', 'LineWidth', 1, 'MarkerFaceColor', 'k');
plot(5, acsfMaxData(:,2), 'ko', 'LineWidth', 1, 'MarkerFaceColor', 'k');
plot(6, acsfMaxDataCombo(:,1), 'k^', 'LineWidth', 1, 'MarkerFaceColor', 'k');

plot(1:3, maxData.dataAvgLecbMecr, 'k.-');
plot(4:5, acsfMaxData, 'k.-');
errorbar(1, maxData.dataAvgLecbMecr(1), maxData.semAvgLecbMecr(1), 's', 'MarkerSize', 15, 'LineWidth',2,'Color', '#289e1e', 'MarkerFaceColor', '#289e1e');
errorbar(2, maxData.dataAvgLecbMecr(2), maxData.semAvgLecbMecr(2), 'r.', 'MarkerSize', 40, 'LineWidth',2);
errorbar(3, maxData.dataAvgLecbMecr(3), maxData.semAvgLecbMecr(3), 'm^', 'MarkerSize', 15, 'MarkerFaceColor', 'm', 'LineWidth',2);
errorbar(4, ioAvg.lecBDual{2}(19,1), ioAvg.lecBDualSem{2}(19,1), 's', 'MarkerSize', 15, 'LineWidth',3,'Color', '#289e1e');
errorbar(5, ioAvg.mecRDual{2}(19,1), ioAvg.mecRDualSem{2}(19,1), 'ro', 'MarkerSize', 15, 'LineWidth', 3);
errorbar(6, ioAvg.combo{2}(19,1), ioAvg.comboSem{2}(19,1), 'm^', 'MarkerSize', 15, 'LineWidth', 3);

xlim([0.5 6.5]); xticks([1 2 3 4 5 6]); xticklabels({'LECB TTX', 'MECR TTX', 'EC TTX', 'LECB ACSF', 'MECR ACSF', 'EC ACSF'});
set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
ylabel('EPSC (pA)');

%EC Average, combined Opsin

rowIdx = find(contains(ttx.sweepEc(:,1), '7_iof'));
rowIdxAvg = find(contains(ttx.sweepEcAvg(:,1), '7_iof'));

figure;
plot(1, minData.dataEc(:,1), 'ks','LineWidth',1,'MarkerFaceColor', 'k'); hold on;
plot(2, minData.dataEc(:,2), 'ko', 'LineWidth',1,'MarkerFaceColor', 'k');
plot(3, ttx.dataEc(rowIdx,1), 'ks', 'LineWidth',1,'MarkerFaceColor', 'k');
plot(4, ttx.dataEc(rowIdx,2), 'ko', 'LineWidth', 1,'MarkerFaceColor','k');
%plot(3, maxData.dataEc(:,1), 'ks','LineWidth',1,'MarkerFaceColor', 'k'); 
%plot(4, maxData.dataEc(:,2), 'ko', 'LineWidth',1,'MarkerFaceColor', 'k');
plot(1:2, minData.dataEc(:,1:2), 'k.-');
plot(3:4, ttx.dataEc(rowIdx,1:2), 'k.-');
%plot(3:4, maxData.dataEc(:,1:2), 'k.-');

errorbar(1, minData.dataAvgEc(1), minData.semAvgEc(1), 's', 'MarkerSize', 15, 'LineWidth',2,'Color', '#289e1e', 'MarkerFaceColor', '#289e1e');
errorbar(2, minData.dataAvgEc(2), minData.semAvgEc(2), 'r.', 'MarkerSize', 40, 'LineWidth',2);
errorbar(3, ttx.dataEcAvg(rowIdxAvg, 1), ttx.semEcAvg(rowIdxAvg, 1),'s', 'MarkerSize', 15, 'LineWidth',2,'Color', '#289e1e', 'MarkerFaceColor', '#289e1e');
errorbar(4, ttx.dataEcAvg(rowIdxAvg, 2), ttx.semEcAvg(rowIdxAvg, 2), 'r.', 'MarkerSize', 40, 'LineWidth',2);
%errorbar(3, maxData.dataAvgEc(1), maxData.semAvgEc(1), 's', 'MarkerSize', 15, 'LineWidth',2,'Color', '#289e1e', 'MarkerFaceColor', '#289e1e');
%errorbar(4, maxData.dataAvgEc(2), maxData.semAvgEc(2), 'r.', 'MarkerSize', 40, 'LineWidth',2);


xlim([0.5 4.5]); xticks([1 2 3 4]); xticklabels({'Min LEC TTX', 'Min MEC TTX', 'Max LEC TTX', 'Max MEC TTX'});
set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
ylabel('EPSC (pA)');


