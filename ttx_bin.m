%ttx binned data 
%bin ttx io data
clear all; close all; 

path = ['/Volumes/KO_Portable/MEC_LEC_CA3/analysis/VC/TTX/io']; cd(path); 
path2 =  ['/Volumes/KO_Portable/MEC_LEC_CA3/data_converted'];
path3 = ['/Volumes/KO_Portable/MEC_LEC_CA3/analysis/mW'];

cd(path);
load('ttx.mat');
cd(path3); 
load('led_power.mat');

%Find Values for Binned Data 
ttx.binCell(:,1) = extractBefore(ttx.sweepCell(:,2), '_');

for i = 1:length(ttx.sweepCell(:,2))
    if contains(ttx.sweepCell(i,2), 'iof') && contains(extractBefore(ttx.sweepCell(i,2), '_'), '10')
       ttx.binCell(i,2) = 1;
    elseif contains(ttx.sweepCell(i,2), 'iof')
       ttx.binCell(i,2) = str2double(strcat('0.',extractBefore(ttx.sweepCell(i,2), '_')));
    elseif contains(ttx.sweepCell(i,2), 'ioc') 
       ttx.binCell(i,2) = extractBefore(ttx.sweepCell(i,2), '_');
    end
    ttx.binCell(i,3) = blue(find(ttx.binCell(i,2) == string(blue(:,1))),4);
    ttx.binCell(i,4) = red(find(ttx.binCell(i,2) == string(red(:,1))),4);
end

ttx.binCellId(:,1) = strcat(string(ttx.binCell(:,3)), 'K', ttx.sweepCell(:,3));
ttx.binCellId(:,2) = strcat(string(ttx.binCell(:,4)), 'K', ttx.sweepCell(:,3));
ttx.binCellId(find(ismissing(ttx.binCellId))) = string(0);

%Bin Data by Cell using ttx.binCell
[u, ~, g] = unique(ttx.binCellId(:,1),'stable');
ttx.dataBinCellBlue = splitapply(@(x) mean(x,1, 'omitnan'), ttx.dataCell(:,1), g);
ttx.sweepBinCellBlue = u; 

[u, ~, g] = unique(ttx.binCellId(:,2),'stable');
ttx.dataBinCellRed = splitapply(@(x) mean(x,1, 'omitnan'), ttx.dataCell(:,2), g);
ttx.sweepBinCellRed = u;

ttx.dataBinCellBlue(find(ttx.sweepBinCellBlue == '0')) = [];
ttx.sweepBinCellBlue(find(ttx.sweepBinCellBlue == '0')) = [];
ttx.dataBinCellRed(find(ttx.sweepBinCellRed == '0')) = [];
ttx.sweepBinCellRed(find(ttx.sweepBinCellRed == '0')) = [];


cd(path2);
fold = dir('*202*');
tempFields = ["sweepBinCellRed", "sweepBinCellBlue"];

for i = 1:2
    anNum = extractBetween(ttx.(tempFields(i)), 'K', '_S');
    
    for ii = 1:length(anNum)
        cd(fold(contains(extractfield(fold,'name'), anNum(ii))).name);
        fold_temp = dir('*202*');
        cd(fold_temp(contains(extractfield(fold_temp,'name'), extractAfter(ttx.(tempFields(i))(ii,1), 'K'))).name);
        pathway = readcell('pathway');
        if contains(tempFields(i), 'Blue')
            ttx.(tempFields(i)){ii,2} = pathway{string(pathway)=='blue',2};
        else 
            ttx.(tempFields(i)){ii,2} = pathway{string(pathway)=='red',2};
        end
        cd(path2);
    end
end

ttx.sweepBinCellBlue(:,2) = strcat(extractBefore(ttx.sweepBinCellBlue(:,1),'K'), ttx.sweepBinCellBlue(:,2));
ttx.sweepBinCellRed(:,2) = strcat(extractBefore(ttx.sweepBinCellRed(:,1),'K'), ttx.sweepBinCellRed(:,2));

[u, ~, g] = unique(ttx.sweepBinCellRed(:,2),'stable');
ttxBinAvg.dataRed = splitapply(@(x) mean(x,1, 'omitnan'), ttx.dataBinCellRed(:,1), g);
ttxBinAvg.sdRed = splitapply(@(x) std(x,1, 'omitnan'), ttx.dataBinCellRed(:,1), g);
ttxBinAvg.numRed = histc(g, unique(g));
ttxBinAvg.semRed = ttxBinAvg.sdRed ./ sqrt(ttxBinAvg.numRed)
ttxBinAvg.sweepRed = u;

[u, ~, g] = unique(ttx.sweepBinCellBlue(:,2),'stable');
ttxBinAvg.dataBlue = splitapply(@(x) mean(x,1, 'omitnan'), ttx.dataBinCellBlue(:,1), g);
ttxBinAvg.sdBlue = splitapply(@(x) std(x,1, 'omitnan'), ttx.dataBinCellBlue(:,1), g);
ttxBinAvg.numBlue = histc(g, unique(g));
ttxBinAvg.semBlue = ttxBinAvg.sdBlue ./ sqrt(ttxBinAvg.numBlue)
ttxBinAvg.sweepBlue = u;

ttx.dataBinCom = [ttx.dataBinCellRed(:,1); ttx.dataBinCellBlue(:,1)];
ttx.sweepBinCellCom = [ttx.sweepBinCellRed(:,2); ttx.sweepBinCellBlue(:,2)]

[u, ~, g] = unique(ttx.sweepBinCellCom(:,1),'stable');
ttxBinAvg.dataCom = splitapply(@(x) mean(x,1, 'omitnan'), ttx.dataBinCom(:,1), g);
ttxBinAvg.sdCom = splitapply(@(x) std(x,1, 'omitnan'), ttx.dataBinCom(:,1), g);
ttxBinAvg.numCom = histc(g, unique(g));
ttxBinAvg.semCom = ttxBinAvg.sdCom ./ sqrt(ttxBinAvg.numCom)
ttxBinAvg.sweepCom = u;

save(['ttx_bin.mat'], 'ttx', 'ttxBinAvg', 'maxData','minData');
%Plot Figures

figure;
errorbar(mW_bin(1:7,4), ttxBinAvg.dataCom(1:7), ttxBinAvg.semCom(1:7), 'r.-', 'MarkerSize', 40, 'LineWidth',2); hold on;
errorbar(mW_bin(1:7,4), ttxBinAvg.dataCom(16:22), ttxBinAvg.semCom(16:22), 's-', 'MarkerSize', 10, 'LineWidth',2,'Color', '#289e1e', 'MarkerFaceColor', '#289e1e')

xlim([1 12]); ylim([-90 0]); ylabel('EPSC (pA)'); xlabel('LED Intensity (%)');
set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);


%Find indices for individual min/max
tempI = find(matches(ttx.sweepBinCellBlue(:,2), '1lec'));
tempI2 = find(matches(ttx.sweepBinCellRed(:,2), '1lec'));

indvBinCellMat = [ttx.dataBinCellBlue(tempI); ttx.dataBinCellRed(tempI2)];
indvBinCellMatSweep = [join(ttx.sweepBinCellBlue(tempI, 1:2)); join(ttx.sweepBinCellRed(tempI2,1:2))];

tempI = find(matches(ttx.sweepBinCellRed(:,2), '1mec'));
tempI2 = find(matches(ttx.sweepBinCellBlue(:,2), '1mec'));
indvBinCellMat = [indvBinCellMat, [ttx.dataBinCellRed(tempI); ttx.dataBinCellBlue(tempI2)]];
indvBinCellMatSweep = [indvBinCellMatSweep, [join(ttx.sweepBinCellRed(tempI, 1:2)); join(ttx.sweepBinCellBlue(tempI2, 1:2))]];

tempI = find(matches(ttx.sweepBinCellBlue(:,2), '7lec'));
tempI2 = find(matches(ttx.sweepBinCellRed(:,2), '7lec'));

indvBinCellMat = [indvBinCellMat, [ttx.dataBinCellBlue(tempI); ttx.dataBinCellRed(tempI2)]];
indvBinCellMatSweep = [indvBinCellMatSweep, [join(ttx.sweepBinCellBlue(tempI, 1:2)); join(ttx.sweepBinCellRed(tempI2,1:2))]];

tempI = find(matches(ttx.sweepBinCellRed(:,2), '7mec'));
tempI2 = find(matches(ttx.sweepBinCellBlue(:,2), '7mec'));

indvBinCellMat = [indvBinCellMat, [ttx.dataBinCellRed(tempI); ttx.dataBinCellBlue(tempI2)]];
indvBinCellMatSweep = [indvBinCellMatSweep, [join(ttx.sweepBinCellRed(tempI, 1:2)); join(ttx.sweepBinCellBlue(tempI2,1:2))]];


figure;
plot(1, indvBinCellMat(:,1), 's', 'LineWidth',1, 'Color', '#808080', 'MarkerFaceColor', '#808080'); hold on;
plot(2, indvBinCellMat(:,2), 'ko', 'LineWidth',1,'Color', '#808080', 'MarkerFaceColor', '#808080');
plot(3, indvBinCellMat(:,3), 'sk', 'LineWidth',1, 'Color', '#808080','MarkerFaceColor', '#808080');
plot(4, indvBinCellMat(:,4), 'ko', 'LineWidth',1, 'Color', '#808080','MarkerFaceColor', '#808080');

plot(1:2, indvBinCellMat(:,1:2), '.-','LineWidth',1, 'Color', '#808080');
plot(3:4, indvBinCellMat(:,3:4), '.-','LineWidth',1,'Color', '#808080');

errorbar(2, ttxBinAvg.dataCom(1), ttxBinAvg.semCom(1), 'ro', 'MarkerSize', 10, 'LineWidth',3);
errorbar(1, ttxBinAvg.dataCom(16), ttxBinAvg.semCom(16), 's', 'MarkerSize', 10, 'LineWidth',3,'Color', '#289e1e', 'MarkerFaceColor', '#289e1e')
errorbar(4, ttxBinAvg.dataCom(7), ttxBinAvg.semCom(7), 'ro', 'MarkerSize', 10, 'LineWidth',3, 'MarkerFaceColor' ,'r'); hold on;
errorbar(3, ttxBinAvg.dataCom(22), ttxBinAvg.semCom(22), 's', 'MarkerSize', 10, 'LineWidth',3,'Color', '#289e1e', 'MarkerFaceColor', '#289e1e')

xlim([0.75 4.25]); xticks([1 2 3 4]); xticklabels({'Min LEC', 'Min MEC', 'Max LEC', 'Max MEC'});
set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
ylabel('EPSC (pA)')