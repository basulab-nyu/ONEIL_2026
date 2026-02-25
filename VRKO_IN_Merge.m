%% Combine VR and KO Data
% 
clear all; close all;
addpath '/Volumes/KO_Portable/MATLAB/EphysAnalysis_2025' ;

path = '/Volumes/KO_Portable/MEC_LEC_CA3/analysis/VC/IN/VR_LEC_Data/';
cd(path);

files = dir('amp*tc*.csv');

for i = 1:length(files)
    tempTbl = table2array(readtable(files(i).name));
    cellId = readtable(files(i).name); cellId = string(cellId.Properties.VariableNames);
    idx = contains(cellId, '_d'); cellId = cellId(idx);
    cellId = extractBetween(cellId, 'C_', '_');
    cellId(2,:) = extractBefore(cellId, 'c');

    lecTc{i}.time = tempTbl(:,1:2:end);
    lecTc{i}.amp = tempTbl(:,2:2:end);

    cd ../ 

    in = string(extractBetween(files(i).name, 'tc_', '.csv'));
    lecTc{i}.in = in;
    temp = dir(strcat(in, '*.mat'));
    temp = load(temp.name); tempF = string(fields(temp)); temp = temp.(tempF);
    idx = find(string(temp.ptwy(:,1)) == "lec");

    tempMat = [string(temp.id(idx)); extractBefore(string(temp.id(idx)), '_')]
    cellId = [cellId, tempMat];
    lecTc{i}.cellId = cellId;
   
   
    for ii = 1:length(idx)
        row = length(cell2mat(temp.amp{idx(ii)}(2:end,1)));
        lecTc{i}.amp(1:row,end+1) = cell2mat(temp.amp{idx(ii)}(2:end,1));
        lecTc{i}.time(1:row, end+1) = temp.time{idx(ii)}*0.25';
    end
    
    [lecTc{i}.cno, ~] = arrayfun(@(c) find(lecTc{i}.time(:,c) == 0, 1, 'first'), 1:size(lecTc{i}.time,2));
    lecTc{i}.normAmp = arrayfun(@(x,y) mean(lecTc{i}.amp(y-20:y,x)), 1:size(lecTc{i}.amp,2), lecTc{i}.cno);
    lecTc{i}.normTime = arrayfun(@(x,y) lecTc{i}.time(y-20:end,x), 1:size(lecTc{i}.amp,2), lecTc{i}.cno, 'UniformOutput', false);
    tempVec = max(cellfun(@(x) length(x), lecTc{i}.normTime)); 
    tempMat = cellfun(@(x) [x; nan(tempVec - numel(x), 1)], lecTc{i}.normTime, 'UniformOutput', false);
    lecTc{i}.normTime = cell2mat(tempMat);
    lecTc{i}.normTime(find(lecTc{i}.normTime == 0)) = NaN; 

    lecTc{i}.normAmpTc = arrayfun(@(x,y) lecTc{i}.amp(y-20:end,x)./lecTc{i}.normAmp(:,x), 1:size(lecTc{i}.amp,2), lecTc{i}.cno,'UniformOutput',false);
    tempMat = cellfun(@(x) [x; nan(tempVec - numel(x), 1)], lecTc{i}.normAmpTc, 'UniformOutput', false);
    lecTc{i}.normAmpTc = cell2mat(tempMat);
    lecTc{i}.normAmpTc(find(lecTc{i}.normAmpTc == 0)) = NaN; 
 
    lecTc{i}.mNormAmpTc = mean(lecTc{i}.normAmpTc, 2, 'omitnan');
    lecTc{i}.mNormAmpTcN = length(lecTc{i}.cno);
    lecTc{i}.mNormAmpStd = std(lecTc{i}.normAmpTc,0,2,'omitnan');
    lecTc{i}.mNormAmpTcSem = lecTc{i}.mNormAmpStd ./ sqrt(lecTc{i}.mNormAmpTcN);

    lecTc{i}.mAmpIndv = lecTc{i}.normAmp;
    lecTc{i}.mAmpIndv(2,:) = arrayfun(@(x,y) mean(lecTc{i}.amp(y+25:y+35,x)), 1:size(lecTc{i}.amp,2), lecTc{i}.cno);
    lecTc{i}.mAmp = mean(lecTc{i}.mAmpIndv,2);
    lecTc{i}.stdAmp = std(lecTc{i}.mAmpIndv,0,2);
    lecTc{i}.nAmp = length(lecTc{i}.mAmpIndv);
    lecTc{i}.semAmp = lecTc{i}.stdAmp ./ sqrt(lecTc{i}.nAmp);

    [row,col] = size(lecTc{i}.normAmpTc);
    if abs(row - mod(row,4)) <= abs(row + (4 - mod(row,4)))
        idx = row - mod(row,4);
    else
        idx = row + (4 - mod(row,4));
    end
    for ii = 1:col
        lecTc{i}.avgAmpNorm(:,ii) = mean(reshape(lecTc{i}.normAmpTc(1:idx,ii), 4,[]), 'omitnan')' ;
    end
    
    lecTc{i}.mAvgAmpNorm= mean(lecTc{i}.avgAmpNorm,2,'omitnan');
    lecTc{i}.mAvgAmpNormStd = std(lecTc{i}.avgAmpNorm,0, 2,'omitnan');
    lecTc{i}.mAvgAmpNormN = cell2mat(arrayfun(@(x) sum(~isnan(lecTc{i}.avgAmpNorm(x,:))),1:size(lecTc{i}.avgAmpNorm,1),'UniformOutput',false));
    lecTc{i}.mAvgAmpNormSem = lecTc{i}.mAvgAmpNormStd ./ sqrt(lecTc{i}.mAvgAmpNormN)';
    
    lecTc{i}.mAvgAmpNormTime = mean(lecTc{i}.normTime(1:4:end,:),2,'omitnan');
    lecTc{i}.mAvgAmpNormTime(6) = 0;

    cd(path);
end

cd ../;

save(['lec_merge_tc'], 'lecTc'); 