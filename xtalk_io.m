close all; clear all;

path = ['/Volumes/KO_Portable/MEC_LEC_CA3/data_converted/']; cd(path);
path2 = ['/Volumes/KO_Portable/MEC_LEC_CA3/analysis/VC/XTALK'];

parent = dir(pwd); parent = parent([parent(:).isdir]); 
parent = parent(~ismember({parent(:).name},{'.','..'}));
xtalk = [];

for i = 1:length(parent)
    cd(parent(i).name) %go to indv animal folder
    child = dir(pwd); child = child([child(:).isdir]); 
    child = child(~ismember({child(:).name},{'.','..'})); %get list of cell folders
       
    for ii = 1:length(child)
        cd(child(ii).name); %navigate to an individual cell
        pathway = readcell('pathway.txt'); 
        if pathway{string(pathway)=='include',2} == 1 && string(pathway{string(pathway)=='blue',2}) == "0"
           xtalk =[xtalk; string(child(ii).name)];
           a = dir('K*ioStruct.mat');
           copyfile(a.name, path2);
           cd ../

        else
           cd ../
           continue
        end       
    end
    cd ../
end

cd(path2);

files = dir('K*ioStruct.mat');
xStruct.blue = [];
xStruct.id = [];
xStruct.path = [];
xStruct.red = [];
xStruct.stim = [];
xStruct.ttxId  = [];
xStruct.pathTtx = [];
xStruct.blueTtx = [];
xStruct.redTtx = [];
xStruct.ttxStim = [];

for i = 1:length(files)
    load(files(i).name)
    xStruct.blue = [xStruct.blue; ioStruct.dataBlueAvg{2}(:,1)];
    xStruct.id = [xStruct.id; string(strcat(ioStruct.Animal, ioStruct.Cell,'_m_', extractAfter(ioStruct.dataAvgSweep{2}(:,1), '_1_')))];
    xStruct.path = [xStruct.path; string(strcat(ioStruct.Animal, ioStruct.Cell))];
    xStruct.stim = [xStruct.stim; ioStruct.dataAvgStim{2}(:,1)];

    xStruct.red = [xStruct.red; ioStruct.dataRedAvg{2}(:,1)];
    %xStruct.redId = [xStruct.redId; extractAfter(ioStruct.dataAvgSweep{2}(:,1), '_1_')];
    %xStruct.redPath = [xStruct.redPath; string(strcat(ioStruct.Animal, ioStruct.Cell))];
    if sum(contains(fieldnames(ioStruct), 'ttx', 'IgnoreCase', 1)) > 0
        xStruct.blueTtx = [xStruct.blueTtx; ioStruct.dataBlueTtx{2}(:,1)];
        xStruct.ttxId = [xStruct.ttxId; string(strcat(ioStruct.Animal, ioStruct.Cell,'_m_', extractAfter(ioStruct.dataAvgSweepTtx{2}(:,1), '_1_')))];
        xStruct.pathTtx = [xStruct.pathTtx; string(strcat(ioStruct.Animal, ioStruct.Cell))];
        xStruct.ttxStim = [xStruct.ttxStim; ioStruct.dataAvgStimTtx{2}(:,1)];

        xStruct.redTtx = [xStruct.redTtx; ioStruct.dataRedTtx{2}(:,1)];

    end
end

[u, ~, g] = unique(extractAfter(xStruct.id, '_m_'),'stable');

xStruct.blueAvg = splitapply(@(x) mean(x, 'omitnan'), xStruct.blue, g);
xStruct.blueAvgStd = splitapply(@(x) std(x, 'omitnan'), xStruct.blue, g);
xStruct.blueAvgSem = xStruct.blueAvgStd ./ length(xStruct.path);
xStruct.avgId = g;

xStruct.redAvg = splitapply(@(x) mean(x, 'omitnan'), xStruct.red, g);
xStruct.redAvgStd = splitapply(@(x) std(x, 'omitnan'), xStruct.red, g);
xStruct.redAvgSem = xStruct.redAvgStd ./ length(xStruct.path);

[u, ~, g] = unique(extractAfter(xStruct.ttxId, '_m_'),'stable');

xStruct.blueAvgTtx =  splitapply(@(x) mean(x, 'omitnan'), xStruct.blueTtx, g);
xStruct.blueAvgTtxStd = splitapply(@(x) std(x, 'omitnan'), xStruct.blueTtx, g);
xStruct.blueAvgTtxSem = xStruct.blueAvgTtxStd ./ length(xStruct.path);

xStruct.redAvgTtx =  splitapply(@(x) mean(x, 'omitnan'), xStruct.redTtx, g);
xStruct.redAvgTtxStd = splitapply(@(x) std(x, 'omitnan'), xStruct.redTtx, g);
xStruct.redAvgTtxSem = xStruct.redAvgTtxStd ./ length(xStruct.path);

save(['xTalk_io'], 'xStruct');
