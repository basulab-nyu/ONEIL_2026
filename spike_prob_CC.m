%% Spike Probability
% Find spike prob of cells that spike only, need to first run
% Compile_Kinetics_Draft and Compile_Kinetics_Binned

clear all; close all;

path = ['/Volumes/KO_Portable/MEC_LEC_CA3/analysis/CC/IO']; cd(path);
load('Compile_Kinetics_Binned.mat');

path2 = ['/Volumes/KO_Portable/MEC_LEC_CA3/analysis/mW']; cd(path2);
load('led_power.mat');


spikeProb.lecAll = [spikeProb.lecBDual(:,1); spikeProb.lecRDual(:,1)];
spikeProb.lecAll(:,3) = [spikeProb.lecBDual(:,4); spikeProb.lecRDual(:,3)];
spikeProb.lecAll(:,2) = [spikeProb.lecBDual(:,2); spikeProb.lecRDual(:,2)];
spikeProb.lecAllPathCell = [spikeProb.lecBDualPathCell; spikeProb.lecRDualPathCell];
spikeProb.lecAllStim = [spikeProb.lecBDualStim; spikeProb.lecRDualStim];


spikeProb.mecAll = [spikeProb.mecRDual(:,1); spikeProb.mecBDual(:,1)];
spikeProb.mecAll(:,3) = [spikeProb.mecRDual(:,3); spikeProb.mecBDual(:,4)];
spikeProb.mecAll(:,2) = [spikeProb.mecRDual(:,2); spikeProb.mecBDual(:,2)];
spikeProb.mecAllPathCell = [spikeProb.mecRDualPathCell; spikeProb.mecBDualPathCell];
spikeProb.mecAllStim = [spikeProb.mecRDualStim; spikeProb.mecBDualStim];


fields = fieldnames(spikeProb);
f = fields(endsWith(fieldnames(spikeProb),'PathCell'));

% 
for i = 1:length(f)
    tempF =  erase(f{i}, 'PathCell');
    tempF2 = strcat(tempF, 'Stim');
    tempNum = strcat(tempF, 'Num');
    tempCell =  strcat(tempNum, 'PathCell');
    [~, col] = size(spikeProb.(tempF));
    if isfield(spikeProb, tempF)
        for ii = 1:col
            [u, ~, g] = unique(extractAfter(spikeProb.(f{i})(:,1), 'K'), 'stable');
            spikeCell.(tempF)(:,ii) = splitapply(@(x) sum(x), spikeProb.(tempF)(:,ii), g); %store summed data to find who spikes and who doesn't
            spikeCell.(f{i})(:,ii) = u; %all the individual cell names
            tempIdx = spikeCell.(tempF)(:,ii) > 0,1; %indices of cells that spike 
            spikeCell.(tempNum)(1:sum(tempIdx),ii) = spikeCell.(tempF)(tempIdx); 
            spikeCell.(tempCell)(1:sum(tempIdx),ii) = u(tempIdx);
            tempIdx = contains(spikeProb.(f{i})(:,1), u(tempIdx));
            spikeCellData.(tempF)(1:sum(tempIdx), ii) = spikeProb.(tempF)(tempIdx,ii);
            spikeCellData.(f{i})(1:sum(tempIdx), (ii*2)-1:ii*2) = [spikeProb.(f{i})(tempIdx,1), spikeProb.(tempF2)(tempIdx,2)];
            
            % spikeCellData.(tempF)(1:sum(contains(spikeProb.(f{i})(:,1),spikeCell.(tempCell)(:,ii))),ii) = ... 
            %     [spikeProb.(tempF)(contains(spikeProb.(f{i})(:,1),spikeCell.(tempCell)(:,ii)),ii) ];
            % spikeCellData.(f{i})(1:sum(contains(spikeProb.(f{i})(:,1), spikeCell.(tempCell)(:,ii))),((ii*2)-1):(ii*2)) = ...
            %     [spikeProb.(f{i})(contains(spikeProb.(f{i})(:,1), spikeCell.(tempCell)(:,ii)),1), spikeProb.(tempF2)(contains(spikeProb.(f{i})(:,1), spikeCell.(tempCell)(:,ii)),2)];
        end
    end
end
% 
% for i = 1:length(f)
%     tempF =  erase(f{i}, 'PathCell');
%     tempF3 = strcat(tempF, 'Stim');
%     tempNum = strcat(tempF, 'Num');
%     tempCell =  strcat(tempNum, 'PathCell');
%     if isfield(spikeProb, tempF)
%         [u, ~, g] = unique(extractAfter(spikeProb.(f{i})(:,1), 'K'), 'stable');
%         spikeCell.(tempF) = splitapply(@(x) sum(x), spikeProb.(tempF)(:,1), g);
%         spikeCell.(f{i}) = u;
%         spikeCell.(tempNum) = spikeCell.(tempF)(spikeCell.(tempF)(:,1) > 0,1);
%         spikeCell.(tempCell) = u(spikeCell.(tempF)(:,1) > 0);
%         spikeCellData.(tempF) = [spikeProb.(tempF)(contains(spikeProb.(f{i})(:,1),spikeCell.(tempCell)),1) ];
%         spikeCellData.(f{i}) = [spikeProb.(f{i})(contains(spikeProb.(f{i})(:,1), spikeCell.(tempCell)),1),spikeProb.(tempF3)(contains(spikeProb.(f{i})(:,1), spikeCell.(tempCell)),2)];
%     end
% end

fields = fieldnames(spikeCell);
f = fields(endsWith(fieldnames(spikeCell),'PathCell'));

% for i = 1:length(f)
%     [~, col] = size(spikeCell.(f{i}));
%     for ii = 1:col
% 
%     end
% 
% end

%Sort the all vectors so they are ordered the same way
[~, tempIdx] = ismember(spikeCell.lecAllPathCell(:,1), spikeCell.comboPathCell(:,1));
spikeCell.comboPathCell = spikeCell.comboPathCell(tempIdx,:);
spikeCell.combo = spikeCell.combo(tempIdx,:);

spikeCell.all = [spikeCell.lecAll(:,1), spikeCell.mecAll(:,1), spikeCell.combo(:,1)];
spikeCellAvg.all = sum(spikeCell.all,2);
spikeCell.allPathCell = [spikeCell.lecAllPathCell, spikeCell.mecAllPathCell, spikeCell.comboPathCell];

spikeCellAvg.allNum = length(spikeCellAvg.all(spikeCellAvg.all(:,1) > 0));
spikeCellAvg.allNumPathCell = spikeCell.allPathCell(spikeCellAvg.all(:,1) > 0);

spikeCell.dual = [spikeCell.lecAll(:,1), spikeCell.mecAll(:,1)];
spikeCellAvg.dual = sum([spikeCell.lecAll(:,1), spikeCell.mecAll(:,1)],2);
spikeCellAvg.dualNum = length(spikeCellAvg.dual(spikeCellAvg.dual(:,1) > 0));
spikeCellAvg.dualPathCell = spikeCell.mecAllPathCell(spikeCellAvg.dual(:,1) > 0);

spikeCellAvg.mec = length(spikeCell.mecAllNum);
spikeCellAvg.lec = length(spikeCell.lecAllNum);

%Change to doing it with spike prob instead of spike bin 
ecOnly = 0; ecDual = 0;
lecOnly = 0; lecDualOnly = 0;
none = 0; allIn = 0;
mecOnly = 0; mecDualOnly = 0;
comboLec = 0; comboMec = 0; comboOnly = 0;

for i = 1:length(spikeCell.all)
    if spikeCell.all(i,3) == 0 & spikeCell.all(i,2) == 0 & spikeCell.all(i,1) == 0
        none = none + 1;
    elseif spikeCell.all(i,3) == 0 & spikeCell.all(i,2) == 0 & spikeCell.all(i,1) > 0
        lecOnly = lecOnly + 1;
    elseif spikeCell.all(i,3) == 0 & spikeCell.all(i,2) > 0 & spikeCell.all(i,1) > 0
        ecOnly = ecOnly + 1;
    elseif spikeCell.all(i,3) > 0 & spikeCell.all(i,2) == 0 & spikeCell.all(i,1) > 0
        comboLec = comboLec + 1;
    elseif spikeCell.all(i,3) > 0 & spikeCell.all(i,2) > 0 & spikeCell.all(i,1) > 0
        allIn = allIn + 1;
    elseif spikeCell.all(i,3) == 0 & spikeCell.all(i,2) > 0 & spikeCell.all(i,1) == 0
        mecOnly = mecOnly + 1;
    elseif spikeCell.all(i,3) > 0 & spikeCell.all(i,2) > 0 & spikeCell.all(i,1) == 0
        comboMec = comboMec + 1;
    elseif spikeCell.all(i,3) > 0 & spikeCell.all(i,2) == 0 & spikeCell.all(i,1) == 0
        comboOnly = comboOnly + 1;
    end
    
    if spikeCell.all(i,2) > 0 & spikeCell.all(i,1) == 0
        mecDualOnly = mecDualOnly + 1;
    elseif spikeCell.all(i,2) > 0 & spikeCell.all(i,1) > 0
        ecDual = ecDual + 1;
    elseif spikeCell.all(i,2) == 0 & spikeCell.all(i,1) > 0
        lecDualOnly = lecDualOnly + 1;
    end
end

spikeCellAvg.ecOnly = ecOnly; 
spikeCellAvg.ecDual = ecDual;
spikeCellAvg.lecOnly = lecOnly; 
spikeCellAvg.lecDualOnly = lecDualOnly;

spikeCellAvg.none = none; 
spikeCellAvg.allIn = allIn;
spikeCellAvg.mecOnly = mecOnly; 
spikeCellAvg.mecDualOnly = mecDualOnly;
spikeCellAvg.comboLec = comboLec; 
spikeCellAvg.comboMec = comboMec; 
spikeCellAvg.comboOnly = comboOnly;

f = ["mecAll", "lecAll","combo"];

for i = 1:length(f)
    tempF = strcat(f(i), 'PathCell'); 
    tempF2 = strcat(f(i), 'Stim');
    tempF3 = strcat(f(i), 'Std');
    tempF4 = strcat(f(i), 'Num');
    tempF5 = strcat(f(i), 'Sem');
    spikeOnly.(f(i)) = [];
    spikeOnly.(tempF) = [];
    spikeOnly.(tempF2) = [];
    %spikeOnly.(tempF6) = [];

    for ii = 1:length(spikeCellAvg.allNumPathCell)
        tempIdx = contains(spikeProb.(tempF)(:,1), spikeCellAvg.allNumPathCell(ii));
        spikeOnly.(f(i)) = [spikeOnly.(f(i)); spikeProb.(f(i))(tempIdx,:)];
        spikeOnly.(tempF) = [spikeOnly.(tempF); spikeProb.(tempF)(tempIdx,:), extractBefore(spikeProb.(tempF)(tempIdx,end),'K')];
        spikeOnly.(tempF2) = [spikeOnly.(tempF2); spikeProb.(tempF2)(tempIdx,2)];
        spikeOnly.(tempF)(ismissing(spikeOnly.(tempF)(:,3)),3) = "0";
    end

    tempIdx = find(spikeOnly.(tempF)(:,end) ~= "0");
    spikeOnly.(f(i)) = spikeOnly.(f(i))(tempIdx,:);
    spikeOnly.(tempF2) = spikeProb.(tempF2)(tempIdx,2);
    spikeOnly.(tempF) = spikeOnly.(tempF)(tempIdx,:);

    [u,~,g] = unique(spikeOnly.(tempF)(:,3), 'stable');
    spikeOnlyBin.(f(i)) = splitapply(@(x) mean(x, 1,'omitnan'),spikeOnly.(f(i)), g);
    spikeOnlyBin.(tempF3) = splitapply(@(x) std(x,'omitnan'),spikeOnly.(f(i)), g);
    spikeOnlyBin.(tempF4) = histc(g, unique(g, 'stable'));
    spikeOnlyBin.(tempF5) = spikeOnlyBin.(tempF3)./sqrt(spikeOnlyBin.(tempF4));
    spikeOnlyBin.(tempF) = u;


    [u,~,g] = unique(spikeOnly.(tempF2), 'stable');
    spikeOnlyAvg.(f(i)) = splitapply(@(x) mean(x, 1,'omitnan'),spikeOnly.(f(i)), g);
    spikeOnlyAvg.(tempF3) = splitapply(@(x) std(x,'omitnan'),spikeOnly.(f(i)), g);
    spikeOnlyAvg.(tempF4) = histc(g, unique(g, 'stable'));
    spikeOnlyAvg.(tempF5) = spikeOnlyAvg.(tempF3)./sqrt(spikeOnlyAvg.(tempF4));
    spikeOnlyAvg.(tempF) = u;

end



cd(path);
save(['Compile_Kinetics_CellSpikeProb.mat'], 'spikeCell', 'spikeCellAvg', 'spikeProb', 'spikeProbAvg', 'spikeBin','spikeBinAvg', 'spikeCellData', 'spikeOnly', 'spikeOnlyBin','spikeOnlyAvg');



