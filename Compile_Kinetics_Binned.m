%% Bin I/O for CC or VC
% Need to run Compile_Kinetics_draft first
% you will also get Spike prob binned and ei binned

clear all; close all;

%What do you want to plot?

%need to NaN really small values in CC 


plt = 0; %Do you want to plot all light intensities (ie mono and polysyaptic for all kinetics files)
plt2 = 0; %Do you want to plot all the ratios

path = '/Volumes/KO_Portable/MEC_LEC_CA3/analysis/VC/IO'; cd(path);
load('Compile_Kinetics.mat');


path2 = ['/Volumes/KO_Portable/MEC_LEC_CA3/analysis/mW']; cd(path2);
load('led_power.mat');

addpath('/Volumes/KO_Portable/MATLAB/EphysAnalysis_2025');

fields = fieldnames(io);
f = fields(endsWith(fieldnames(io),'Stim'));

for i = 1:length(f)
    if contains(f{i}, 'B')
        temp = cellfun(@(y) arrayfun(@(x) find(x==blue(:,1)), y), io.(f{i}), 'UniformOutput', false);
        temp2 = cellfun(@(x) blue(x,4), temp, 'UniformOutput', false);
        io.(f{i}) = cellfun(@(x, y, z) [x, y, z], io.(f{i}), temp, temp2, 'UniformOutput', false);
        ppr.(f{i}) = io.(f{i});
        if contains(path, 'CC','IgnoreCase',true)
            spikeProb.(f{i}) = io.(f{i}){1};
        end
    elseif contains(f{i}, 'R')
        temp = cellfun(@(y) arrayfun(@(x) find(x==red(:,1)), y), io.(f{i}), 'UniformOutput', false);
        temp2 = cellfun(@(x) red(x,4), temp, 'UniformOutput', false);
        io.(f{i}) = cellfun(@(x, y, z) [x, y, z], io.(f{i}), temp, temp2, 'UniformOutput', false);
        ppr.(f{i}) = io.(f{i});
        if contains(path, 'CC','IgnoreCase',true)
            spikeProb.(f{i}) = io.(f{i}){1};
        end
    elseif contains(f{i}, 'combo')
        temp = cellfun(@(y) arrayfun(@(x) find(x==combo(:,1)), y), io.(f{i}), 'UniformOutput', false);
        temp2 = cellfun(@(x) combo(x,4), temp, 'UniformOutput', false);
        io.(f{i}) = cellfun(@(x, y, z) [x, y, z], io.(f{i}), temp, temp2, 'UniformOutput', false);
        ppr.(f{i}) = io.(f{i});
        if contains(path, 'CC','IgnoreCase',true)
            spikeProb.(f{i}) = io.(f{i}){1};
        end
    end
end


%Bin to get average for individual cells
f = fields(contains(fields, 'PathCell'));

for i = 1:length(f)
    tempF = string(replace(f(i), 'PathCell', 'Stim'));
    if sum(cellfun(@(x) isempty(x), io.(tempF)) > 0)
        continue
    end
    tempC = cellfun(@(x, y) strcat(string(x(:,end)), y), io.(tempF), io.(f{i}), 'UniformOutput', false);
    tempC = cellfun(@(x) fillmissing(x, 'constant', "0"), tempC, 'UniformOutput', false);
    io.(f{i}) = cellfun(@(x, y) [x, y], io.(f{i}), tempC, 'UniformOutput', false);
    tempC = cellfun(@(x, y) strcat(string(x(:,2)), y(:,1)), io.(tempF), io.(f{i}), 'UniformOutput', false);
    ei.(f{i}) = tempC;
    ppr.(f{i}) = tempC;
    if contains(path, 'CC','IgnoreCase',true)
        spikeProb.(f{i}) = io.(f{i}){1};
    end
    % if sum(cellfun(@(x) isempty(x), io.(tempF)) > 0)
    %     continue
    % end
    if contains(path, 'VC','IgnoreCase',true)
        tempF = erase(tempF, 'Stim');
        tempCell = find(contains(kinetics, {'amp_io10', 'sum_io10'})); 
        for ii = 1:length(tempCell)
            tempIdx = find(io.(tempF){tempCell(ii)}(:,1) <= 0);
            tempCell2 = find(contains(kinetics, '10'));
            for iii = 1:length(tempCell2)
                io.(tempF){tempCell2(iii)}(tempIdx) = NaN;
            end
        end
        tempCell = find(contains(kinetics, {'amp_io80', 'sum_io80'})); 
        for ii = 1:length(tempCell)
            tempIdx = io.(tempF){tempCell(ii)}(:,1) >= 0;
            tempCell2 = find(contains(kinetics, '80'));
            for iii = 1:length(tempCell2)
                io.(tempF){tempCell2(iii)}(tempIdx) = NaN;
            end
        end
    end

    if contains(path, 'CC','IgnoreCase',true)
        tempF = erase(tempF, 'Stim');
        tempCell = find(contains(kinetics, {'amp_io70', 'sum_io70'})); 
        for ii = 1:length(tempCell)
            tempIdx = io.(tempF){tempCell(ii)}(:,1) <= 0;
            tempCell2 = find(contains(kinetics, '70'));
            for iii = 1:length(tempCell2)
                io.(tempF){tempCell2(iii)}(tempIdx) = NaN;
            end
        end
    end
    
    if contains(path, 'VC')
        tempF = erase(tempF, 'Stim');
        tempF2 = strcat(tempF, 'Stim');
        [idxV, idx] = ismember(ei.(f{i}){1}, ei.(f{i}){2});
        ei.(tempF) = abs(io.(tempF){2}(idx(idxV),:) ./io.(tempF){1}(idxV,:)); 
        if all(matches(ei.(f{i}){1}(idxV), ei.(f{i}){2}(idx(idxV))))
            ei.(f{i}) = [ei.(f{i}){1}(idxV,:), io.(f{i}){2}(idx(idxV),end)];
            ei.(tempF2) = io.(tempF2){1}(idxV,:);
        else
             fprintf([f{i} 'something is wrong w/ ei'])
        end
        
        [u, ~, g] = unique(ei.(f{i})(:,1), 'stable');
        eiCell.(f{i}) = u;
        eiCell.(tempF) = splitapply(@(x) mean(x, 1, 'omitnan'), ei.(tempF),g); 
        eiCell.(tempF2) = splitapply(@(x) mean(x, 1, 'omitnan'), ei.(tempF2),g);
        
        [u, ~, g] = unique(ei.(tempF2)(:,2), 'stable');
        eiAvg.(tempF) = splitapply(@(x) mean(x, 1, 'omitnan'), ei.(tempF),g);
        eiAvg.(tempF2) = u; 
        eiAvg.(f{i}) = unique(extractAfter(ei.(f{i})(:,1),'K'));
        
        tempSd = strcat(tempF, 'Sd');
        tempSem = strcat(tempF, 'Sem');
        tempN = strcat(tempF, 'Num');
        
        eiAvg.(tempSd) = splitapply(@(x) std(x, 1, 'omitnan'), ei.(tempF),g);
        eiAvg.(tempN) = histc(g, unique(g));
        eiAvg.(tempSem) = eiAvg.(tempSd) ./ sqrt(eiAvg.(tempN));

        [u, ~, g] = unique(ei.(f{i})(:,end), 'stable');
        eiBin.(tempF) = splitapply(@(x) mean(x, 1, 'omitnan'), ei.(tempF),g);
        eiBin.(tempF) = eiBin.(tempF)((u~= "0"),:);
        eiBin.(f{i}) = [u(u~= "0"), extractBefore(u(u~= "0"), 'K')];
    end

    tempF = erase(tempF, 'Stim');
    [u, ~, g] = cellfun(@(x) unique(x(:,end), 'stable'), io.(f{i}), 'UniformOutput', false);
    ioBin.(tempF) = cellfun(@(data, groups) splitapply(@(x) mean(x,1, 'omitnan'), data, groups), io.(tempF), g, 'UniformOutput', false);
    ioBin.(tempF) = cellfun(@(x, y) y(~(x == "0"),:), u, ioBin.(tempF), 'UniformOutput', false); 
    ioBin.(f{i}) = cellfun(@(x) x(~(x== "0")), u, 'UniformOutput', false);
    ioBin.(f{i}) = cellfun(@(x) [x, extractBefore(x, 'K')], ioBin.(f{i}), 'UniformOutput', false);

    tempI = find(cellfun(@(x) isempty(x), ppr.(tempF)) == 0); 
    pprBin.(tempF) =  cellfun(@(data, groups) splitapply(@(x) mean(x,1, 'omitnan'), data, groups), ppr.(tempF)(tempI), g(tempI), 'UniformOutput', false);
    pprBin.(tempF) = cellfun(@(x, y) y(~(x == "0"),:), u(tempI), pprBin.(tempF), 'UniformOutput', false);
    pprBin.(f{i}) = cellfun(@(x) x(~(x== "0")), u(tempI), 'UniformOutput', false);
    pprBin.(f{i}) = cellfun(@(x) [x, extractBefore(x, 'K')], pprBin.(f{i}), 'UniformOutput', false);

    if contains(path, 'CC')
        [u, ~, g] = unique(spikeProb.(f{i})(:,end),'stable');
        spikeBin.(tempF) = splitapply(@(x) mean(x,1, 'omitnan'), spikeProb.(tempF), g);
        spikeBin.(tempF) = spikeBin.(tempF)(~(u == "0"),:); 
        spikeBin.(f{i}) = u(~(u== "0"));
        spikeBin.(f{i}) = [spikeBin.(f{i}), extractBefore(spikeBin.(f{i}), 'K')];
    end

end


% Make field with lec/mec (individually labelled) and lec/mec all 
% (Cell Averages)

if contains(path, 'VC','IgnoreCase',true)
    ioBin.lec = cellfun(@(x,y) [x(:,1);y(:,1)], ioBin.lecB, ioBin.lecR,'UniformOutput', false);
    ioBin.lecPathCell = cellfun(@(x,y) [x;y], ioBin.lecBPathCell, ioBin.lecRPathCell,'UniformOutput', false);
    
    ioBin.mec = cellfun(@(x,y) [x(:,1);y(:,1)], ioBin.mecB, ioBin.mecR,'UniformOutput', false);
    ioBin.mecPathCell = cellfun(@(x,y) [x;y], ioBin.mecBPathCell, ioBin.mecRPathCell,'UniformOutput', false);
    
    ioBin.lecDual = cellfun(@(x,y) [x(:,1);y(:,1)], ioBin.lecBDual, ioBin.lecRDual,'UniformOutput', false);
    ioBin.lecDualPathCell = cellfun(@(x,y) [x;y], ioBin.lecBDualPathCell, ioBin.lecRDualPathCell,'UniformOutput', false);
    
    ioBin.mecDual = cellfun(@(x,y) [x(:,1);y(:,1)], ioBin.mecRDual, ioBin.mecBDual,'UniformOutput', false);
    ioBin.mecDualPathCell = cellfun(@(x,y) [x;y], ioBin.mecRDualPathCell, ioBin.mecBDualPathCell,'UniformOutput', false);
   
    ioBin.lecAll = cellfun(@(x,y,z,a) [x(:,1);y(:,1);z(:,1);a(:,1)], ioBin.lecB, ioBin.lecR, ioBin.lecBDual, ioBin.lecRDual, 'UniformOutput', false);
    ioBin.lecAllPathCell = cellfun(@(x,y,z,a) [x;y;z;a], ioBin.lecBPathCell, ioBin.lecRPathCell, ioBin.lecBDualPathCell, ioBin.lecRDualPathCell, 'UniformOutput', false);
    
    ioBin.mecAll = cellfun(@(x,y,z,a) [x(:,1);y(:,1);z(:,1);a(:,1)], ioBin.mecB, ioBin.mecR, ioBin.mecRDual, ioBin.mecBDual, 'UniformOutput', false);
    ioBin.mecAllPathCell = cellfun(@(x,y,z,a) [x;y;z;a], ioBin.mecBPathCell, ioBin.mecRPathCell, ioBin.mecRDualPathCell, ioBin.mecBDualPathCell, 'UniformOutput', false);
    %pprBin.lec = cellfun(@(x,y) [x(:,1);y(:,1)], pprBin.lecB, pprBin.lecR,'UniformOutput', false);
    
    % pprBin.lecPathCell = cellfun(@(x,y) [x;y], pprBin.lecBPathCell, pprBin.lecRPathCell,'UniformOutput', false);
    % 
    % pprBin.mec = cellfun(@(x,y) [x;y], pprBin.mecB, pprBin.mecR,'UniformOutput', false);
    % pprBin.mecPathCell = cellfun(@(x,y) [x;y], pprBin.mecBPathCell, pprBin.mecRPathCell,'UniformOutput', false);
    % 
    % pprBin.lecDual = cellfun(@(x,y) [x(:,1:2);y(:,1:2)], pprBin.lecBDual, pprBin.lecRDual,'UniformOutput', false);
    % pprBin.lecDualPathCell = cellfun(@(x,y) [x;y], pprBin.lecBDualPathCell, pprBin.lecRDualPathCell,'UniformOutput', false);
    % 
    % pprBin.mecDual = cellfun(@(x,y) [x(:,1);y(:,1)], pprBin.mecRDual, pprBin.mecBDual,'UniformOutput', false);
    % pprBin.mecDualPathCell = cellfun(@(x,y) [x;y], pprBin.mecRDualPathCell, pprBin.mecBDualPathCell,'UniformOutput', false);
    
    % pprBin.lecAll = cellfun(@(x,y,z,a) [x(:,2);y(:,2);z(:,2);a(:,2)], pprBin.lecB, pprBin.lecR, pprBin.lecBDual, pprBin.lecRDual, 'UniformOutput', false);
    % pprBin.lecAllPathCell = cellfun(@(x,y,z,a) [x;y;z;a], pprBin.lecBPathCell, pprBin.lecRPathCell, pprBin.lecBDualPathCell, pprBin.lecRDualPathCell, 'UniformOutput', false);
    % 
    % pprBin.mecAll = cellfun(@(x,y,z,a) [x(:,1);y(:,1);z(:,1);a(:,1)], pprBin.mecB, pprBin.mecR, pprBin.mecRDual, pprBin.mecBDual, 'UniformOutput', false);
    % pprBin.mecAllPathCell = cellfun(@(x,y,z,a) [x;y;z;a], pprBin.mecBPathCell, pprBin.mecRPathCell, pprBin.mecRDualPathCell, pprBin.mecBDualPathCell, 'UniformOutput', false);

    t = [1, 2, 27, 28];
    for n = 1:length(t)
        ioBin.lecLec{n} = [ioBin.lecB{t(n)}(:,end); ioBin.lecR{t(n)}(:,end); ioBin.lecBDual{t(n)}(:,2); ioBin.lecRDual{t(n)}(:,2)];
        ioBin.lecMec{n} = [ioBin.lecBDual{t(n)}(:,4); ioBin.lecRDual{t(n)}(:,3)];
        ioBin.lecLecPathCell{n} = [ioBin.lecBPathCell{t(n)}; ioBin.lecRPathCell{t(n)}; ioBin.lecBDualPathCell{t(n)}; ioBin.lecRDualPathCell{t(n)}];
        ioBin.lecMecPathCell{n} = [ioBin.lecBDualPathCell{t(n)}; ioBin.lecRDualPathCell{t(n)}];
    
        ioBin.mecMec{n} = [ioBin.mecB{t(n)}(:,end); ioBin.mecR{t(n)}(:,end); ioBin.mecRDual{t(n)}(:,2); ioBin.mecBDual{t(n)}(:,2)];
        ioBin.mecLec{n} = [ioBin.mecRDual{t(n)}(:,3); ioBin.mecBDual{t(n)}(:,4)];
        ioBin.mecMecPathCell{n} = [ioBin.mecBPathCell{t(n)}; ioBin.mecRPathCell{t(n)}; ioBin.mecRDualPathCell{t(n)}; ioBin.mecBDualPathCell{t(n)}];
        ioBin.mecLecPathCell{n} = [ioBin.mecRDualPathCell{t(n)}; ioBin.mecBDualPathCell{t(n)}];

        pprBin.lecLec{n} = [pprBin.lecB{t(n)}(:,1); pprBin.lecR{t(n)}(:,1); pprBin.lecBDual{t(n)}(:,1); pprBin.lecRDual{t(n)}(:,1)];
        pprBin.lecMec{n} = [pprBin.lecBDual{t(n)}(:,2); pprBin.lecRDual{t(n)}(:,2)];
        pprBin.lecLecPathCell{n} = [pprBin.lecBPathCell{t(n)}; pprBin.lecRPathCell{t(n)}; pprBin.lecBDualPathCell{t(n)}; pprBin.lecRDualPathCell{t(n)}];
        pprBin.lecMecPathCell{n} = [pprBin.lecBDualPathCell{t(n)}; pprBin.lecRDualPathCell{t(n)}];
        
        pprBin.mecMec{n} = [pprBin.mecB{t(n)}(:,1); pprBin.mecR{t(n)}(:,1); pprBin.mecRDual{t(n)}(:,1); pprBin.mecBDual{t(n)}(:,1)];
        pprBin.mecLec{n} = [pprBin.mecRDual{t(n)}(:,2); pprBin.mecBDual{t(n)}(:,2)];
        pprBin.mecMecPathCell{n} = [pprBin.mecBPathCell{t(n)}; pprBin.mecRPathCell{t(n)}; pprBin.mecRDualPathCell{t(n)}; pprBin.mecBDualPathCell{t(n)}];
        pprBin.mecLecPathCell{n} = [pprBin.mecRDualPathCell{t(n)}; pprBin.mecBDualPathCell{t(n)}];
    end
    eiBin.mecAll = [eiBin.mecB(:,1); eiBin.mecR(:,1); eiBin.mecRDual(:,1); eiBin.mecBDual(:,1)];
    eiBin.mecAllPathCell = [eiBin.mecBPathCell; eiBin.mecRPathCell; eiBin.mecRDualPathCell; eiBin.mecBDualPathCell];

    eiBin.lecAll = [eiBin.lecB(:,1); eiBin.lecR(:,1); eiBin.lecRDual(:,1); eiBin.lecBDual(:,1)];
    eiBin.lecAllPathCell = [eiBin.lecBPathCell; eiBin.lecRPathCell; eiBin.lecRDualPathCell; eiBin.lecBDualPathCell];

    eiBin.mecDual = [eiBin.mecRDual(:,1); eiBin.mecBDual(:,1)];
    eiBin.mecDualPathCell = [eiBin.mecRDualPathCell; eiBin.mecBDualPathCell];

    eiBin.lecDual = [eiBin.lecBDual(:,1); eiBin.lecRDual(:,1)];
    eiBin.lecDualPathCell = [eiBin.lecBDualPathCell; eiBin.lecRDualPathCell];

    eiBin.mec = [eiBin.mecB(:,1); eiBin.mecR(:,1)];
    eiBin.mecPathCell = [eiBin.mecBPathCell; eiBin.mecRPathCell];

    eiBin.lec = [eiBin.lecB(:,1); eiBin.lecR(:,1)];
    eiBin.lecPathCell = [eiBin.lecBPathCell; eiBin.lecRPathCell];
    

elseif contains(path, 'CC','IgnoreCase',true)
    ioBin.lecAll = cellfun(@(x,y) [x(:,1);y(:,1)], ioBin.lecBDual, ioBin.lecRDual, 'UniformOutput', false);
    ioBin.lecAllPathCell = cellfun(@(x,y) [x;y], ioBin.lecBDualPathCell, ioBin.lecRDualPathCell, 'UniformOutput', false);
    
    ioBin.mecAll = cellfun(@(x,y) [x(:,1);y(:,1)], ioBin.mecRDual, ioBin.mecBDual, 'UniformOutput', false);
    ioBin.mecAllPathCell = cellfun(@(x,y) [x;y], ioBin.mecRDualPathCell, ioBin.mecBDualPathCell, 'UniformOutput', false);
    
    % pprBin.lecAll = cellfun(@(x,y) [x(:,1);y(:,1)], pprBin.lecBDual, pprBin.lecRDual, 'UniformOutput', false);
    % pprBin.lecAllPathCell = cellfun(@(x,y) [x;y], pprBin.lecBDualPathCell, pprBin.lecRDualPathCell, 'UniformOutput', false);
    % 
    % pprBin.mecAll = cellfun(@(x,y) [x(:,1);y(:,1)], pprBin.mecRDual, pprBin.mecBDual, 'UniformOutput', false);
    % pprBin.mecAllPathCell = cellfun(@(x,y) [x;y], pprBin.mecRDualPathCell, pprBin.mecBDualPathCell, 'UniformOutput', false);
    % 
    t = [1, 14];
    for n = 1:length(t)
        ioBin.lecLec{n} = [ioBin.lecBDual{t(n)}(:,2); ioBin.lecRDual{t(n)}(:,2)];
        ioBin.lecMec{n} = [ioBin.lecBDual{t(n)}(:,4); ioBin.lecRDual{t(n)}(:,3)];
        ioBin.lecLecPathCell{n} = [ioBin.lecBDualPathCell{t(n)}; ioBin.lecRDualPathCell{t(n)}];
        ioBin.lecMecPathCell{n} = [ioBin.lecBDualPathCell{t(n)}; ioBin.lecRDualPathCell{t(n)}];
    
        ioBin.mecMec{n} = [ioBin.mecRDual{t(n)}(:,2); ioBin.mecBDual{t(n)}(:,2)];
        ioBin.mecLec{n} = [ioBin.mecRDual{t(n)}(:,3); ioBin.mecBDual{t(n)}(:,4)];
        ioBin.mecMecPathCell{n} = [ioBin.mecRDualPathCell{t(n)}; ioBin.mecBDualPathCell{t(n)}];
        ioBin.mecLecPathCell{n} = [ioBin.mecRDualPathCell{t(n)}; ioBin.mecBDualPathCell{t(n)}];

        pprBin.lecLec{n} = [pprBin.lecBDual{n}(:,1); pprBin.lecRDual{n}(:,1)];
        pprBin.lecMec{n} = [pprBin.lecBDual{n}(:,2); pprBin.lecRDual{n}(:,2)];
        pprBin.lecLecPathCell{n} = [pprBin.lecBDualPathCell{n}; pprBin.lecRDualPathCell{n}];
        pprBin.lecMecPathCell{n} = [pprBin.lecBDualPathCell{n}; pprBin.lecRDualPathCell{n}];
        
        pprBin.mecMec{n} = [pprBin.mecRDual{n}(:,1); pprBin.mecBDual{n}(:,1)];
        pprBin.mecLec{n} = [pprBin.mecRDual{n}(:,2); pprBin.mecBDual{n}(:,2)];
        pprBin.mecMecPathCell{n} = [pprBin.mecRDualPathCell{n}; pprBin.mecBDualPathCell{n}];
        pprBin.mecLecPathCell{n} = [pprBin.mecRDualPathCell{n}; pprBin.mecBDualPathCell{n}];
    end
    % t = [1, 14];
    % for i = 1:length(t)
    %     ioBin.lecLec{i} = [ioBin.lecB{t(i)}(:,end); ioBin.lecR{t(i)}(:,end); ioBin.lecBDual{t(i)}(:,2); ioBin.lecRDual{t(i)}(:,2)];
    %     ioBin.lecMec{i} = [ioBin.lecBDual{t(i)}(:,4); ioBin.lecRDual{t(i)}(:,3)];
    %     ioBin.lecLecPathCell{i} = [ioBin.lecBPathCell{t(i)}; ioBin.lecRPathCell{t(i)}; ioBin.lecBDualPathCell{t(i)}; ioBin.lecRDualPathCell{t(i)}];
    %     ioBin.lecMecPathCell{i} = [ioBin.lecBDualPathCell{t(i)}; ioBin.lecRDualPathCell{t(i)}];
    % 
    %     ioBin.mecMec{i} = [ioBin.mecB{t(i)}(:,end); ioBin.mecR{t(i)}(:,end); ioBin.mecRDual{t(i)}(:,2); ioBin.mecBDual{t(i)}(:,2)];
    %     ioBin.mecLec{i} = [ioBin.mecRDual{t(i)}(:,3); ioBin.mecBDual{t(i)}(:,4)];
    %     ioBin.mecMecPathCell{i} = [ioBin.mecBPathCell{t(i)}; ioBin.mecRPathCell{t(i)}; ioBin.mecRDualPathCell{t(i)}; ioBin.mecBDualPathCell{t(i)}];
    %     ioBin.mecLecPathCell{i} = [ioBin.mecRDualPathCell{t(i)}; ioBin.mecBDualPathCell{t(i)}];
    % end

    spikeBin.lecAll =  [spikeBin.lecBDual;spikeBin.lecRDual];
    spikeBin.lecAllPathCell = [spikeBin.lecBDualPathCell; spikeBin.lecRDualPathCell];
    
    spikeBin.mecAll = [spikeBin.mecRDual; spikeBin.mecBDual];
    spikeBin.mecAllPathCell = [spikeBin.mecRDualPathCell;spikeBin.mecBDualPathCell];

end




fields = fieldnames(ioBin);
f = fields(contains(fields, 'PathCell'));

for i = 1:length(f)
    tempF = string(erase(f{i}, 'PathCell'));
    tempSd = string(strcat(tempF, 'Sd'));
    tempN = string(strcat(tempF, 'N'));
    tempSem = string(strcat(tempF, 'Sem'));

    [u, ~, g] = cellfun(@(x) unique(x(:,end), 'stable'), ioBin.(f{i}), 'UniformOutput', false);
    ioBinAvg.(tempF) = cellfun(@(data, groups) splitapply(@(x) mean(x,1, 'omitnan'), data, groups), ioBin.(tempF), g, 'UniformOutput', false);
    ioBinAvg.(tempSd) = cellfun(@(data, groups) splitapply(@(x) std(x,1, 'omitnan'), data, groups), ioBin.(tempF), g, 'UniformOutput', false);
    ioBinAvg.(tempN) = cellfun(@(data, groups) histc(data, unique(groups,'stable')), g, g, 'UniformOutput', false);
    ioBinAvg.(tempSem) = cellfun(@(data, groups) data./sqrt(groups), ioBinAvg.(tempSd), ioBinAvg.(tempN), 'UniformOutput', false);
    ioBinAvg.(f{i}) = u;
    
    if sum(matches(fieldnames(pprBin), tempF)) > 0
        [u, ~, g] = cellfun(@(x) unique(x(:,end), 'stable'), pprBin.(f{i}), 'UniformOutput', false);
        pprBinAvg.(tempF) = cellfun(@(data, groups) splitapply(@(x) mean(x,1, 'omitnan'), data, groups), pprBin.(tempF), g, 'UniformOutput', false);
        pprBinAvg.(tempSd) = cellfun(@(data, groups) splitapply(@(x) std(x,1, 'omitnan'), data, groups), pprBin.(tempF), g, 'UniformOutput', false);
        pprBinAvg.(tempN) = cellfun(@(data, groups) histc(data, unique(groups,'stable')), g, g, 'UniformOutput', false);
        pprBinAvg.(tempSem) = cellfun(@(data, groups) data./sqrt(groups), pprBinAvg.(tempSd), pprBinAvg.(tempN), 'UniformOutput', false);
        pprBinAvg.(f{i}) = u;
    end
   
    if contains(path, 'CC') && isfield(spikeBin, f(i))
        [u, ~, g] = unique(spikeBin.(f{i})(:,end), 'stable');
        spikeBinAvg.(tempF) = splitapply(@(x) mean(x,1, 'omitnan'), spikeBin.(tempF), g);
        spikeBinAvg.(tempSd) = splitapply(@(x) std(x,1, 'omitnan'), spikeBin.(tempF), g);
        spikeBinAvg.(tempN) = histc(g, unique(g, 'stable'));
        spikeBinAvg.(tempSem) = spikeBinAvg.(tempSd)./sqrt(spikeBinAvg.(tempN));
        spikeBinAvg.(f{i}) = u;
    end
end

if contains(path, 'VC')
    fields = fieldnames(eiBin);
    f = fields(contains(fields, 'PathCell'));
    for i = 1:length(f)
        tempF = string(erase(f{i}, 'PathCell'));
        tempSd = string(strcat(tempF, 'Sd'));
        tempN = string(strcat(tempF, 'N'));
        tempSem = string(strcat(tempF, 'Sem'));

        [u, ~, g] = unique(eiBin.(f{i})(:,end), 'stable');
        eiBinAvg.(tempF) = splitapply(@(x) mean(x,1, 'omitnan'), eiBin.(tempF), g);
        eiBinAvg.(tempSd) = splitapply(@(x) std(x,1, 'omitnan'), eiBin.(tempF), g);
        eiBinAvg.(tempN) = histc(g, unique(g, 'stable'));
        eiBinAvg.(tempSem) = eiBinAvg.(tempSd)./sqrt(eiBinAvg.(tempN));
        eiBinAvg.(f{i}) = u;

    end
end

fields = fieldnames(ioBin);
f = fields(contains(fields, 'DualPathCell'));
for i = 1:length(f)
    if contains(f{i}, "mecB") || contains(f{i}, "lecR")
        ratio.(f{i}) = cellfun(@(x) extractBefore(x(:,1), "mec"), ioBin.(f{i}), 'UniformOutput',false); 
        ratio.(f{i}) = cellfun(@(x,y) [x, y(:,end)], ratio.(f{i}), ioBin.(f{i}), 'UniformOutput',false);
        diff.(f{i}) = cellfun(@(x) extractBefore(x(:,1), "mec"), ioBin.(f{i}), 'UniformOutput',false); 
        diff.(f{i}) = cellfun(@(x,y) [x, y(:,end)], ratio.(f{i}), ioBin.(f{i}), 'UniformOutput',false);

    elseif contains(f{i}, "mecR") || contains(f{i}, "lecB")
        ratio.(f{i}) = cellfun(@(x) extractBefore(x(:,1), "lec"), ioBin.(f{i}), 'UniformOutput',false);
        ratio.(f{i}) = cellfun(@(x,y) [x, y(:,end)], ratio.(f{i}), ioBin.(f{i}), 'UniformOutput',false);
        diff.(f{i}) = cellfun(@(x) extractBefore(x(:,1), "lec"), ioBin.(f{i}), 'UniformOutput',false);
        diff.(f{i}) = cellfun(@(x,y) [x, y(:,end)], ratio.(f{i}), ioBin.(f{i}), 'UniformOutput',false);
    end
end

[tempIdx, tempIdx2] = cellfun(@(x,y) ismember(x(:,1), y(:,1)), ratio.lecBDualPathCell, ratio.mecRDualPathCell, 'UniformOutput', false);
tempIdx2 = cellfun(@(x) replaceZeros(x), tempIdx2, 'UniformOutput', false);
ratio.mecRlecB = cellfun(@(x,y,z,a) x(z,:)./y(a,:), ioBin.mecRDual, ioBin.lecBDual, tempIdx2, tempIdx, 'UniformOutput', false);
ratio.mecRlecBPathCell = cellfun(@(x,y) (x(y,:)), ratio.lecBDualPathCell, tempIdx, 'UniformOutput',false);

[tempIdx, tempIdx2] = cellfun(@(x,y) ismember(x(:,1),y(:,1)), ratio.lecRDualPathCell, ratio.mecBDualPathCell, 'UniformOutput', false);
tempIdx2 = cellfun(@(x) replaceZeros(x), tempIdx2, 'UniformOutput', false);
ratio.mecBlecR = cellfun(@(x,y,z,a) x(z,:)./y(a,:), ioBin.mecBDual, ioBin.lecRDual, tempIdx2, tempIdx, 'UniformOutput', false);
ratio.mecBlecRPathCell = cellfun(@(x,y) (x(y,:)), ratio.lecRDualPathCell, tempIdx, 'UniformOutput',false);

ratio.allDual  = cellfun(@(x, y) [x;y], ratio.mecBlecR, ratio.mecRlecB, 'UniformOutput',false);
ratio.allDualPathCell  = cellfun(@(x, y) [x;y], ratio.mecBlecRPathCell, ratio.mecRlecBPathCell, 'UniformOutput',false);

%%Difference (m-l/m+1)
[tempIdx, tempIdx2] = cellfun(@(x,y) ismember(x(:,1), y(:,1)), diff.lecBDualPathCell, diff.mecRDualPathCell, 'UniformOutput', false);
tempIdx2 = cellfun(@(x) replaceZeros(x), tempIdx2, 'UniformOutput', false);
diff.mecRlecB = cellfun(@(x,y,z,a) (abs(x(z,:))-abs(y(a,:)))./(abs(x(z,:))+abs(y(a,:))), ioBin.mecRDual, ioBin.lecBDual, tempIdx2, tempIdx, 'UniformOutput', false);
diff.mecRlecBPathCell = cellfun(@(x,y) (x(y,:)), diff.lecBDualPathCell, tempIdx, 'UniformOutput',false);

[tempIdx, tempIdx2] = cellfun(@(x,y) ismember(x(:,1),y(:,1)), diff.lecRDualPathCell, diff.mecBDualPathCell, 'UniformOutput', false);
tempIdx2 = cellfun(@(x) replaceZeros(x), tempIdx2, 'UniformOutput', false);
diff.mecBlecR = cellfun(@(x,y,z,a) (abs(x(z,:))-abs(y(a,:)))./(abs(x(z,:))+abs(y(a,:))), ioBin.mecBDual, ioBin.lecRDual, tempIdx2, tempIdx, 'UniformOutput', false);
diff.mecBlecRPathCell = cellfun(@(x,y) (x(y,:)), diff.lecRDualPathCell, tempIdx, 'UniformOutput',false);

diff.allDual  = cellfun(@(x, y) [x;y], diff.mecBlecR, diff.mecRlecB, 'UniformOutput',false);
diff.allDualPathCell  = cellfun(@(x, y) [x;y], diff.mecBlecRPathCell, diff.mecRlecBPathCell, 'UniformOutput',false);

fields = fieldnames(ratio);
f = fields(~endsWith(fields, 'Cell'));

for i = 1:length(f)
    tempF = strcat(f{i}, "PathCell");
    tempSd = strcat(f{i}, "Sd");
    tempN =  strcat(f{i}, "N");
    tempSem = strcat(f{i}, "Sem");

    [u, ~, g] = cellfun(@(x) unique(x(:,end), 'stable'), ratio.(tempF), 'UniformOutput', false);
    ratioAvg.(f{i}) = cellfun(@(data, groups) splitapply(@(x) mean(x,1, 'omitnan'), data, groups), ratio.(f{i}), g, 'UniformOutput', false);
    ratioAvg.(tempSd) = cellfun(@(data, groups) splitapply(@(x) std(x,1, 'omitnan'), data, groups),ratio.(f{i}), g, 'UniformOutput', false);
    ratioAvg.(tempN) = cellfun(@(data, groups) histc(data, unique(groups, 'stable')), g, g, 'UniformOutput', false);
    ratioAvg.(tempSem) = cellfun(@(data, groups) data./sqrt(groups), ratioAvg.(tempSd), ratioAvg.(tempN), 'UniformOutput', false);
    ratioAvg.(tempF) = u;

    [u, ~, g] = cellfun(@(x) unique(x(:,end), 'stable'), diff.(tempF), 'UniformOutput', false);
    diffAvg.(f{i}) = cellfun(@(data, groups) splitapply(@(x) mean(x,1, 'omitnan'), data, groups), diff.(f{i}), g, 'UniformOutput', false);
    diffAvg.(tempSd) = cellfun(@(data, groups) splitapply(@(x) std(x,1, 'omitnan'), data, groups),diff.(f{i}), g, 'UniformOutput', false);
    diffAvg.(tempN) = cellfun(@(data, groups) histc(data, unique(groups, 'stable')), g, g, 'UniformOutput', false);
    diffAvg.(tempSem) = cellfun(@(data, groups) data./sqrt(groups), diffAvg.(tempSd), diffAvg.(tempN), 'UniformOutput', false);
    diffAvg.(tempF) = u;
end




%% Plotting
fields = fieldnames(ioBinAvg);
f = fields(contains(fields, 'PathCell'));

for i = 1:length(f)
    if contains(f{i}, 'combo')
        ioBinAvg.(f{i}) = cellfun(@(x) [x,comboBin(:,end)], ioBinAvg.(f{i}), 'UniformOutput', false);
    else
        ioBinAvg.(f{i}) = cellfun(@(x) [x,mW_bin(:,end)], ioBinAvg.(f{i}), 'UniformOutput', false);
    end
end

if plt && contains(path, 'VC','IgnoreCase',true)
    for i = 1:length(kinetics)
        %Non Combined LEC v MEC Plots
        figure; 
        subplot(3,4, 1);
        errorbar(str2double(ioBinAvg.mecRPathCell{i}(:,2)), ioBinAvg.mecR{i}, ioBinAvg.mecRSem{i}, 'ro-','MarkerSize',10, 'MarkerFaceColor', 'r', 'LineWidth',2); hold on;
        errorbar(str2double(ioBinAvg.lecRPathCell{i}(:,2)), ioBinAvg.lecR{i}, ioBinAvg.lecRSem{i}, 'o-','MarkerSize',10,'Color','#289e1e', 'MarkerFaceColor', '#289e1e', 'LineWidth',2);
        legend('MECR', 'LECR'); ylabel(kinetics(i)); xlabel('LED Intensity (%)');
        set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3); title(['I/O LECR v. MEC-R', kinetics(i)]); 
        
        subplot(3,4,5); 
        plot(str2double(ioBin.mecRPathCell{i}(:,2)), ioBin.mecR{i}(:,1), 'r.', 'LineWidth', 1);
        
        subplot(3,4,9);
        plot(str2double(ioBin.lecRPathCell{i}(:,2)), ioBin.lecR{i}(:,1), '.','Color','#289e1e', 'LineWidth', 1);
        
        subplot(3,4,2);
        errorbar(str2double(ioBinAvg.mecBPathCell{i}(:,2)), ioBinAvg.mecB{i}(:,1), ioBinAvg.mecBSem{i}(:,1), 'ro-','MarkerSize',10, 'MarkerFaceColor', 'r', 'LineWidth',2); hold on;
        errorbar(str2double(ioBinAvg.lecBPathCell{i}(:,2)), ioBinAvg.lecB{i}(:,1), ioBinAvg.lecBSem{i}(:,1), 'o-','MarkerSize',10,'Color','#289e1e', 'MarkerFaceColor', '#289e1e', 'LineWidth',2);
        legend('MECB', 'LECB', 'FontSize', 16); xlabel('LED Intensity (%)');
        set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3); title(['I/O LECB v. MEC-B', kinetics(i)]); 
        
        subplot(3,4,6);
        plot(str2double(ioBin.mecBPathCell{i}(:,2)), ioBin.mecB{i}(:,1), 'r.', 'LineWidth', 1); 
        
        subplot(3,4,10);
        plot(str2double(ioBin.lecBPathCell{i}(:,2)), ioBin.lecB{i}(:,1), '.','Color','#289e1e', 'LineWidth', 1);
        
        subplot(3,4,3);
        errorbar(str2double(ioBinAvg.mecBDualPathCell{i}(:,2)), ioBinAvg.mecBDual{i}(:,1), ioBinAvg.mecBDualSem{i}(:,1), 'ro-','MarkerSize',10, 'MarkerFaceColor', 'r', 'LineWidth',2); hold on;
        errorbar(str2double(ioBinAvg.lecRDualPathCell{i}(:,2)), ioBinAvg.lecRDual{i}(:,1), ioBinAvg.lecRDualSem{i}(:,1), 'o-','MarkerSize',10,'Color','#289e1e', 'MarkerFaceColor', '#289e1e', 'LineWidth',2);
        legend('MECB', 'LECR', 'FontSize', 16); xlabel('LED Intensity (%)');
        set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3); title(['I/O MECBxLECR Dual', kinetics(i)]); 
        
        subplot(3,4,7);
        plot(str2double(ioBin.mecBDualPathCell{i}(:,2)), ioBin.mecBDual{i}(:,1), 'r.', 'LineWidth', 1); 
        
        subplot(3,4,11);
        plot(str2double(ioBin.lecRDualPathCell{i}(:,2)), ioBin.lecRDual{i}(:,1), '.','Color','#289e1e', 'LineWidth', 1);
        
        subplot(3,4,4);
        errorbar(str2double(ioBinAvg.mecRDualPathCell{i}(:,2)), ioBinAvg.mecRDual{i}(:,1), ioBinAvg.mecRDualSem{i}(:,1), 'ro-','MarkerSize',10, 'MarkerFaceColor', 'r', 'LineWidth',2); hold on;
        errorbar(str2double(ioBinAvg.lecBDualPathCell{i}(:,2)), ioBinAvg.lecBDual{i}(:,1), ioBinAvg.lecBDualSem{i}(:,1), 'o-','MarkerSize',10,'Color','#289e1e', 'MarkerFaceColor', '#289e1e', 'LineWidth',2);
        legend('MECR', 'LECB', 'FontSize', 16); xlabel('LED Intensity (%)');
        set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3); title(['I/O MECRxLECB Dual', kinetics(i)]); 
        
        subplot(3,4,8);
        plot(str2double(ioBin.mecRDualPathCell{i}(:,2)), ioBin.mecRDual{i}(:,1), 'r.', 'LineWidth', 1); 
        
        subplot(3,4,12);
        plot(str2double(ioBin.lecBDualPathCell{i}(:,2)), ioBin.lecBDual{i}(:,1), '.','Color','#289e1e', 'LineWidth', 1);
        
        %Combined LEC Indv or LEC Dual labelled plots
        figure;
        subplot(3,2,1);
        errorbar(str2double(ioBinAvg.mecPathCell{i}(:,2)), ioBinAvg.mec{i}(:,1), ioBinAvg.mecSem{i}(:,1), 'ro-','MarkerSize',10, 'MarkerFaceColor', 'r', 'LineWidth',2); hold on;
        errorbar(str2double(ioBinAvg.lecPathCell{i}(:,2)), ioBinAvg.lec{i}(:,1), ioBinAvg.lecSem{i}(:,1), 'o-','MarkerSize',10,'Color','#289e1e', 'MarkerFaceColor', '#289e1e', 'LineWidth',2);
        legend('MEC R or B', 'LEC R or B', 'FontSize', 16); xlabel('LED Intensity (%)'); ylabel(kinetics(i));
        set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3); title(['I/O MEC LEC (indv Labelled)', kinetics(i)]); 
        
        subplot(3,2,3);
        plot(str2double(ioBin.mecPathCell{i}(:,2)), ioBin.mec{i}(:,1), 'r.', 'LineWidth', 1); 
        
        subplot(3,2,5)
        plot(str2double(ioBin.lecPathCell{i}(:,2)), ioBin.lec{i}(:,1), '.','Color','#289e1e', 'LineWidth', 1);
        
        subplot(3,2,2);
        errorbar(str2double(ioBinAvg.mecAllPathCell{i}(:,2)), ioBinAvg.mecAll{i}(:,1), ioBinAvg.mecAllSem{i}(:,1), 'ro-','MarkerSize',10, 'MarkerFaceColor', 'r', 'LineWidth',2); hold on;
        errorbar(str2double(ioBinAvg.lecAllPathCell{i}(:,2)), ioBinAvg.lecAll{i}(:,1), ioBinAvg.lecAllSem{i}(:,1), 'o-','MarkerSize',10,'Color','#289e1e', 'MarkerFaceColor', '#289e1e', 'LineWidth',2);
        legend('MEC R or B', 'LEC R or B', 'FontSize', 16); xlabel('LED Intensity (%)'); 
        set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3); title(['I/O MEC LEC (All Combined)', kinetics(i)]); 
        
        subplot(3,2,4);
        plot(str2double(ioBin.mecAllPathCell{i}(:,2)), ioBin.mecAll{i}(:,1), 'r.', 'LineWidth', 1); 
        
        subplot(3,2,6)
        plot(str2double(ioBin.lecAllPathCell{i}(:,2)), ioBin.lecAll{i}(:,1), '.','Color','#289e1e', 'LineWidth', 1);
    
    end
    if plt2
        k = find(contains(kinetics, {'amp_io', 'sum_io'}));
        for i = 1:length(k)
            i = k(i); 

            figure;
            subplot(2,3,1)
            errorbar(str2double(ratioAvg.mecRlecBPathCell{i}(:,1)), ratioAvg.mecRlecB{i}(:,1), ratioAvg.mecRlecBSem{i}(:,1),'.-k','MarkerSize', 10, 'LineWidth',2);
            xlabel('LED Intensity (%)'); ylabel('Ratio MEC/LEC');
            set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3); title(['Ratio MECR/LECB', kinetics(i)], 'Interpreter','none'); 
            
            subplot(2,3,4)
            plot(str2double(ratio.mecRlecBPathCell{i}(:,2)), ratio.mecRlecB{i}(:,1), '.','Color', [0.7 0.7 0.7], 'LineWidth', 1); 
            
            subplot(2,3,2)
            errorbar(str2double(ratioAvg.mecBlecRPathCell{i}(:,1)), ratioAvg.mecBlecR{i}(:,1), ratioAvg.mecBlecRSem{i}(:,1),'.-k','MarkerSize', 10, 'LineWidth',2);
            xlabel('LED Intensity (%)'); ylabel('Ratio MEC/LEC');
            set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3); title(['Ratio MECB/LECR', kinetics(i)], 'Interpreter','none'); 
            
            subplot(2,3,5)
            plot(str2double(ratio.mecBlecRPathCell{i}(:,2)), ratio.mecBlecR{i}(:,1), '.','Color', [0.7 0.7 0.7], 'LineWidth', 1); 
            
            subplot(2,3,3)
            errorbar(str2double(ratioAvg.allDualPathCell{i}(:,1)), ratioAvg.allDual{i}(:,1), ratioAvg.allDualSem{i}(:,1),'.-k','MarkerSize', 10, 'LineWidth',2);
            xlabel('LED Intensity (%)'); ylabel('Ratio MEC/LEC');
            set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3); title(['Ratio MEC/LEC', kinetics(i)], 'Interpreter','none'); 
            
            subplot(2,3,6)
            plot(str2double(ratio.allDualPathCell{i}(:,2)), ratio.allDual{i}(:,1), '.','Color', [0.7 0.7 0.7], 'LineWidth', 1); 
        end
      
        for i = 1:length(k)
            i = k(i); 

            figure;
            subplot(2,3,1)
            errorbar(str2double(diffAvg.mecRlecBPathCell{i}(:,1)), diffAvg.mecRlecB{i}(:,1), diffAvg.mecRlecBSem{i}(:,1),'.-k','MarkerSize', 10, 'LineWidth',2);
            ylim([-1, 1]);
            xlabel('LED Intensity (%)'); ylabel('Ratio m-l/m+1');
            set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3); title(['Ratio MECR-LECB/MECR+LECB', kinetics(i)], 'Interpreter','none'); 
            
            subplot(2,3,4)
            plot(str2double(diff.mecRlecBPathCell{i}(:,2)), diff.mecRlecB{i}(:,1), '.','Color', [0.7 0.7 0.7], 'LineWidth', 1); 
            
            subplot(2,3,2)
            errorbar(str2double(diffAvg.mecBlecRPathCell{i}(:,1)), diffAvg.mecBlecR{i}(:,1), diffAvg.mecBlecRSem{i}(:,1),'.-k','MarkerSize', 10, 'LineWidth',2);
            ylim([-1, 1]);
            xlabel('LED Intensity (%)'); ylabel('Ratio m-l/m+1');
            set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3); title(['Ratio MECB-LECR/MECB+LECR', kinetics(i)],'Interpreter','none'); 
            
            subplot(2,3,5)
            plot(str2double(diff.mecBlecRPathCell{i}(:,2)), diff.mecBlecR{i}(:,1), '.','Color', [0.7 0.7 0.7], 'LineWidth', 1); 
            
            subplot(2,3,3)
            errorbar(str2double(diffAvg.allDualPathCell{i}(:,1)), diffAvg.allDual{i}(:,1), diffAvg.allDualSem{i}(:,1),'.-k','MarkerSize', 10, 'LineWidth',2);
            ylim([-1, 1]);
            xlabel('LED Intensity (%)'); ylabel('Ratio m-l/m+1');
            set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3); title(['Ratio MEC-LEC/MEC+LEC', kinetics(i)], 'Interpreter','none'); 
            
            subplot(2,3,6)
            plot(str2double(ratio.allDualPathCell{i}(:,2)), diff.allDual{i}(:,1), '.','Color', [0.7 0.7 0.7], 'LineWidth', 1); 
        end
    end
end

pprBin = rmfield(pprBin, 'comboExp'); pprBin = rmfield(pprBin, 'comboExpPathCell'); 

cd(path);
if contains(path, 'VC','IgnoreCase',true)
    save('Compile_Kinetics_Binned', 'io', 'ioAvg', 'ioBin', 'ioBinAvg','ratio','ratioAvg', 'diff', 'diffAvg', 'kinetics', 'mW_bin', 'ei', 'eiAvg', 'eiBin', 'eiBinAvg', 'eiCell',...
        'ppr', 'pprAvg', 'pprBin','pprBinAvg');
    %save('ioBinAvg', 'ioBinAvg');
elseif contains(path, 'CC','IgnoreCase',true)
    save('Compile_Kinetics_Binned', 'io', 'ioAvg', 'ioBin', 'ioBinAvg','ratio','ratioAvg', 'diff', 'diffAvg', 'kinetics', 'mW_bin','spikeProbAvg','spikeProb','spikeBin','spikeBinAvg',...
        'ppr', 'pprAvg', 'pprBin','pprBinAvg');
end
%Add lec on one axis v mec on other axis plot
