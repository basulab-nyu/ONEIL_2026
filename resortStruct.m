%% Reorder For Stats/Plotting 
% rows = cell
% columns = intensity or other variable
clear all; close all;


%SORTING WRONG getting error in CC check what is going on - 09/15

path = ['/Volumes/KO_Portable/MEC_LEC_CA3/analysis/CC/IO']; cd(path);
load('Compile_Kinetics_Binned.mat');
%load('ioBinAvg.mat');

path2 = ['/Volumes/KO_Portable/MEC_LEC_CA3/analysis/mW']; cd(path2);
load('led_power.mat');


if contains(path, 'CC')
    load('Compile_Kinetics_CellSpikeProb.mat');
    s = ["ioBin", "spikeBin", "io", "diff","ratio","spikeProb", "spikeOnly", "pprBin"]; %took out spikecelldata "spikeCellData",
    pattern = ["Bin", "diff", "ratio", "Cell","Only"]; %this is for structures that don't have "stim information"
elseif contains(path,'VC')
    s = ["ioBin", "io", "diff","ratio","eiCell","eiBin", "pprBin"];
    pattern = ["Bin", "diff", "ratio"]; %this is for structures that don't have "stim information"
end

p = ["Path", "Stim", "Std", "Num", "Sem"];
for i = 1:length(s)
    tempS = eval(s(i));
    f = fieldnames(tempS);
    f = f(~contains(f, p));
    for ii = 1:length(f)
        tempF = string(strcat(f(ii),'PathCell'));
        if iscell(tempS.(tempF)) && ~isempty(tempS.(tempF){1})
            tempIds = cellfun(@(x) extractAfter(x(:,1), 'K'), tempS.(tempF), 'UniformOutput', false);
            [u,~,g] = cellfun(@(x) unique(x,'stable'), tempIds,'UniformOutput', false);
            for iii = 1:numel(tempS.(tempF))
                if ~isempty(tempS.(tempF){iii}) & contains(s(i), pattern)
                    tempIdx = str2double(tempS.(tempF){iii}(:,2));
                    out.(tempF){iii} = string();
                    out.(f{ii}){iii} = [];
                    for j = 1:length(g{iii})
                        out.(tempF){iii}(g{iii}(j),tempIdx(j)) = tempS.(tempF){iii}(j,1);
                        out.(f{ii}){iii}(g{iii}(j),tempIdx(j)) = tempS.(f{ii}){iii}(j,1);
                    end
                elseif ~isempty(tempS.(tempF){iii})
                    tempF2 = strcat(f{ii}, 'Stim');
                    tempIdx = tempS.(tempF2){iii}(:,2);
                    out.(tempF){iii} = string();
                    out.(f{ii}){iii} = [];
                    for j = 1:length(g{iii})
                        out.(tempF){iii}(g{iii}(j),tempIdx(j)) = strcat(string(tempIdx(j)), tempS.(tempF){iii}(j,1));
                        out.(f{ii}){iii}(g{iii}(j),tempIdx(j)) = tempS.(f{ii}){iii}(j,1);
                    end
                end
            end
        elseif ~isempty(tempS.(tempF){1})
            tempIds = extractAfter(tempS.(tempF)(:,1), 'K');
            [u,~,g] = unique(tempIds,'stable');
            if ~isempty(tempS.(tempF)) & contains(s(i), pattern)
                tempIdx = str2double(tempS.(tempF)(:,end));
                out.(tempF) = string();
                out.(f{ii}) = [];
                for j = 1:length(g)
                    out.(tempF)(g(j), tempIdx(j)) = tempS.(tempF)(j,1);
                    out.(f{ii})(g(j), tempIdx(j)) = tempS.(f{ii})(j,1);
                end
            elseif ~isempty(tempS.(tempF))
                tempF2 = strcat(f{ii}, 'Stim');
                tempIdx = tempS.(tempF2)(:,2);
                out.(tempF) = string();
                out.(f{ii}) = [];
                for j = 1:length(g)
                    out.(tempF)(g(j), tempIdx(j)) = strcat(string(tempIdx(j)), tempS.(tempF)(j,1));
                    out.(f{ii})(g(j), tempIdx(j)) = tempS.(f{ii})(j,1);    
                end
            end
        end
        if isempty(out) || ~isfield(out, f{ii})
            continue
        elseif  isa(out.(f{ii}), 'cell')   
            tempMiss = cellfun(@(x) ismissing(x), out.(tempF), 'UniformOutput', false);
            out.(f{ii}) = cellfun(@(x,y) setNaN(x,y), out.(f{ii}), tempMiss, 'UniformOutput', false);
        elseif ~isempty(out) && isa(out.(f{ii}), 'double')
            tempMiss = ismissing(out.(tempF));
            out.(f{ii})(tempMiss) = NaN;
        end
    end
    assignin('base', strcat(s(i), 'Sort'), out);
    out = [];
end

vars = who;
vars = vars(endsWith(vars,'Sort'));

%% Heatmap Values
if contains(path, 'CC')
    p = ["Cell", "Stim"];
    f = fieldnames(spikeOnly);
    f = f(~contains(f, p));
    for i = 1:length(f)
        tempF = strcat(f{i}, 'PathCell');
        tempF2 = strcat(f{i}, 'Val');
        for ii = 1:spikeCellAvg.allNum
           if ~contains(tempF, 'combo', 'IgnoreCase',true)
               tempIdx = logical(spikeOnlySort.(f{i})(ii,:));
               spikeOnlySort.(tempF2)(ii,1:sum(tempIdx)) = mW_bin(tempIdx,4);
           elseif contains(tempF, 'combo', 'IgnoreCase',true)
               tempIdx = logical(spikeOnlySort.(f{i})(ii,comboBin(:,2)));
               temp(ii,1:length(tempIdx)) = spikeOnlySort.(f{i})(ii,comboBin(:,2));
               spikeOnlySort.(tempF2)(ii,1:sum(tempIdx)) = comboBin(tempIdx,4)';
           end
        end
        spikeOnlySort.(tempF2)(find(spikeOnlySort.(tempF2) == 0)) = NaN;
    end
    spikeOnlySort.combo = temp;
end

cd(path); save('Kinetics_Resort.mat', vars{:});