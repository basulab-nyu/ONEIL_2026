%% After running indv_kinetics, compile and average
% Last Edit: 6/6/2025

%Something weird happening w/ LECR Dual in both VC and CC

clear all; close all; 

path_2 = ['/Volumes/KO_Portable/MEC_LEC_CA3/analysis/CC/IO']; cd(path_2);

files = dir('2*ioStruct.mat');
plt  = 0;

io.lecB = [];
io.lecBStim = [];
io.lecBStimId = [];
io.lecR = [];
io.lecRStim = [];
io.lecRStimId = [];
io.lecBDual = [];
io.lecBDualStim = [];
io.lecBDualStimId = [];
io.lecRDual = [];
io.lecRDualStim = [];
io.lecRDualStimId = [];

io.mecB = [];
io.mecBStim = [];
io.mecBStimId = [];
io.mecR = [];
io.mecRStim = [];
io.mecRStimId = [];
io.mecBDual = [];
io.mecBDualStim = [];
io.mecBDualStimId = [];
io.mecRDual = [];
io.mecRDualStim = [];
io.mecRDualStimId = [];

io.mecRDualPath = [];
io.mecRPath = [];
io.mecBDualPath = [];
io.lecRDualPath = [];
io.mecBPath = [];
io.lecBDualPath = [];
io.lecRPath = [];
io.lecBPath = [];

io.mecRDualPathCell = [];
io.mecRPathCell = [];
io.mecBDualPathCell = [];
io.lecRDualPathCell = [];
io.mecBPathCell = [];
io.lecBDualPathCell = [];
io.lecRPathCell = [];
io.lecBPathCell = [];


io.combo = [];
io.comboPath = [];
io.comboPathCell = [];
io.comboStim = [];
io.comboStimId = [];

io.comboExp = [];
io.comboExpStim = [];

if contains(path_2, 'vc', 'IgnoreCase',true)
    kinetics = ["amp_io10", "amp_io80", "ar19_io10", "ar19_io80", "bsT_io10", "bsT_io80", "dc19_io10", ...
    "dc19_io80", "dc28_io10", "dc28_io80", "hw_io10", "hw_io80", "lagB_io10", "lagB_io80",...
    "lagR_io10", "lagR_io80", "ostB_io10", "ostB_io80", "ostR_io10", "ostR_io80", "pkT_io10", ...
    "pkT_io80", "rs19_io10", "rs19_io80", "rs28_io10", "rs28_io80", "sum_io10", "sum_io80"];
end

if contains(path_2, 'cc', 'IgnoreCase',true)
    kinetics = ["amp_io70", "ar19_io70", "bsT_io70", "dc19_io70", "dc28_io70","hw_io70", "lagB_io70",...
    "lagR_io70", "ostB_io70", "ostR_io70", "pkT_io70", "rs19_io70", "rs28_io70", "sum_io70"];
end

a = string(fieldnames(io));
for i = 1:length(kinetics)
    for ii = 1:length(fieldnames(io))
        io.(a(ii)){i} = [];
        ppr.(a(ii)){i} = [];
    end
end



%% Get all files together

for i = 1:length(files)
    load(files(i).name);

    if length(ioStruct.Animal) ~= 5 %Add a zero to animals <100 to make same size 
        ioStruct.Animal = [ioStruct.Animal(1:2), '0', ioStruct.Animal(3:4)];
    end
    if sum(contains(ioStruct.kinetics, '_io_')) && contains(path_2, 'cc', 'IgnoreCase', true)
       ioStruct.kinetics = insertAfter(ioStruct.kinetics, "io", "70")
    end    

    %uStim = unique(ioStruct.kinetics);
    if isfield(ioStruct, 'dataAvg')
        if sum(contains(ioStruct.kinetics, 'io70')) > 1 && contains(path_2, 'vc', 'IgnoreCase', true)
            ioStruct.kinetics = strrep(ioStruct.kinetics, 'io70', 'io80');
        end
        if string(ioStruct.blueId) == "lec"
            if ioStruct.redId == 0
                for ii = 1:length(ioStruct.kinetics)
                    temp_idx = find(kinetics == extractBefore(ioStruct.kinetics(ii), '_2'), 1);
                    io.lecB{temp_idx} = [io.lecB{temp_idx}; ioStruct.dataBlueAvg{ii}];
                    io.lecBStim{temp_idx} = [io.lecBStim{temp_idx}; ioStruct.dataAvgStim{ii}(:,1)];
                    io.lecBStimId{temp_idx} = [io.lecBStimId{temp_idx}; string(ioStruct.dataBlueStimId{ii}(1:2))];
                    io.lecBPath{temp_idx} = [io.lecBPath{temp_idx}; ioStruct.Animal, ioStruct.Cell, ioStruct.blueId, ioStruct.redId, ioStruct.clamp];
                    io.lecBPathCell{temp_idx} = [ io.lecBPathCell{temp_idx}; string(repmat([ioStruct.Animal, ioStruct.Cell, ioStruct.blueId, ioStruct.redId],length(ioStruct.dataAvgStim{ii}(:,1)),1))];

                    if ~isempty(ioStruct.pprAvg{temp_idx})
                        t = contains(ioStruct.pprStimId{temp_idx}, 'b');
                        ppr.lecB{temp_idx} = [ppr.lecB{temp_idx}; ioStruct.pprAvg{temp_idx}(:,t)]; 
                        ppr.lecBStimId{temp_idx} = [ppr.lecBStimId{temp_idx}; ioStruct.pprStimId{temp_idx}(:,t)]; 
                        ppr.lecBPath{temp_idx} = [ppr.lecBPath{temp_idx}; io.lecBPath{temp_idx}] ;
                        ppr.lecBPathCell{temp_idx} =[ppr.lecBPathCell{temp_idx};  io.lecBPathCell{temp_idx}];   
                    end
                end
                %io.lecBPath = [io.lecBPath; ioStruct.Animal, ioStruct.Cell, ioStruct.blueId, ioStruct.redId, ioStruct.clamp];
                
    
            elseif string(ioStruct.redId) == "mec"
               for ii = 1:length(ioStruct.kinetics)
                    temp_idx = find(kinetics == extractBefore(ioStruct.kinetics(ii), '_2'), 1);
                    io.lecBDual{temp_idx} = [io.lecBDual{temp_idx}; ioStruct.dataBlueAvg{ii}, ioStruct.dataBlue{ii}(:,2:end)];
                    io.lecBDualStim{temp_idx} = [io.lecBDualStim{temp_idx}; ioStruct.dataAvgStim{ii}(:,1)];
                    io.lecBDualStimId{temp_idx} = [io.lecBDualStimId{temp_idx}; string(ioStruct.dataBlueStimId{ii})];
    
                    io.mecRDual{temp_idx} = [io.mecRDual{temp_idx}; ioStruct.dataRedAvg{ii}, ioStruct.dataRed{ii}(:,2:end)];
                    io.mecRDualStim{temp_idx} = [io.mecRDualStim{temp_idx}; ioStruct.dataAvgStim{ii}(:,1)];
                    io.mecRDualStimId{temp_idx} = [io.mecRDualStimId{temp_idx}; string(ioStruct.dataRedStimId{ii})];
                    io.lecBDualPath{temp_idx} = [io.lecBDualPath{temp_idx}; ioStruct.Animal, ioStruct.Cell, ioStruct.blueId, ioStruct.redId, ioStruct.clamp];
                    io.mecRDualPath{temp_idx} = [io.mecRDualPath{temp_idx}; ioStruct.Animal, ioStruct.Cell, ioStruct.blueId, ioStruct.redId, ioStruct.clamp];
                    io.lecBDualPathCell{temp_idx} = [io.lecBDualPathCell{temp_idx}; string(repmat([ioStruct.Animal, ioStruct.Cell, ioStruct.blueId, ioStruct.redId],length(ioStruct.dataAvgStim{ii}(:,1)),1))];
                    io.mecRDualPathCell{temp_idx} = [io.mecRDualPathCell{temp_idx}; string(repmat([ioStruct.Animal, ioStruct.Cell, ioStruct.blueId, ioStruct.redId],length(ioStruct.dataAvgStim{ii}(:,1)),1))];
                    
                    io.combo{temp_idx} = [io.combo{temp_idx}; ioStruct.dataCombo{ii}];
                    io.comboStim{temp_idx} = [io.comboStim{temp_idx}; ioStruct.dataAvgStim{ii}(:,1)];
                    io.comboPath{temp_idx} = [io.comboPath{temp_idx}; ioStruct.Animal, ioStruct.Cell, ioStruct.blueId, ioStruct.redId, ioStruct.clamp];
                    io.comboPathCell{temp_idx} = [io.comboPathCell{temp_idx}; string(repmat([ioStruct.Animal, ioStruct.Cell, ioStruct.blueId, ioStruct.redId],length(ioStruct.dataAvgStim{ii}(:,1)),1))];
                    
                    io.comboExp{temp_idx} = [io.comboExp{temp_idx}; (ioStruct.dataBlueAvg{ii}+ioStruct.dataRedAvg{ii})];
                    io.comboExpStim{temp_idx} = [io.comboExpStim{temp_idx}; ioStruct.dataAvgStim{ii}(:,1)];
                    io.comboExpPathCell{temp_idx} = io.comboPathCell{temp_idx};
                    
                    if ~isempty(ioStruct.pprAvg{temp_idx})
                        t = matches(ioStruct.pprStimId{temp_idx}, 'bb');
                        t2 = matches(ioStruct.pprStimId{temp_idx}, 'rb');

                        ppr.lecBDual{temp_idx} = [ppr.lecBDual{temp_idx}; ioStruct.pprAvg{temp_idx}(:,t), ioStruct.pprAvg{temp_idx}(:,t2)]; 
                        ppr.lecBDualStimId{temp_idx} = [ppr.lecBDualStimId{temp_idx}; ioStruct.pprStimId{temp_idx}(:,t),ioStruct.pprStimId{temp_idx}(:,t2)]; 
                    
                        ppr.lecBDualPath{temp_idx} = [ppr.lecBDualPath{temp_idx}; io.lecBPath{temp_idx}];
                        ppr.lecBDualPathCell{temp_idx} = [ppr.lecBDualPathCell{temp_idx}; io.lecBPathCell{temp_idx}]; 

                        t = matches(ioStruct.pprStimId{temp_idx}, 'rr');
                        t2 = matches(ioStruct.pprStimId{temp_idx}, 'br');
                        ppr.mecRDual{temp_idx} = [ppr.mecRDual{temp_idx}; ioStruct.pprAvg{temp_idx}(:,t), ioStruct.pprAvg{temp_idx}(:,t2)]; 
                        ppr.mecRDualStimId{temp_idx} = [ppr.mecRDualStimId{temp_idx}; ioStruct.pprStimId{temp_idx}(:,t),ioStruct.pprStimId{temp_idx}(:,t2)]; 
                    
                        ppr.mecRDualPath{temp_idx} = [ppr.mecRDualPath{temp_idx}; io.mecRPath{temp_idx}];
                        ppr.mecRDualPathCell{temp_idx} = [ppr.mecRDualPathCell{temp_idx}; io.mecRPathCell{temp_idx}]; 
                       
                        t = matches(ioStruct.pprStimId{temp_idx}, 'brbr');
                        ppr.combo{temp_idx} = [ppr.combo{temp_idx}; ioStruct.pprAvg{temp_idx}(:,t)];
                        ppr.comboStimId{temp_idx} = [ppr.comboStimId{temp_idx}; ioStruct.pprStimId{temp_idx}(:,t)]; 
                        ppr.comboPath{temp_idx} = [ppr.comboPath{temp_idx}; io.comboPath{temp_idx}];
                        ppr.comboPathCell{temp_idx} = [ppr.comboPathCell{temp_idx}; io.comboPathCell{temp_idx}];
                    end
               end

             
            end
        elseif string(ioStruct.blueId) == "mec"
            if ioStruct.redId == 0
                for ii = 1:length(ioStruct.kinetics)
                    temp_idx = find(kinetics == extractBefore(ioStruct.kinetics(ii), '_2'), 1);
                    io.mecB{temp_idx} = [io.mecB{temp_idx}; ioStruct.dataBlueAvg{ii}];
                    io.mecBStim{temp_idx} = [io.mecBStim{temp_idx}; ioStruct.dataAvgStim{ii}(:,1)];
                    io.mecBStimId{temp_idx} = [io.mecBStimId{temp_idx}; string(ioStruct.dataBlueStimId{ii}(1:2))];
                    io.mecBPath{temp_idx} = [io.mecBPath{temp_idx}; ioStruct.Animal, ioStruct.Cell, ioStruct.blueId, ioStruct.redId, ioStruct.clamp];
                    io.mecBPathCell{temp_idx} = [io.mecBPathCell{temp_idx}; string(repmat([ioStruct.Animal, ioStruct.Cell, ioStruct.blueId, ioStruct.redId],length(ioStruct.dataAvgStim{ii}(:,1)),1))];
                   
                    if ~isempty(ioStruct.pprAvg{temp_idx})
                        t = contains(ioStruct.pprStimId{temp_idx}, 'b');
                        ppr.mecB{temp_idx} = [ppr.mecB{temp_idx}; ioStruct.pprAvg{temp_idx}(:,t)]; 
                        ppr.mecBStimId{temp_idx} = [ppr.mecBStimId{temp_idx}; ioStruct.pprStimId{temp_idx}(:,t)]; 
                        ppr.mecBPath{temp_idx} = [ppr.mecBPath{temp_idx}; io.mecBPath{temp_idx}] ;
                        ppr.mecBPathCell{temp_idx} =[ppr.mecBPathCell{temp_idx};  io.mecBPathCell{temp_idx}];   
                    end
                   
                end
                %io.mecBPath = [io.mecBPath; ioStruct.Animal, ioStruct.Cell, ioStruct.blueId, ioStruct.redId, ioStruct.clamp];
                
    
            elseif string(ioStruct.redId) == "lec"
               for ii = 1:length(ioStruct.kinetics)
                    temp_idx = find(kinetics == extractBefore(ioStruct.kinetics(ii), '_2'), 1);
                    io.mecBDual{temp_idx} = [io.mecBDual{temp_idx}; ioStruct.dataBlueAvg{ii}, ioStruct.dataBlue{ii}(:,2:end)];
                    io.mecBDualStim{temp_idx} = [io.mecBDualStim{temp_idx}; ioStruct.dataAvgStim{ii}(:,1)];
                    io.mecBDualStimId{temp_idx} = [io.mecBDualStimId{temp_idx}; string(ioStruct.dataBlueStimId{ii})];
    
                    io.lecRDual{temp_idx} = [io.lecRDual{temp_idx}; ioStruct.dataRedAvg{ii}, ioStruct.dataRed{ii}(:,2:end)];
                    io.lecRDualStim{temp_idx} = [io.lecRDualStim{temp_idx}; ioStruct.dataAvgStim{ii}(:,1)];
                    io.lecRDualStimId{temp_idx} = [io.lecRDualStimId{temp_idx}; string(ioStruct.dataRedStimId{ii})];
                   
                    io.lecRDualPath{temp_idx} = [io.lecRDualPath{temp_idx}; ioStruct.Animal, ioStruct.Cell, ioStruct.blueId, ioStruct.redId, ioStruct.clamp];
                    io.mecBDualPath{temp_idx} = [io.mecBDualPath{temp_idx}; ioStruct.Animal, ioStruct.Cell, ioStruct.blueId, ioStruct.redId, ioStruct.clamp];
                    io.lecRDualPathCell{temp_idx} = [io.lecRDualPathCell{temp_idx}; string(repmat([ioStruct.Animal, ioStruct.Cell, ioStruct.blueId, ioStruct.redId],length(ioStruct.dataAvgStim{ii}(:,1)),1))];
                    io.mecBDualPathCell{temp_idx} = [io.mecBDualPathCell{temp_idx}; string(repmat([ioStruct.Animal, ioStruct.Cell, ioStruct.blueId, ioStruct.redId],length(ioStruct.dataAvgStim{ii}(:,1)),1))];
                    

                    io.combo{temp_idx} = [io.combo{temp_idx}; ioStruct.dataCombo{ii}];
                    io.comboStim{temp_idx} = [io.comboStim{temp_idx}; ioStruct.dataAvgStim{ii}(:,1)];
                    io.comboPath{temp_idx} = [io.comboPath{temp_idx}; ioStruct.Animal, ioStruct.Cell, ioStruct.blueId, ioStruct.redId, ioStruct.clamp];
                    io.comboPathCell{temp_idx} = [io.comboPathCell{temp_idx}; string(repmat([ioStruct.Animal, ioStruct.Cell, ioStruct.blueId, ioStruct.redId],length(ioStruct.dataAvgStim{ii}(:,1)),1))];

                    io.comboExp{temp_idx} = [io.comboExp{temp_idx}; (ioStruct.dataBlueAvg{ii}+ioStruct.dataRedAvg{ii})];
                    io.comboExpStim{temp_idx} = [io.comboExpStim{temp_idx}; ioStruct.dataAvgStim{ii}(:,1)];
                    io.comboExpPathCell{temp_idx} = io.comboPathCell{temp_idx};
                     if ~isempty(ioStruct.pprAvg{temp_idx})
                        t = matches(ioStruct.pprStimId{temp_idx}, 'bb');
                        t2 = matches(ioStruct.pprStimId{temp_idx}, 'rb');

                        ppr.mecBDual{temp_idx} = [ppr.mecBDual{temp_idx}; ioStruct.pprAvg{temp_idx}(:,t), ioStruct.pprAvg{temp_idx}(:,t2)]; 
                        ppr.mecBDualStimId{temp_idx} = [ppr.mecBDualStimId{temp_idx}; ioStruct.pprStimId{temp_idx}(:,t),ioStruct.pprStimId{temp_idx}(:,t2)]; 
                    
                        ppr.mecBDualPath{temp_idx} = [ppr.mecBDualPath{temp_idx}; io.mecBPath{temp_idx}];
                        ppr.mecBDualPathCell{temp_idx} = [ppr.mecBDualPathCell{temp_idx}; io.mecBPathCell{temp_idx}]; 

                        t = matches(ioStruct.pprStimId{temp_idx}, 'rr');
                        t2 = matches(ioStruct.pprStimId{temp_idx}, 'br');
                        ppr.lecRDual{temp_idx} = [ppr.lecRDual{temp_idx}; ioStruct.pprAvg{temp_idx}(:,t), ioStruct.pprAvg{temp_idx}(:,t2)]; 
                        ppr.lecRDualStimId{temp_idx} = [ppr.lecRDualStimId{temp_idx}; ioStruct.pprStimId{temp_idx}(:,t),ioStruct.pprStimId{temp_idx}(:,t2)]; 
                    
                        ppr.lecRDualPath{temp_idx} = [ppr.lecRDualPath{temp_idx}; io.lecRPath{temp_idx}];
                        ppr.lecRDualPathCell{temp_idx} = [ppr.lecRDualPathCell{temp_idx}; io.lecRPathCell{temp_idx}]; 
                       
                        t = matches(ioStruct.pprStimId{temp_idx}, 'brbr');
                        ppr.combo{temp_idx} = [ppr.combo{temp_idx}; ioStruct.pprAvg{temp_idx}(:,t)];
                        ppr.comboStimId{temp_idx} = [ppr.comboStimId{temp_idx}; ioStruct.pprStimId{temp_idx}(:,t)]; 
                        ppr.comboPath{temp_idx} = [ppr.comboPath{temp_idx}; io.comboPath{temp_idx}];
                        ppr.comboPathCell{temp_idx} = [ppr.comboPathCell{temp_idx}; io.comboPathCell{temp_idx}];
                    end
               end
            end
        elseif string(ioStruct.redId) == "mec" && ioStruct.blueId == 0
           for ii = 1:length(ioStruct.kinetics)
                temp_idx = find(kinetics == extractBefore(ioStruct.kinetics(ii), '_2'), 1);
                io.mecR{temp_idx} = [io.mecR{temp_idx}; [ioStruct.dataRedAvg{ii}, ioStruct.dataRed{ii}(:,2)]];
                io.mecRStim{temp_idx} = [io.mecRStim{temp_idx}; ioStruct.dataAvgStim{ii}(:,1)];
                io.mecRStimId{temp_idx} = [io.mecRStimId{temp_idx}; string(ioStruct.dataRedStimId{ii})];
                io.mecRPath{temp_idx} = [io.mecRPath{temp_idx}; ioStruct.Animal, ioStruct.Cell, ioStruct.blueId, ioStruct.redId, ioStruct.clamp];
                io.mecRPathCell{temp_idx} = [io.mecRPathCell{temp_idx}; string(repmat([ioStruct.Animal, ioStruct.Cell, ioStruct.blueId, ioStruct.redId],length(ioStruct.dataAvgStim{ii}(:,1)),1))];
                if ~isempty(ioStruct.pprAvg{temp_idx})
                    t = contains(ioStruct.pprStimId{temp_idx}, 'b');
                    ppr.mecR{temp_idx} = [ppr.mecR{temp_idx}; ioStruct.pprAvg{temp_idx}(:,t)]; 
                    ppr.mecRStimId{temp_idx} = [ppr.mecRStimId{temp_idx}; ioStruct.pprStimId{temp_idx}(:,t)]; 
                    ppr.mecRPath{temp_idx} = [ppr.mecRPath{temp_idx}; io.mecRPath{temp_idx}] ;
                    ppr.mecRPathCell{temp_idx} =[ppr.mecRPathCell{temp_idx};  io.mecRPathCell{temp_idx}];   
                end
           end
            %io.mecRPath = [io.mecRPath; ioStruct.Animal, ioStruct.Cell, ioStruct.blueId, ioStruct.redId, ioStruct.clamp];
        elseif string(ioStruct.redId) == "lec" && ioStruct.blueId == 0
            for ii = 1:length(ioStruct.kinetics)
                temp_idx = find(kinetics == extractBefore(ioStruct.kinetics(ii), '_2'), 1);
                io.lecR{temp_idx} = [io.lecR{temp_idx}; [ioStruct.dataRedAvg{ii}, ioStruct.dataRed{ii}(:,2)]]; 
                io.lecRStim{temp_idx} = [io.lecRStim{temp_idx}; ioStruct.dataAvgStim{ii}(:,1)];
                io.lecRStimId{temp_idx} = [io.lecRStimId{temp_idx}; string(ioStruct.dataRedStimId{ii})];
                io.lecRPath{temp_idx} = [io.lecRPath{temp_idx}; ioStruct.Animal, ioStruct.Cell, ioStruct.blueId, ioStruct.redId, ioStruct.clamp];
                io.lecRPathCell{temp_idx} = [io.lecRPathCell{temp_idx}; string(repmat([ioStruct.Animal, ioStruct.Cell, ioStruct.blueId, ioStruct.redId],length(ioStruct.dataAvgStim{ii}(:,1)),1))];
                if ~isempty(ioStruct.pprAvg{temp_idx})
                    t = contains(ioStruct.pprStimId{temp_idx}, 'b');
                    ppr.lecR{temp_idx} = [ppr.lecR{temp_idx}; ioStruct.pprAvg{temp_idx}(:,t)]; 
                    ppr.lecRStimId{temp_idx} = [ppr.lecRStimId{temp_idx}; ioStruct.pprStimId{temp_idx}(:,t)]; 
                    ppr.lecRPath{temp_idx} = [ppr.lecRPath{temp_idx}; io.lecRPath{temp_idx}] ;
                    ppr.lecRPathCell{temp_idx} =[ppr.lecRPathCell{temp_idx};  io.lecRPathCell{temp_idx}];   
                end
            end
        end
    end

end

%% Spike Probability if CC 
f = ["lecB", "lecR", "lecRDual", "lecBDual", "mecB", "mecR", "mecRDual", "mecBDual", "combo", "comboExp"];
for i = 1:length(f)
    if contains(path_2, 'cc', 'IgnoreCase', true) && ~isempty(io.(f(i)){1})
%         if sum(find(io.(f(i)){1} > 60)) > 0
        temp_idx = find(io.(f(i)){1} > 40);
        temp_idx2 = find(io.(f(i)){1} < 40);
        temp_idx3 = find(isnan(io.(f(i)){1}));
        for ii = 1:length(length(io.(f(i))))
         spikeProb.(f(i)) = io.(f(i)){ii};
         spikeProb.(f(i))(temp_idx) = 1;
         spikeProb.(f(i))(temp_idx2) = 0;
         spikeProb.(f(i))(temp_idx3) = 0;
        end
%         end
        
        temp_field = strcat(f(i), 'Stim');
        temp_idx = find(isempty(spikeProb.(f(i))) == 0);
        %temp_field5 = strcat(f(i), 'kin'); ioAvg.(temp_field5) = kinetics(temp_idx)
        [groupIdx, categories] = grp2idx(cell2mat(io.(temp_field)(temp_idx)));
        spikeProbAvg.(f(i)) = splitapply(@(x) mean(x, 'omitnan'), spikeProb.(f(i)), groupIdx);
        spikeProbAvg.(temp_field) = categories; 

        temp_field2 = strcat(f(i), 'Std');
        spikeProbAvg.(temp_field2) = splitapply(@(x) std(x, 'omitnan'), spikeProb.(f(i)), groupIdx);
     
        temp_field3 = strcat(f(i), 'Num');
        spikeProbAvg.(temp_field3) = histc(io.(temp_field){1}, unique(io.(temp_field){1}, 'stable'));

        temp_field4 = strcat(f(i), 'Sem');
        temp_idx = 1:length(temp_idx);
        spikeProbAvg.(temp_field4) = spikeProbAvg.(temp_field2)./sqrt(spikeProbAvg.(temp_field3));
        
    end
end

if contains(path_2, 'cc', 'IgnoreCase', true)
    spikeProbAvg.mecAvg = mean([spikeProbAvg.mecRDual(:,1), spikeProbAvg.mecBDual(:,1)],2);
    spikeProbAvg.lecAvg = mean([spikeProbAvg.lecRDual(:,1), spikeProbAvg.lecBDual(:,1)],2);
    spikeProbAvg.mecSemAvg = mean([spikeProbAvg.mecRDualSem(:,1), spikeProbAvg.mecBDualSem(:,1)],2);
    spikeProbAvg.lecSemAvg = mean([spikeProbAvg.lecRDualSem(:,1), spikeProbAvg.lecBDualSem(:,1)],2);
end



%% Averaging Across Stim Values

f = ["lecB", "lecR", "lecRDual", "lecBDual", "mecB", "mecR", "mecRDual", "mecBDual", "combo", "comboExp"];
g = [1, 14];
for i = 1:length(f)
    for ii = 1:length(g)
        if sum(find(io.(f(i)){g(ii)} > 40)) > 0 && contains(path_2, 'cc','IgnoreCase',true)
            temp_idx = find(io.(f(i)){g(ii)} > 40);
            %for iii = 1:length(length(io.(f(i))))
             io.(f(i)){g(ii)}(temp_idx) = NaN;
            %end
        end
    
        if sum(find(io.(f(i)){g(ii)} <= 0)) > 0 && contains(path_2, 'cc','IgnoreCase',true)
            temp_idx = find(io.(f(i)){g(ii)} < 0.05);
            %for iii = 1:length(length(io.(f(i))))
                io.(f(i)){g(ii)}(temp_idx) = 0.1;
            %end
        end
    end
end

%% Make PPR Structure, integrated above 
% f = ["lecB", "lecR", "lecRDual", "lecBDual", "mecB", "mecR", "mecRDual", "mecBDual", "combo", "comboExp"];
% 
% for i = 1:length(f)
%     if cellfun(@(x) ~isempty(x), io.(f{i}))
%         if contains(f{i}, 'BDual')
%             ppr.(f{i}) = cellfun(@(x) x(:,1:2)./x(:,1), io.(f{i}), 'UniformOutput', false);
%             tempP = f{i}; tempP(tempP =='B') = 'R'; 
%             if tempP(1) == 'l', tempP(1) = 'm';  elseif tempP(1) == 'm', tempP(1) = 'l'; end      
%             ppr.(f{i}) = cellfun(@(z, w) [z, w], ppr.(f{i}), cellfun(@(x,y) x(:,3)./x(:,3), io.(f{i}), 'UniformOutput', false), 'UniformOutput', false);
%             ppr.(f{i}) = cellfun(@(z, w) [z, w], ppr.(f{i}), cellfun(@(x,y) x(:,4)./y(:,4), io.(f{i}), io.(tempP), 'UniformOutput', false), 'UniformOutput', false);
% 
%         elseif contains(f{i}, 'RDual')
%             ppr.(f{i}) = cellfun(@(x) x(:,1:2)./x(:,1), io.(f{i}), 'UniformOutput', false);
%             tempP = f{i}; tempP(tempP =='R') = 'B'; 
%             if tempP(1) == 'l', tempP(1) = 'm';  elseif tempP(1) == 'm', tempP(1) = 'l'; end      
%             ppr.(f{i}) = cellfun(@(z, w) [z, w], ppr.(f{i}), cellfun(@(x,y) x(:,3)./y(:,3), io.(f{i}), io.(tempP), 'UniformOutput', false), 'UniformOutput', false);
%             ppr.(f{i}) = cellfun(@(z, w) [z, w], ppr.(f{i}), cellfun(@(x,y) x(:,4)./x(:,4), io.(f{i}), 'UniformOutput', false), 'UniformOutput', false);
% 
%         else
%             ppr.(f{i}) = cellfun(@(x) x./x(:,1), io.(f{i}), 'UniformOutput', false);
%             tempIdx = cellfun(@(x) x > 50, ppr.(f{i}), 'UniformOutput', false);
%             ppr.(f{i}) = cellfun(@(x,y) setNaN(x, y), ppr.(f{i}), tempIdx,'UniformOutput', false);
%             ppr.(strcat(f{i}, 'Stim')) = io.(strcat(f{i}, 'Stim'));
%             %ppr.(strcat(f{i}, 'StimId')) = io.(strcat(f{i}, 'StimId'));
%             ppr.(strcat(f{i}, 'PathCell')) = io.(strcat(f{i}, 'PathCell'));
% 
%         end
%         tempIdx = cellfun(@(x) x > 50, ppr.(f{i}), 'UniformOutput', false);
%         ppr.(f{i}) = cellfun(@(x,y) setNaN(x, y), ppr.(f{i}), tempIdx,'UniformOutput', false);
%         ppr.(strcat(f{i}, 'Stim')) = io.(strcat(f{i}, 'Stim'));
%         %ppr.(strcat(f{i}, 'StimId')) = io.(strcat(f{i}, 'StimId'));
%         ppr.(strcat(f{i}, 'PathCell')) = io.(strcat(f{i}, 'PathCell'));
%     end
% end


for i = 1:length(f)
    if sum(cellfun(@(x) isempty(x), io.(f(i)))) < 28
        temp_field = strcat(f(i), 'Stim');
        temp_idx = find(cellfun(@(x) isempty(x), io.(f(i))) == 0);
        temp_field5 = strcat(f(i), 'kin'); 
        ioAvg.(temp_field5) = kinetics(temp_idx);

        [groupIdx, categories] = cellfun(@(x) grp2idx(x), io.(temp_field)(temp_idx),'UniformOutput',false);
        ioAvg.(f(i)) = cellfun(@(data, groups) splitapply(@(x) mean(x, 'omitnan'), data, groups), io.(f(i))(temp_idx), groupIdx, 'UniformOutput', false);
        ioAvg.(temp_field) = categories; 

        temp_field2 = strcat(f(i), 'Std');
        ioAvg.(temp_field2) = cellfun(@(data, groups) splitapply(@(x) std(x, 'omitnan'), data, groups), io.(f(i))(temp_idx), groupIdx, 'UniformOutput', false);
        
        temp_field3 = strcat(f(i), 'Num');
        ioAvg.(temp_field3) = cellfun(@(data, groups) histc(data, unique(groups, 'stable')), io.(temp_field)(temp_idx), io.(temp_field)(temp_idx), 'UniformOutput', false);

        %temp_field5 = strcat(f(i), 'kin'); ioAvg.(temp_field5) = kinetics(temp_idx);

        temp_field4 = strcat(f(i), 'Sem');
        temp_idx = 1:length(temp_idx);
        ioAvg.(temp_field4) = cellfun(@(data, groups) data./sqrt(groups), ioAvg.(temp_field2)(temp_idx), ioAvg.(temp_field3)(temp_idx), 'UniformOutput', false);

        temp_idx = find(cellfun(@(x) isempty(x), ppr.(f(i))) == 0);
        if sum(matches(fieldnames(ppr), f{i})) > 0
           temp = cellfun(@(x) x > 50, ppr.(f{i}), 'UniformOutput', false);
           ppr.(f{i}) = cellfun(@(x,y) setNaN(x, y), ppr.(f{i}), temp,'UniformOutput', false);

           pprAvg.(temp_field5) = kinetics(temp_idx);
           pprAvg.(f(i)) = cellfun(@(data, groups) splitapply(@(x) mean(x, 'omitnan'), data, groups), ppr.(f(i))(temp_idx), groupIdx(temp_idx), 'UniformOutput', false);
           pprAvg.(temp_field) = categories; 
           pprAvg.(temp_field2) = cellfun(@(data, groups) splitapply(@(x) std(x, 'omitnan'), data, groups), ppr.(f(i))(temp_idx), groupIdx(temp_idx), 'UniformOutput', false);
           pprAvg.(temp_field3) = cellfun(@(data, groups) histc(data, unique(groups, 'stable')), io.(temp_field)(temp_idx), io.(temp_field)(temp_idx), 'UniformOutput', false);
           temp_idx = 1:length(temp_idx);
           pprAvg.(temp_field4) = cellfun(@(data, groups) data./sqrt(groups), pprAvg.(temp_field2)(temp_idx), pprAvg.(temp_field3)(temp_idx), 'UniformOutput', false);
        end
%         n = cellfun(@(data, groups) histc(@(x) std(x, 'omitnan'), data, groups), io.(temp_field(i))(temp_idx), unique(io.(temp_field(i))(temp_idx)), 'UniformOutput', false);
%         
%         uStim = unique(io.(temp_field){ii});
%         for iii = 1:length(uStim)
%            nu_vec(iii) = sum(io.(temp_field){ii} == uStim(iii))
%         end
    else
      fprintf([mouse, id])  
    end
end

%Extra Stuff  - PPR





%% Variables for Plotting - Manual Add

ioAvg.mWBlue = [0.0048; 0.0127; 0.021; 0.0291; 0.0375; 0.0457; 0.0458; 0.0623; ...
    0.0704; 0.0787; 0.158; 0.235; 0.309; 0.381; 0.45; 0.518; 0.584; 0.648; 0.711; ... 
    1.283; 1.78; 2.22; 2.63; 2.99; 3.32; 3.6; 3.84; 4.06];
ioAvg.mWBlue(:,2) = (ioAvg.mWBlue./max(ioAvg.mWBlue))*100;
ioAvg.mWBlue(:,3) = [(0.1:0.1:1), (2:1:10), (20:10:100)]';

ioAvg.mWRed = [0.002; 0.0059; 0.01; 0.015; 0.0198; 0.0247; 0.0247; 0.0348; 0.04; 0.0451;...
    0.0984; 0.153; 0.208; 0.263; 0.318; 0.374; 0.429; 0.484; 0.54; 1.09; 1.63;...
    2.16; 2.66;3.1; 3.47;3.79;4.07;4.33];
ioAvg.mWRed(:,2) = (ioAvg.mWRed./max(ioAvg.mWRed))*100;
ioAvg.mWRed(:,3) = [(0.1:0.1:1), (2:1:10), (20:10:100)]';


%% Manual Plotting
if plt

if contains(path_2, 'vc', 'IgnoreCase',true)
    for i = 1:length(ioAvg.lecB)
        figure;
        errorbar(ioAvg.mWBlue(10:end,2), ioAvg.lecB{i}(:,1), ioAvg.lecBSem{i}(:,1), '.-','MarkerSize', 30, 'Color', '#289e1e',... 
            'LineWidth',2 ); hold on;
        errorbar(ioAvg.mWBlue(10:end,2), ioAvg.mecB{i}(:,1), ioAvg.mecBSem{i}(:,1), 'r.-', 'MarkerSize', 30,'LineWidth',2);
        xline(18, '--', 'LineWidth',3);
        legend('LEC(14)', 'MEC(50)', 'FontSize', 16); ylabel(extractBefore(kinetics(i), '_'), 'FontSize', 16); xlabel('LED Intensity (%)', 'FontSize', 16);
        set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3); title(['I/O', kinetics(i), 'Single Labelled Chronos VC']);

        figure;
        errorbar(ioAvg.mWBlue(10:end,2), ioAvg.lecB{i}(:,2), ioAvg.lecBSem{i}(:,2), '.-','MarkerSize', 30, 'Color', '#289e1e','LineWidth',2); hold on;
        errorbar(ioAvg.mWBlue(10:end,2), ioAvg.mecB{i}(:,2), ioAvg.mecBSem{i}(:,2), 'r.-', 'MarkerSize', 30,'LineWidth',2);
        xline(18, '--', 'LineWidth',3);
        legend('LEC(14)', 'MEC(50)', 'FontSize', 16); ylabel(extractBefore(kinetics(i), '_'), 'FontSize', 16); xlabel('LED Intensity (%)', 'FontSize', 16);
        set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3); title(['I/O', kinetics(i), 'Single Labelled Chronos VC Pulse 2']);
        
    end
    for i = 1:length(ioAvg.lecBDual)
        figure;
        errorbar(ioAvg.mWBlue(10:end,2), ioAvg.lecBDual{i}(:,1), ioAvg.lecBDualSem{i}(:,1), '.-','MarkerSize', 30, 'Color', '#289e1e','LineWidth',2); hold on;
        errorbar(ioAvg.mWRed(10:end,2), ioAvg.mecRDual{i}(:,1), ioAvg.mecRDualSem{i}(:,1), 'r.-', 'MarkerSize', 30,'LineWidth',2);
        xline(18, '--', 'LineWidth',3);
        legend('LECB(8)', 'MECR(8)', 'FontSize', 16); ylabel(extractBefore(ioAvg.lecBDualkin(i), '_'), 'FontSize', 16); xlabel('LED Intensity (%)', 'FontSize', 16);
        set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3); title(['I/O', io.lecBDual{2}(i), 'Dual Labelled Chronos VC Pulse 1']); 
    end
    for i = 1:length(ioAvg.mecBDual)
        figure;
        errorbar(ioAvg.mWRed(10:end,2), ioAvg.lecRDual{i}(:,1), ioAvg.lecRDualSem{i}(:,1), '.-','MarkerSize', 30, 'Color', '#289e1e','LineWidth',2); hold on;
        errorbar(ioAvg.mWBlue(10:end,2), ioAvg.mecBDual{i}(:,1), ioAvg.mecBDualSem{i}(:,1), 'r.-', 'MarkerSize', 30,'LineWidth',2);
        xline(18, '--', 'LineWidth',3);
        legend('LECR(4)', 'MECB(4)', 'FontSize', 16); ylabel(extractBefore(ioAvg.lecBDualkin(i), '_'), 'FontSize', 16); xlabel('LED Intensity (%)', 'FontSize', 16);
        set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3); title(['I/O', io.lecBDual{2}(i), 'Dual Labelled Chronos VC Pulse 1']); 
    end
end

%%

% if contains(path, 'vc', 'IgnoreCase', true)
%     
%     figure;
%     errorbar(str2double(ioAvg.lecBStim{1}), ioAvg.lecB{1}(:,1), ioAvg.lecBSem{1}(:,1), '.-','MarkerSize', 30, 'Color', '#289e1e'); hold on;
%     errorbar( str2double(ioAvg.mecBStim{1}), ioAvg.mecB{1}(:,1), ioAvg.mecBSem{1}(:,1), 'r.-', 'MarkerSize', 30);
%     legend('LEC (14)', 'MEC(50)', 'FontSize', 16); ylabel('Amplitude (mV)', 'FontSize', 16); xlabel('LED Intensity (%)', 'FontSize', 16);
%     set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
%     
%     figure;
%     errorbar(str2double(ioAvg.lecBStim{2}), ioAvg.lecB{2}(:,1), ioAvg.lecBSem{2}(:,1), '.-','MarkerSize', 30, 'Color', '#289e1e'); hold on;
%     errorbar( str2double(ioAvg.mecBStim{2}), ioAvg.mecB{2}(:,1), ioAvg.mecBSem{2}(:,1), 'r.-', 'MarkerSize', 30); 
%     legend('LEC (14)', 'MEC(50)', 'FontSize', 16); ylabel('Amplitude (mV)', 'FontSize', 16); xlabel('LED Intensity (%)', 'FontSize', 16);
%     set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
%     
%     figure;
%     errorbar(str2double(ioAvg.lecBStim{2}), ioAvg.lecB{2}(:,1), ioAvg.lecBSem{2}(:,1), '.-','MarkerSize', 30, 'Color', '#289e1e'); hold on;
%     errorbar( str2double(ioAvg.mecBStim{2}), ioAvg.mecB{2}(:,1), ioAvg.mecBSem{2}(:,1), 'r.-', 'MarkerSize', 30); 
%     errorbar(str2double(ioAvg.lecBStim{1}), ioAvg.lecB{1}(:,1), ioAvg.lecBSem{1}(:,1), '.-','MarkerSize', 30, 'Color', '#289e1e');
%     errorbar( str2double(ioAvg.mecBStim{1}), ioAvg.mecB{1}(:,1), ioAvg.mecBSem{1}(:,1), 'r.-', 'MarkerSize', 30);
%     legend('LEC (17)', 'MEC(61)', 'FontSize', 16); ylabel('Amplitude (mV)', 'FontSize', 16); xlabel('LED Intensity (%)', 'FontSize', 16);
%     set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
%     
%     %
%     
%     figure;
%     errorbar(str2double(ioAvg.lecBStim{3}), ioAvg.lecB{3}(:,1), ioAvg.lecBSem{3}(:,1), '.-','MarkerSize', 30, 'Color', '#289e1e'); hold on;
%     errorbar( str2double(ioAvg.mecBStim{3}), ioAvg.mecB{3}(:,1), ioAvg.mecBSem{3}(:,1), 'r.-', 'MarkerSize', 30);
%     legend('LEC (17)', 'MEC(61)', 'FontSize', 16); ylabel('AUC 10-90', 'FontSize', 16); xlabel('LED Intensity (%)', 'FontSize', 16);
%     set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
%     
%     figure;
%     errorbar(str2double(ioAvg.lecBStim{4}), ioAvg.lecB{4}(:,1), ioAvg.lecBSem{4}(:,1), '.-','MarkerSize', 30, 'Color', '#289e1e'); hold on;
%     errorbar( str2double(ioAvg.mecBStim{4}), ioAvg.mecB{4}(:,1), ioAvg.mecBSem{4}(:,1), 'r.-', 'MarkerSize', 30); 
%     legend('LEC (17)', 'MEC(61)', 'FontSize', 16); ylabel('AUC 10-90', 'FontSize', 16); xlabel('LED Intensity (%)', 'FontSize', 16);
%     set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
%     
%     %
%     
%     figure;
%     errorbar(str2double(ioAvg.lecBStim{7}), ioAvg.lecB{7}(:,1), ioAvg.lecBSem{7}(:,1), '.-','MarkerSize', 30, 'Color', '#289e1e'); hold on;
%     errorbar( str2double(ioAvg.mecBStim{7}), ioAvg.mecB{7}(:,1), ioAvg.mecBSem{7}(:,1), 'r.-', 'MarkerSize', 30);
%     legend('LEC (17)', 'MEC(61)', 'FontSize', 16); ylabel('dc 10-90', 'FontSize', 16); xlabel('LED Intensity (%)', 'FontSize', 16);
%     set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
%     
%     figure;
%     errorbar(str2double(ioAvg.lecBStim{8}), ioAvg.lecB{8}(:,1), ioAvg.lecBSem{8}(:,1), '.-','MarkerSize', 30, 'Color', '#289e1e'); hold on;
%     errorbar( str2double(ioAvg.mecBStim{8}), ioAvg.mecB{8}(:,1), ioAvg.mecBSem{8}(:,1), 'r.-', 'MarkerSize', 30); 
%     legend('LEC (17)', 'MEC(61)', 'FontSize', 16); ylabel('dc 10-90', 'FontSize', 16); xlabel('LED Intensity (%)', 'FontSize', 16);
%     set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
%     
%     % 
%     
%     figure;
%     errorbar(str2double(ioAvg.lecBStim{9}), ioAvg.lecB{9}(:,1), ioAvg.lecBSem{9}(:,1), '.-','MarkerSize', 30, 'Color', '#289e1e'); hold on;
%     errorbar( str2double(ioAvg.mecBStim{9}), ioAvg.mecB{9}(:,1), ioAvg.mecBSem{9}(:,1), 'r.-', 'MarkerSize', 30);
%     legend('LEC (17)', 'MEC(61)', 'FontSize', 16); ylabel('dc 20-80', 'FontSize', 16); xlabel('LED Intensity (%)', 'FontSize', 16);
%     set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
%     
%     figure;
%     errorbar(str2double(ioAvg.lecBStim{10}), ioAvg.lecB{10}(:,1), ioAvg.lecBSem{10}(:,1), '.-','MarkerSize', 30, 'Color', '#289e1e'); hold on;
%     errorbar( str2double(ioAvg.mecBStim{10}), ioAvg.mecB{10}(:,1), ioAvg.mecBSem{10}(:,1), 'r.-', 'MarkerSize', 30); 
%     legend('LEC (17)', 'MEC(61)', 'FontSize', 16); ylabel('dc 20-80', 'FontSize', 16); xlabel('LED Intensity (%)', 'FontSize', 16);
%     set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
%     
%     %
%     
%     figure;
%     errorbar(str2double(ioAvg.lecBStim{11}), ioAvg.lecB{11}(:,1), ioAvg.lecBSem{11}(:,1), '.-','MarkerSize', 30, 'Color', '#289e1e'); hold on;
%     errorbar( str2double(ioAvg.mecBStim{11}), ioAvg.mecB{11}(:,1), ioAvg.mecBSem{11}(:,1), 'r.-', 'MarkerSize', 30);
%     legend('LEC (17)', 'MEC(61)', 'FontSize', 16); ylabel('hw', 'FontSize', 16); xlabel('LED Intensity (%)', 'FontSize', 16);
%     set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
%     
%     figure;
%     errorbar(str2double(ioAvg.lecBStim{12}), ioAvg.lecB{12}(:,1), ioAvg.lecBSem{12}(:,1), '.-','MarkerSize', 30, 'Color', '#289e1e'); hold on;
%     errorbar( str2double(ioAvg.mecBStim{12}), ioAvg.mecB{12}(:,1), ioAvg.mecBSem{12}(:,1), 'r.-', 'MarkerSize', 30); 
%     legend('LEC (17)', 'MEC(61)', 'FontSize', 16); ylabel('hw', 'FontSize', 16); xlabel('LED Intensity (%)', 'FontSize', 16);
%     set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
%     
%     %
%     
%     figure;
%     errorbar(str2double(ioAvg.lecBStim{23}), ioAvg.lecB{23}(:,1), ioAvg.lecBSem{23}(:,1), '.-','MarkerSize', 30, 'Color', '#289e1e'); hold on;
%     errorbar( str2double(ioAvg.mecBStim{23}), ioAvg.mecB{23}(:,1), ioAvg.mecBSem{23}(:,1), 'r.-', 'MarkerSize', 30);
%     legend('LEC (17)', 'MEC(61)', 'FontSize', 16); ylabel('rs 10-90', 'FontSize', 16); xlabel('LED Intensity (%)', 'FontSize', 16);
%     set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
%     
%     figure;
%     errorbar(str2double(ioAvg.lecBStim{24}), ioAvg.lecB{24}(:,1), ioAvg.lecBSem{24}(:,1), '.-','MarkerSize', 30, 'Color', '#289e1e'); hold on;
%     errorbar( str2double(ioAvg.mecBStim{24}), ioAvg.mecB{24}(:,1), ioAvg.mecBSem{24}(:,1), 'r.-', 'MarkerSize', 30); 
%     legend('LEC (17)', 'MEC(61)', 'FontSize', 16); ylabel('rs 10-90', 'FontSize', 16); xlabel('LED Intensity (%)', 'FontSize', 16);
%     set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
%     
%     %
%     
%     figure;
%     errorbar(str2double(ioAvg.lecBStim{25}), ioAvg.lecB{25}(:,1), ioAvg.lecBSem{25}(:,1), '.-','MarkerSize', 30, 'Color', '#289e1e'); hold on;
%     errorbar( str2double(ioAvg.mecBStim{25}), ioAvg.mecB{25}(:,1), ioAvg.mecBSem{25}(:,1), 'r.-', 'MarkerSize', 30);
%     legend('LEC (17)', 'MEC(61)', 'FontSize', 16); ylabel('rs 20-80', 'FontSize', 16); xlabel('LED Intensity (%)', 'FontSize', 16);
%     set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
%     
%     figure;
%     errorbar(str2double(ioAvg.lecBStim{26}), ioAvg.lecB{26}(:,1), ioAvg.lecBSem{26}(:,1), '.-','MarkerSize', 30, 'Color', '#289e1e'); hold on;
%     errorbar( str2double(ioAvg.mecBStim{26}), ioAvg.mecB{26}(:,1), ioAvg.mecBSem{26}(:,1), 'r.-', 'MarkerSize', 30); 
%     legend('LEC (17)', 'MEC(61)', 'FontSize', 16); ylabel('rs 20-80', 'FontSize', 16); xlabel('LED Intensity (%)', 'FontSize', 16);
%     set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
%     
%     %
%     
%     figure;
%     errorbar(str2double(ioAvg.lecBStim{27}), ioAvg.lecB{27}(:,2), ioAvg.lecBSem{27}(:,2), '.-','MarkerSize', 30, 'Color', '#289e1e'); hold on;
%     errorbar( str2double(ioAvg.mecBStim{27}), ioAvg.mecB{27}(:,2), ioAvg.mecBSem{27}(:,2), 'r.-', 'MarkerSize', 30);
%     legend('LEC (17)', 'MEC(61)', 'FontSize', 16); ylabel('sum', 'FontSize', 16); xlabel('LED Intensity (%)', 'FontSize', 16);
%     set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
%     
%     figure;
%     errorbar(str2double(ioAvg.lecBStim{28}), ioAvg.lecB{28}(:,2), ioAvg.lecBSem{28}(:,2), '.-','MarkerSize', 30, 'Color', '#289e1e'); hold on;
%     errorbar( str2double(ioAvg.mecBStim{28}), ioAvg.mecB{28}(:,2), ioAvg.mecBSem{28}(:,2), 'r.-', 'MarkerSize', 30); 
%     legend('LEC (17)', 'MEC(61)', 'FontSize', 16); ylabel('sum', 'FontSize', 16); xlabel('LED Intensity (%)', 'FontSize', 16);
%     set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
%     
%     %
%     
%     figure;
%     errorbar(str2double(ioAvg.lecBStim{13}), ioAvg.lecB{13}(:,1), ioAvg.lecBSem{13}(:,1), '.-','MarkerSize', 30, 'Color', '#289e1e'); hold on;
%     errorbar( str2double(ioAvg.mecBStim{13}), ioAvg.mecB{13}(:,1), ioAvg.mecBSem{13}(:,1), 'r.-', 'MarkerSize', 30);
%     legend('LEC (17)', 'MEC(61)', 'FontSize', 16); ylabel('latency', 'FontSize', 16); xlabel('LED Intensity (%)', 'FontSize', 16);
%     set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
%     
%     figure;
%     errorbar(str2double(ioAvg.lecBStim{14}), ioAvg.lecB{14}(:,1), ioAvg.lecBSem{14}(:,1), '.-','MarkerSize', 30, 'Color', '#289e1e'); hold on;
%     errorbar( str2double(ioAvg.mecBStim{14}), ioAvg.mecB{14}(:,1), ioAvg.mecBSem{14}(:,1), 'r.-', 'MarkerSize', 30); 
%     legend('LEC (17)', 'MEC(61)', 'FontSize', 16); ylabel('latency', 'FontSize', 16); xlabel('LED Intensity (%)', 'FontSize', 16);
%     set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
%     
%     %
%     
%     figure;
%     errorbar(str2double(ioAvg.lecBStim{15}), ioAvg.lecB{15}(:,1), ioAvg.lecBSem{15}(:,1), '.-','MarkerSize', 30, 'Color', '#289e1e'); hold on;
%     errorbar( str2double(ioAvg.mecBStim{15}), ioAvg.mecB{15}(:,1), ioAvg.mecBSem{15}(:,1), 'r.-', 'MarkerSize', 30);
%     legend('LEC (17)', 'MEC(61)', 'FontSize', 16); ylabel('latency', 'FontSize', 16); xlabel('LED Intensity (%)', 'FontSize', 16);
%     set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
%     
%     figure;
%     errorbar(str2double(ioAvg.lecBStim{16}), ioAvg.lecB{16}(:,1), ioAvg.lecBSem{16}(:,1), '.-','MarkerSize', 30, 'Color', '#289e1e'); hold on;
%     errorbar( str2double(ioAvg.mecBStim{16}), ioAvg.mecB{16}(:,1), ioAvg.mecBSem{16}(:,1), 'r.-', 'MarkerSize', 30); 
%     legend('LEC (17)', 'MEC(61)', 'FontSize', 16); ylabel('latency', 'FontSize', 16); xlabel('LED Intensity (%)', 'FontSize', 16);
%     set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
% 
%     %
%     figure;
%     errorbar(str2double(ioAvg.lecBStim{15}), ioAvg.lecB{15}(:,1), ioAvg.lecBSem{15}(:,1), '.-','MarkerSize', 30, 'Color', '#289e1e'); hold on;
%     errorbar( str2double(ioAvg.mecBStim{15}), ioAvg.mecB{15}(:,1), ioAvg.mecBSem{15}(:,1), 'r.-', 'MarkerSize', 30);
%     legend('LEC (17)', 'MEC(61)', 'FontSize', 16); ylabel('latency', 'FontSize', 16); xlabel('LED Intensity (%)', 'FontSize', 16);
%     set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
% end


%%
if contains(path_2, 'cc','IgnoreCase',true)
    
    figure;
    errorbar(ioAvg.mWBlue(10:end, 2), ioAvg.lecBDual{1}(:,1), ioAvg.lecBDualSem{1}(:,1), '.-','MarkerSize', 30, 'Color', '#289e1e', ...
        'LineWidth', 2); hold on;
    errorbar(ioAvg.mWRed(10:end, 2), ioAvg.mecRDual{1}(:,1), ioAvg.mecRDualSem{1}(:,1), 'r.-', 'MarkerSize', 30, 'LineWidth', 2);
    errorbar(ioAvg.mWBlue(10:end, 2), ioAvg.combo{1}(:,1), ioAvg.comboSem{1}(:,1), 'm.-', 'MarkerSize', 30, 'LineWidth', 2);
    legend('LECB (30)', 'MECR (30)','EC (47)', 'FontSize', 16); ylabel('Amplitude', 'FontSize', 16); xlabel('LED Intensity (%)', 'FontSize', 16);
    set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);

    figure; 
    errorbar(ioAvg.mWBlue(10:19, 2), ioAvg.lecBDual{1}(1:10,1), ioAvg.lecBDualSem{1}(1:10,1), '.-','MarkerSize', 30, 'Color', '#289e1e', ...
        'LineWidth', 2); hold on;
    errorbar(ioAvg.mWRed(10:19, 2), ioAvg.mecRDual{1}(1:10,1), ioAvg.mecRDualSem{1}(1:10,1), 'r.-', 'MarkerSize', 30, 'LineWidth', 2);
    errorbar(ioAvg.mWBlue(10:19, 2), ioAvg.combo{1}(1:10,1), ioAvg.comboSem{1}(1:10,1), 'm.-', 'MarkerSize', 30, 'LineWidth', 2);
    legend('LECB (30)', 'MECR (30)', 'EC (47)', 'FontSize', 16); ylabel('Amplitude (mV)', 'FontSize', 16); xlabel('LED Intensity (%)', 'FontSize', 16);
    set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);

    figure;
    errorbar(ioAvg.mWRed(10:end, 2), ioAvg.lecRDual{1}(:,1), ioAvg.lecRDualSem{1}(:,1), '.-','MarkerSize', 30, 'Color', '#289e1e', ...
        'LineWidth', 2); hold on;
    errorbar(ioAvg.mWBlue(10:end, 2), ioAvg.mecBDual{1}(:,1), ioAvg.mecBDualSem{1}(:,1), 'r.-', 'MarkerSize', 30, 'LineWidth', 2);
    errorbar(ioAvg.mWBlue(10:end, 2), ioAvg.combo{1}(:,1), ioAvg.comboSem{1}(:,1), 'm.-', 'MarkerSize', 30, 'LineWidth', 2);
    legend('LECR (17)', 'MECB (17)', 'EC (47)', 'FontSize', 16); ylabel('Amplitude (mV)', 'FontSize', 16); xlabel('LED Intensity (%)', 'FontSize', 16);
    set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);

    figure; 
    errorbar(ioAvg.mWRed(10:19, 2), ioAvg.lecRDual{1}(1:10,1), ioAvg.lecRDualSem{1}(1:10,1), '.-','MarkerSize', 30, 'Color', '#289e1e', ...
        'LineWidth', 2); hold on;
    errorbar(ioAvg.mWBlue(10:19, 2), ioAvg.mecBDual{1}(1:10,1), ioAvg.mecBDualSem{1}(1:10,1), 'r.-', 'MarkerSize', 30, 'LineWidth', 2);
    errorbar(ioAvg.mWBlue(10:19, 2), ioAvg.combo{1}(1:10,1), ioAvg.comboSem{1}(1:10,1), 'm.-', 'MarkerSize', 30, 'LineWidth', 2);
    legend('LECR (17)', 'MECB (17)', 'EC (47)', 'FontSize', 16); ylabel('Amplitude (mV)', 'FontSize', 16); xlabel('LED Intensity (%)', 'FontSize', 16);
    set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);

    figure; 
    errorbar(ioAvg.mWRed(10:end, 2), spikeProbAvg.lecRDual(:,1), spikeProbAvg.lecRDualSem(:,1), '.-','MarkerSize', 30, 'Color', '#289e1e', ...
        'LineWidth', 2); hold on;
    errorbar(ioAvg.mWBlue(10:end, 2), spikeProbAvg.mecBDual(:,1), spikeProbAvg.mecBDualSem(:,1), 'r.-', 'MarkerSize', 30, 'LineWidth', 2);
    errorbar(ioAvg.mWBlue(10:end, 2), spikeProbAvg.combo(:,1), spikeProbAvg.comboSem(:,1), 'm.-', 'MarkerSize', 30, 'LineWidth', 2);
    legend('LECR (17)', 'MECB (17)', 'EC (47)', 'FontSize', 16); ylabel('Spike Probability', 'FontSize', 16); xlabel('LED Intensity (%)', 'FontSize', 16);
    set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);

    figure; 
    errorbar(ioAvg.mWBlue(10:end, 2), spikeProbAvg.lecBDual(:,1), spikeProbAvg.lecBDualSem(:,1), '.-','MarkerSize', 30, 'Color', '#289e1e', ...
        'LineWidth', 2); hold on;
    errorbar(ioAvg.mWRed(10:end, 2), spikeProbAvg.mecRDual(:,1), spikeProbAvg.mecRDualSem(:,1), 'r.-', 'MarkerSize', 30, 'LineWidth', 2);
    errorbar(ioAvg.mWBlue(10:end, 2), spikeProbAvg.combo(:,1), spikeProbAvg.comboSem(:,1), 'm.-', 'MarkerSize', 30, 'LineWidth', 2);
    legend('LECB (30)', 'MECR (30)', 'EC (47)', 'FontSize', 16); ylabel('Spike Probability', 'FontSize', 16); xlabel('LED Intensity (%)', 'FontSize', 16);
    set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);

    figure; 
    errorbar(ioAvg.mWBlue(10:end, 2), spikeProbAvg.lecAvg(:,1), spikeProbAvg.lecSemAvg(:,1), '.-','MarkerSize', 30, 'Color', '#289e1e', ...
        'LineWidth', 2); hold on;
    errorbar(ioAvg.mWBlue(10:end, 2), spikeProbAvg.mecAvg(:,1), spikeProbAvg.mecSemAvg(:,1), 'r.-', 'MarkerSize', 30, 'LineWidth', 2);
    errorbar(ioAvg.mWBlue(10:end, 2), spikeProbAvg.combo(:,1), spikeProbAvg.comboSem(:,1), 'm.-', 'MarkerSize', 30, 'LineWidth', 2);
    legend('LEC (47)', 'MEC (47)', 'EC (47)', 'FontSize', 16); ylabel('Spike Probability', 'FontSize', 16); xlabel('LED Intensity (%)', 'FontSize', 16);
    set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
    
end


for i = 1:length(ioAvg.lecB)
        figure;
        errorbar(ioAvg.mWBlue(10:end,2), ioAvg.lecB{i}(:,1), ioAvg.lecBSem{i}(:,1), '.-','MarkerSize', 30, 'Color', '#289e1e', 'LineWidth',2); hold on;
        errorbar(ioAvg.mWBlue(10:end,2), ioAvg.mecB{i}(:,1), ioAvg.mecBSem{i}(:,1), 'r.-', 'MarkerSize', 30,'LineWidth',2);
        xline(18, '--', 'LineWidth',3);
        legend('LEC(14)', 'MEC(56)', 'FontSize', 16); ylabel(extractBefore(kinetics(i), '_'), 'FontSize', 16); xlabel('LED Intensity (%)', 'FontSize', 16);
        set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3); title([kinetics(i), 'Single Labelled Chronos VC']);
end
end 

if contains(path_2, 'cc','IgnoreCase',true)
    save(['Compile_Kinetics'], 'io', 'ioAvg', 'kinetics', 'spikeProb', 'spikeProbAvg', 'ppr', 'pprAvg');
end

if contains(path_2, 'vc','IgnoreCase',true)
    save(['Compile_Kinetics'], 'io', 'ioAvg', 'kinetics','ppr', 'pprAvg');
end