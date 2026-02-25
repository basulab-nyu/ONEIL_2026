%% Graph Individual Mouse Kinetics - Updated 04/24/2025

close all; clear all;

mouse = 'KO148';
date = '20230809';
id = 'S3C3';

path = ['/Volumes/KO_Portable/MEC_LEC_CA3/data_converted/', date, '_', mouse,'/',date,'_',mouse,'_',id,'/']; cd(path);
path2 = ['/Volumes/KO_Portable/MEC_LEC_CA3/data_converted/', date, '_', mouse,'/',date,'_',mouse,'_',id,'/io'];

his = 0; %Do you want to plot histogram of all responses
plt = 0;

files = dir([path, 'io/*.csv']);
load([path, mouse, id, '_dataStruct'], 'dataStruct');
pathway = readcell('pathway'); 
kinetics = ["sum", "rs", "pkT", "lagB", "lagR", "dc", "hw", "amp", "ar19", "bsT", "ost"];

%% Do not run the rest of the script if the cell is labelled do not include in the pathway file
if sum(pathway{string(pathway)=='include',2}) < 1
    error('do not include cell')
    fprintf(mouse)
    fprintf(id)
    return
end

%% Get Kinetics Files - list of files wanted can be edited using variable kinetics, need to edit if want something else (might be able to take this out)
files_name = [];
for i = 1:length(files)
    if contains(files(i).name,kinetics) && ~startsWith(files(i).name, '._')
       files_name = [files_name; string(files(i).name)];
    end
end


%% Find Stim Information - restrict to kinetics files only
cd(path2);
for i = 1:length(files_name)
    io = readtable(files_name(i)); io.Properties.VariableNames{1} = 'start';

    if all(contains(io.start, strcat(extractBetween(date, 5,6), '_',extractAfter(date, 6), '_0')))
        t = strcat(extractBetween(date, 5,6), '_',extractAfter(date, 6), '_0');
        file_num = unique(cellfun(@str2num, extractBetween(io.start,t, '_1')),'stable'); 
        call_num = file_num + 1-(cellfun(@str2num,  extractBetween(dataStruct{1}.name,t, '.abf'))); 

    elseif ~all(contains(io.start, date))
        t = strcat(extractBetween(string(io.start(1)), 'w', '_0'), '_0');
        file_num = unique(cellfun(@str2num, extractBetween(io.start,t, '_1')),'stable'); 
        call_num = file_num + 1-(cellfun(@str2num,  extractBetween(dataStruct{1}.name,t, '.abf'))); 

    elseif all(contains(io.start, date))
         t = strcat(date,'_0');
         file_num = unique(cellfun(@str2num, extractBetween(io.start,t, '_1')),'stable'); 
         call_num = file_num + 1-(cellfun(@str2num,  extractBetween(dataStruct{1}.name,t, '.abf'))); 
    end
    %file_num = unique(cellfun(@str2num, extractBetween(io.start,12, 14)),'stable', 'UniformOutput', false); %need to change for older files!!
    %normally _00 but change for long file numbers

    %call_num = file_num + 1-(cellfun(@str2num,  extractBetween(dataStruct{1}.name,'_00', '.abf'))); 
    stimIntensity{i} = [];
    for j = 1:length(file_num)
        if round(max(dataStruct{call_num(j)}.chan{3}),1) == round(max(dataStruct{call_num(j)}.chan{2}),1)
            sI = round(max(dataStruct{call_num(j)}.chan{2}),1);
            sI_temp = stimIntensity{i};
            stimIntensity{i} = [sI_temp sI];
        else
            fprintf('stims at different mV, need to verify')
        end   
    end
    
end

%% Extract Data from Table and Average the Proper Trials
% output is ioStruct with the following fields:
% Data, Sweep (sweep number), stimIntensity, Name (what stim are you
% looking at ex: blue or red or blue blue etc. The previous fields have
% n = length(files_name) cells in them
% Animal, Date, Cell

%% Extract Data from all desired files
idx_drug = ["cno", "drug", "apv", "sr"];
% for i = 1:length(files_name)
%     io = readtable(files_name(i));  io.Properties.VariableNames{1} = ['start'];
%     ioStruct.stimIdAll{i} = extractAfter(extractBetween(io.Properties.VariableNames(:,2:end), '_', '___'), '_');
% 
% 
%     if any(contains(io.start, idx_drug))
%         tempIdx = contains(io.start, idx_drug);
%         ioStruct.drugSweep{i} = extractAfter(io.start(tempIdx), '_0');
%         ioStruct.drugStim{i} = stimIntensity{i}(tempIdx);
%     end
% 
%     tempIdx = ~contains(io.start, idx_drug);
%     ioStruct.Data{i} = io{tempIdx,2:end}; 
%     ioStruct.stimIntensity{i} = stimIntensity{i}(tempIdx);
% end
% ioStruct.Animal = mouse; ioStruct.Date = date; ioStruct.Cell = id;
% ioStruct.kinetics = files_name;
% ioStruct.blueId = pathway{string(pathway)=='blue',2};
% ioStruct.redId = pathway{string(pathway)=='red',2};
% ioStruct.in = pathway{string(pathway)=='in',2};
% ioStruct.clamp = pathway{string(pathway)=='clamp',2};
% 
% %%
% if contains(ioStruct.clamp, 'cc', 'IgnoreCase', true)
%     temp_idx = ioStruct.Data{contains(files_name, 'sum')} <= 0;
%     ioStruct.Data = cellfun(@(x) setNaN(x, temp_idx), ioStruct.Data, 'UniformOutput', false);
% 
%     temp_idx = ioStruct.Data{contains(files_name, 'amp')} <= 0;
%     ioStruct.Data = cellfun(@(x) setNaN(x, temp_idx), ioStruct.Data, 'UniformOutput', false);
% elseif contains(ioStruct.clamp, 'vc', 'IgnoreCase', true) 
%     temp_idx = ioStruct.Data{contains(files_name, 'amp_io10')} <= 0;
%     temp_idx2 = find(contains(files_name, '_io10'));
%     for i = 1:length(temp_idx2)
%         ioStruct.Data{temp_idx2(i)}(temp_idx) = NaN;
%     end
%     temp_idx = ioStruct.Data{contains(files_name, 'sum_io80')} >= 0;
%     for i = 1:length(temp_idx2)
%         ioStruct.Data{temp_idx2(i)}(temp_idx) = NaN;
%     end
% 
%     temp_idx = ioStruct.Data{contains(files_name, 'amp_io80')} >= 0;
%     temp_idx2 = find(contains(files_name, '_io80'));
%     for i = 1:length(temp_idx2)
%         ioStruct.Data{temp_idx2(i)}(temp_idx) = NaN;
%     end
% 
%     temp_idx = ioStruct.Data{contains(files_name, 'sum_io80')} >= 0;
%     for i = 1:length(temp_idx2)
%         ioStruct.Data{temp_idx2(i)}(temp_idx) = NaN;
%     end
% end
% 
% [u,~, g] = cellfun(@(x) unique(x,'stable'), ioStruct.stimIntensity, 'UniformOutput',false);
% ioStruct.dataAvg = cellfun(@(data, groups) splitapply(@(x) mean(x,1, 'omitnan'), data, groups), ioStruct.Data, g,'UniformOutput', false);
% ioStruct.dataAvgStim = u;


for i = 1:length(files_name)
   
    cd(path2);
    io = readtable(files_name(i)); io.Properties.VariableNames{1} = ['start'];
    partStr = {};
    %ioStruct.stimIdAll{i} = extractAfter(extractBetween(io.Properties.VariableNames(:,2:end), '_', '___'), '_');
    ioStruct.stimIdAll{i} = extractAfter(extractBetween(io.Properties.VariableNames(:,2:end), '_', '___'), '_');
    
    for j = 1:length(io.start)
        %% sort out files that are w/ ttx or drug (CNO)
        if contains(io.start(j),'ttx', 'IgnoreCase', true) || contains(io.start(j),'drug', 'IgnoreCase', true) || contains(io.start(j),'sr', 'IgnoreCase', true) || contains(io.start(j),'apv', 'IgnoreCase', true) 
            %strcmp('ttx',extractBetween(amp_io.start,strlength(amp_io.start)-2, strlength(amp_io.start)))
            if contains(io.start(j), strcat(extractBetween(date, 5,6), '_',extractAfter(date, 6), '_0'))
                %partStr(j) = extractBetween(io.start(j),strlength(io.start(j))-12, strlength(io.start(j)));
                partStr(j) = extractAfter(io.start(j), strcat(extractBetween(date, 5,6), '_',extractAfter(date, 6), '_0'))
            elseif ~contains(io.start(j), date)
                partStr(j) = extractAfter(io.start(j), '_0');
            else
                partStr(j) = extractAfter(io.start(j), strcat(date, '_0'));
            end
        else
            if contains(io.start(j), strcat(extractBetween(date, 5,6), '_',extractAfter(date, 6), '_0'))
                %partStr(j) = extractBetween(io.start(j),strlength(io.start(j))-12, strlength(io.start(j)));
                partStr(j) = extractAfter(io.start(j), strcat(extractBetween(date, 5,6), '_',extractAfter(date, 6), '_0'));
            elseif ~contains(io.start(j), date)
                partStr(j) = extractAfter(io.start(j), '_0');
            else
                partStr(j) = extractAfter(io.start(j), strcat(date, '_0'));
            end
        end
    end 
    u_partStr = unique(partStr, 'stable');
    for ii = 1:length(u_partStr)
        compStr = strcmp(u_partStr{ii}, partStr);
        if sum(compStr) > 1
           ioStruct.Data{i}(ii,:) = mean(io{compStr,2:end},'omitnan');
           ioStruct.Sweep{i}(ii,1:sum(compStr)) = cellfun(@string, io.start(compStr))';
        else
          amp_find = io(compStr,:);
          ioStruct.Data{i}(ii,:) = table2array(amp_find(:,2:end));
          a = io.start(compStr); ioStruct.Sweep{i}(ii,1) = string(a(1));
        end
        ioStruct.stimIntensity{i}(ii,1:sum(compStr)) = stimIntensity{i}(compStr);
    end
    
    %ioStruct.Name{i} = io.Properties.VariableNames(:,2:end);
    ioStruct.Animal = mouse; ioStruct.Date = date; ioStruct.Cell = id;
end

%Get Stim ID Information
ioStruct.stimId = extractAfter(extractBetween(io.Properties.VariableNames(:,2:end), '_', '___'), '_');

% if sum(contains(io.Properties.VariableNames(:,2:end), 'b1')) > 0
% else if sum(contains(io.Properties.VariableNames(:,2:end), 'bb')) > 0
% 
% else
%     fprintf('cannot find stimId')
% end

ioStruct.kinetics = files_name;
ioStruct.blueId = pathway{string(pathway)=='blue',2};
ioStruct.redId = pathway{string(pathway)=='red',2};
ioStruct.in = pathway{string(pathway)=='in',2};
ioStruct.clamp = pathway{string(pathway)=='clamp',2};

if contains(ioStruct.clamp, 'cc', 'IgnoreCase', true)
    temp_idx = ioStruct.Data{contains(files_name, 'sum')} <= 0;
    ioStruct.Data{contains(files_name, 'sum')}(temp_idx) = NaN;

    temp_idx = ioStruct.Data{contains(files_name, 'amp')} <= 0;
    ioStruct.Data{contains(files_name, 'amp')}(temp_idx) = NaN;
end

cd(path);
save([path, mouse, id, '_ioStruct'], 'ioStruct', 'mouse', 'id', 'path', 'date');



%% Average across the same stim 
% average all sweeps with the same stim value and store in new part of
% structure (dataAvg), separate out ttx from non ttx conditions

% Average based on stim value 

for i = 1:length(files_name)  
    
    if sum(contains(ioStruct.Sweep{i}(:,1), 'ttx', 'IgnoreCase', true) > 0)
       temp_vec = contains(ioStruct.Sweep{i}(:,1), 'ttx', 'IgnoreCase', true);
       temp_data = ioStruct.Data{i}(temp_vec, :);
       u_stim = unique(ioStruct.stimIntensity{i}(temp_vec,1));
       %ampStruct.dataAvgStimTtx{i} = NaN(length(u_stim),50); ampStruct.dataAvgSweepTtx{i} = strings(length(u_stim),50);
       
       for j = 1:length(u_stim)
           comp_vec = ioStruct.stimIntensity{i}(temp_vec,1) == u_stim(j);
           if sum(comp_vec) > 1
               ioStruct.dataAvgTtx{i}(j,:) = mean(temp_data(comp_vec,:),'omitnan');
           else
               ioStruct.dataAvgTtx{i}(j,:) = temp_data(comp_vec,:);
           end
           temp_stim = ioStruct.stimIntensity{i}(comp_vec,:);
           ioStruct.dataAvgStimTtx{i}(j,1:length(temp_stim(:))) = temp_stim(:)';
           ts = ioStruct.Sweep{i}(temp_vec,:); temp_sweep = ts(comp_vec,:);
           ioStruct.dataAvgSweepTtx{i}(j,1:length(temp_sweep(:))) = temp_sweep(:)';
       end
       for k = 1:2:length(ioStruct.stimIdAll{i})
           ioStruct.pprAvgTtx{i}(:,k) = (ioStruct.dataAvgTtx{i}(:,k+1)./ioStruct.dataAvgTtx{i}(:,k));
       end
       temp_vec = ~contains(ioStruct.Sweep{i}(:,1), 'ttx', 'IgnoreCase', true);
       temp_data = ioStruct.Data{i}(temp_vec, :);
       u_stim = unique(ioStruct.stimIntensity{i}(temp_vec,1));
       ioStruct.dataAvgStim{i} = NaN(length(u_stim),10); 
       ioStruct.dataAvgSweep{i} = strings(length(u_stim),10);

       for j = 1:length(u_stim)
           comp_vec = ioStruct.stimIntensity{i}(temp_vec,1) == u_stim(j);
           temp_stim = ioStruct.stimIntensity{i}(comp_vec,:);
           ioStruct.dataAvgStim{i}(j,1:length(temp_stim(:))) = temp_stim(:)';
           ts = ioStruct.Sweep{i}(temp_vec,:); 
           temp_sweep = ts(comp_vec,:);
           ioStruct.dataAvgSweep{i}(j,1:length(temp_sweep(:))) = temp_sweep(:)';
           
           if sum(comp_vec) > 1
               ioStruct.dataAvg{i}(j,:) = mean(temp_data(comp_vec,:),'omitnan');
           else
               ioStruct.dataAvg{i}(j,:) = temp_data(comp_vec,:);
           end
       end  
    elseif sum(contains(ioStruct.Sweep{i}(:,1), idx_drug, 'IgnoreCase', true) > 0) 

       temp_vec = find(contains(ioStruct.Sweep{i}(:,1), idx_drug, 'IgnoreCase', true) == true);
       temp_data = ioStruct.Data{i}(temp_vec, :);
       u_stim = unique(ioStruct.stimIntensity{i}(temp_vec,1));
       %ampStruct.dataAvgStimTtx{i} = NaN(length(u_stim),50); ampStruct.dataAvgSweepTtx{i} = strings(length(u_stim),50);
       
       for j = 1:length(u_stim)
           comp_vec = ioStruct.stimIntensity{i}(temp_vec,1) == u_stim(j);
           if sum(comp_vec) > 1
               ioStruct.dataAvgDrug{i}(j,:) = mean(temp_data(comp_vec,:),'omitnan');
           else
               ioStruct.dataAvgDrug{i}(j,:) = temp_data(comp_vec,:);
           end
           temp_stim = ioStruct.stimIntensity{i}(comp_vec,:);
           ioStruct.dataAvgStimDrug{i}(j,1:length(temp_stim(:))) = temp_stim(:)';
           ts = ioStruct.Sweep{i}(temp_vec,:); temp_sweep = ts(comp_vec,:);
           ioStruct.dataAvgSweepDrug{i}(j,1:length(temp_sweep(:))) = temp_sweep(:)';
       end
       %for k = 1:2:10 need to be 1:2:length of number of columns
           %ioStruct.pprAvgDrug{i}(:,k) = (ioStruct.dataAvgDrug{i}(:,k+1)./ioStruct.dataAvgDrug{i}(:,k))
       %end
       temp_vec = ~contains(ioStruct.Sweep{i}(:,1), idx_drug, 'IgnoreCase', true);
       temp_data = ioStruct.Data{i}(temp_vec, :);
       u_stim = unique(ioStruct.stimIntensity{i}(temp_vec,1));
       ioStruct.dataAvgStim{i} = NaN(length(u_stim),10); ioStruct.dataAvgSweep{i} = strings(length(u_stim),10);

       for j = 1:length(u_stim)
           comp_vec = ioStruct.stimIntensity{i}(temp_vec,1) == u_stim(j);
           temp_stim = ioStruct.stimIntensity{i}(comp_vec,:);
           ioStruct.dataAvgStim{i}(j,1:length(temp_stim(:))) = temp_stim(:)';
           ts = ioStruct.Sweep{i}(temp_vec,:); temp_sweep = ts(comp_vec,:);
           ioStruct.dataAvgSweep{i}(j,1:length(temp_sweep(:))) = temp_sweep(:)';
           
           if sum(comp_vec) > 1
               ioStruct.dataAvg{i}(j,:) = mean(temp_data(comp_vec,:),'omitnan');
           else
               ioStruct.dataAvg{i}(j,:) = temp_data(comp_vec,:);
           end
       end

    else
       temp_vec = ~contains(ioStruct.Sweep{i}(:,1), 'ttx', 'IgnoreCase', true);
       temp_data = ioStruct.Data{i}(temp_vec, :);
       u_stim = unique(ioStruct.stimIntensity{i}(temp_vec,1));
       %ampStruct.dataAvgStim{i} = NaN(length(u_stim),10); ampStruct.dataAvgSweep{i} = strings(length(u_stim),10);

       for j = 1:length(u_stim)
           comp_vec = ioStruct.stimIntensity{i}(:,1) == u_stim(j);
           temp_stim = ioStruct.stimIntensity{i}(comp_vec,:);
           ioStruct.dataAvgStim{i}(j,1:length(temp_stim(:))) = temp_stim(:)';
           temp_sweep = ioStruct.Sweep{i}(comp_vec,:);
           ioStruct.dataAvgSweep{i}(j,1:length(temp_sweep(:))) = temp_sweep(:)';
           
           if sum(comp_vec) > 1
               ioStruct.dataAvg{i}(j,:) = mean(ioStruct.Data{i}(comp_vec,:),'omitnan');
           else
               ioStruct.dataAvg{i}(j,:) = ioStruct.Data{i}(comp_vec,:);
           end
       end  
    end
end

%failsafe to make an extra cell for drug io otherwise will error out if
%there is only drug io10 

if isfield(ioStruct, 'dataAvgDrug') && length(ioStruct.dataAvgDrug) ~= length(files_name)
   ioStruct.dataAvgDrug{length(files_name)} = [];
end

%% Calculate Various Metrics for Amplitude and Summation Values

%Find index for amp and sum
temp_idx = [find(contains(files_name, 'amp')); find(contains(files_name, 'sum'))];

%% PPR
if isfield(ioStruct, 'dataAvg')
    for i = 1:length(temp_idx)
        j = temp_idx(i);
        %cd(path2); io = readtable(files_name(temp_idx(i))); temp_num = length(extractAfter(extractBetween(io.Properties.VariableNames(:,2:end), '_', '___'), '_'))
        %ioStruct.stimIdAll{temp_idx(i)}
        [u, ~, g] = unique(ioStruct.stimIdAll{j}, 'stable');
        ioStruct.pprStimId{j} = u(accumarray(g, 1) > 1);
        %ioStruct.pprStimId{i} = ioStruct.stimIdAll{temp_idx(i)}(temp_idx2);
        % for ii = 1:2:length(ioStruct.stimIdAll{j})
        %     ioStruct.ppr{j}(:,ii) = abs((ioStruct.Data{j}(:,ii+1)./ioStruct.Data{j}(:,ii)));
        %     ioStruct.pprAvg{j}(:,ii) = abs((ioStruct.dataAvg{j}(:,ii+1)./ioStruct.dataAvg{j}(:,ii)));
        % end
        for ii = 1:length(ioStruct.pprStimId{j})
            tempIdx2 = find(string(ioStruct.stimIdAll{j}) == ioStruct.pprStimId{j}{ii}, 1,'first');
            if contains(ioStruct.clamp, 'cc', 'IgnoreCase', true)
                tempM = ioStruct.Data{j};
                tempM(find(ioStruct.Data{j} < 0.5)) = NaN; 
            else
                tempM = ioStruct.Data{j};
            end

            %ioStruct.ppr{j}(:,ii) = abs((ioStruct.Data{j}(:,tempIdx2+1)./ioStruct.Data{j}(:,tempIdx2)));
            ioStruct.ppr{j}(:,ii) = abs((tempM(:,tempIdx2+1)./tempM(:,tempIdx2)));
            if j <= length(ioStruct.dataAvg) && ~isempty(ioStruct.dataAvg{j})
                %ioStruct.pprAvg{j}(:,ii) = abs((ioStruct.dataAvg{j}(:,tempIdx2+1)./ioStruct.dataAvg{j}(:,tempIdx2)));
                if contains(ioStruct.clamp, 'cc', 'IgnoreCase', true)
                    tempM = ioStruct.dataAvg{j};
                    tempM(find(ioStruct.dataAvg{j} < 0.5)) = NaN; 
                else
                    tempM = ioStruct.dataAvg{j};
                end
                ioStruct.pprAvg{j}(:,ii) = abs((tempM(:,tempIdx2+1)./tempM(:,tempIdx2)));
            end
        end
    end
end

%cd ../

% FOR EVALUATING CROSS TALK, don't use this
% for i = 1:length(temp_idx)
%     temp_num(1) = length(ioStruct.dataAvg{1}(:,1)); %temp_num(2) = length(ampStruct.dataAvg{2}(:,1));
%     if ei & temp_num(1) ~= temp_num(2)
%         temp_vec = zeros(max([length(ioStruct.dataAvg{1}(:,1)), length(ioStruct.dataAvg{2}(:,1))]),1); temp_vec(temp_vec==0) = NaN;
%         ioStruct.pprAvgXTalk(:,i) = temp_vec;
%         ioStruct.pprAvgXTalk(1:temp_num(i),i) = (ioStruct.dataAvg{i}(:,4)./ioStruct.dataAvg{i}(:,9));
%     else 
%         %if empty(length(ampStruct.dataAvg{i}(:,4)) == length(ampStruct.dataAvg{i}(:,9))
%         ioStruct.pprAvgXTalk(:,i) = (ioStruct.dataAvg{i}(:,4)./ioStruct.dataAvg{i}(:,9));
%         %else
%            % diff_length = abs(length(ampStruct.dataAvg{i}(:,4)) - length(ampStruct.dataAvg{i}(:,9)))
%             
%     end 
%     save([path, mouse, id, '_ioStruct'], 'ioStruct', 'mouse', 'id', 'path', 'date');
% end

%% EI Ratio
if sum(contains(files_name, 'io10'))>0 && sum(contains(files_name, 'io80'))>0
    temp_idx = [find(contains(files_name, 'amp_io10')); find(contains(files_name, 'amp_io80'))];
    if sum(size(ioStruct.dataAvg{temp_idx(2)}) ==  size(ioStruct.dataAvg{temp_idx(1)})) > 1
        ioStruct.ei = abs(ioStruct.dataAvg{temp_idx(2)})./abs(ioStruct.dataAvg{temp_idx(1)});
        save([path, mouse, id, '_ioStruct'], 'ioStruct', 'mouse', 'id', 'path', 'date');
%     else
%       
%       for i = 1:length(temp_idx)  
%        temp_vec = ;
%        temp_data = ioStruct.Data{i}(temp_vec, :);
%        u_stim = unique(ioStruct.dataAvgStim{i}(temp_vec,1));
%        %ampStruct.dataAvgStim{i} = NaN(length(u_stim),10); ampStruct.dataAvgSweep{i} = strings(length(u_stim),10);
% 
%        for j = 1:length(u_stim)
%            comp_vec = ioStruct.stimIntensity{i}(:,1) == u_stim(j);
%            temp_stim = ioStruct.stimIntensity{i}(comp_vec,:);
%            ioStruct.dataAvgStim{i}(j,1:length(temp_stim(:))) = temp_stim(:)';
%            temp_sweep = ioStruct.Sweep{i}(comp_vec,:);
%            ioStruct.dataAvgSweep{i}(j,1:length(temp_sweep(:))) = temp_sweep(:)';
%        end
%       end

    end


elseif sum(contains(files_name, 'io10'))>0 && sum(contains(files_name, 'io70'))>0
    temp_idx = [find(contains(files_name, 'amp_io10')); find(contains(files_name, 'amp_io70'))];
    ioStruct.ei = abs(ioStruct.dataAvg{temp_idx(2)})./abs(ioStruct.dataAvg{temp_idx(1)});
    save([path, mouse, id, '_ioStruct'], 'ioStruct', 'mouse', 'id', 'path', 'date');
end


%% Average LED Responses Together

if sum(contains(ioStruct.stimId, 'b1')) > 0
    
%     if sum(contains(ioStruct.stimId, 'b1b2'))<1 || sum(contains(ioStruct.stimId, 'b3'))<1 || sum(contains(ioStruct.stimId, 'b4'))<1 || sum(contains(ioStruct.stimId, 'b5b6'))<1
%       error('something wrong with stim labels check file names');
%       return
%     end

    idx = find(contains(ioStruct.kinetics, 'sum') == false);
    idx2 = find(contains(ioStruct.kinetics, 'sum') == true);

    for j = 1:length(idx)
        i = idx(j); 
        ioStruct.dataBlue{i} = ioStruct.dataAvg{i}(:,contains(ioStruct.stimIdAll{i}, 'b')); 

        if sum(contains(ioStruct.stimId, 'dg'))>0
            ioStruct.dataDg{i} = ioStruct.dataAvg{i}(:,contains(ioStruct.stimIdAll{i}, 'dg'));
        end
        if sum(contains(ioStruct.stimId, 'rc'))>0
            ioStruct.dataRc{i} = ioStruct.dataAvg{i}(:,contains(ioStruct.stimIdAll{i}, 'rc'));
        end

        ti = [find(contains(ioStruct.stimIdAll{i}, 'b1b2') == true, 1, 'first'), find(contains(ioStruct.stimIdAll{i}, 'b3') == true, 1), find(contains(ioStruct.stimIdAll{i}, 'b4') == true, 1), find(contains(ioStruct.stimIdAll{i}, 'b5b6') == true, 1, 'first')];
        ti2 = [find(contains(ioStruct.stimIdAll{i}, 'b1b2') == true, 1, 'last'), find(contains(ioStruct.stimIdAll{i}, 'b5b6') == true, 1, 'last')];
        ioStruct.dataBlueAvg{i}(:,1) = mean(ioStruct.dataAvg{i}(:,ti), 2, 'omitnan');
        ioStruct.dataBlueAvg{i}(:,2) = mean(ioStruct.dataAvg{i}(:,ti2), 2, 'omitnan');
        ioStruct.dataBlueStimId{i} = ioStruct.stimIdAll{i}(contains(ioStruct.stimId, 'b'));

        %ioStruct.dataBlueAvg{i}(:,1) = mean([ioStruct.dataBlue{i}(:,1) ioStruct.dataBlue{i}(:,3) ioStruct.dataBlue{i}(:,4) ioStruct.dataBlue{i}(:,5)], 2, 'omitnan');
        %ioStruct.dataBlueAvg{i}(:,2) = mean([ioStruct.dataBlue{i}(:,2) ioStruct.dataBlue{i}(:,6)], 2, 'omitnan');  
        %ioStruct.dataBlueStimId{i} = ioStruct.stimIdAll{i}(contains(ioStruct.stimId, 'b'))
        if isfield(ioStruct, 'dataAvgDrug') && ~isempty(ioStruct.dataAvgDrug{i})
            ioStruct.dataBlueDrug{i} = ioStruct.dataAvgDrug{i}(:,contains(ioStruct.stimIdAll{i}, 'b'));
            ioStruct.dataBlueAvgDrug{i} = mean(ioStruct.dataAvgDrug{i}(:,ti), 2, 'omitnan');
            if sum(contains(ioStruct.stimId, 'dg'))>0
                ioStruct.dataDgDrug{i} = ioStruct.dataAvgDrug{i}(:,contains(ioStruct.stimIdAll{i}, 'dg'));
            end
            if sum(contains(ioStruct.stimId, 'rc'))>0
                ioStruct.dataRcDrug{i} = ioStruct.dataAvgDrug{i}(:,contains(ioStruct.stimIdAll{i}, 'rc'));
            end
        end

    end
   
    for j = 1:length(idx2)
        i = idx2(j);
        
        ioStruct.dataBlue{i} = ioStruct.dataAvg{i}(:,contains(ioStruct.stimIdAll{i}, 'b'));

        if sum(contains(ioStruct.stimId, 'dg'))>0
            ioStruct.dataDg{i} = ioStruct.dataAvg{i}(:,find(contains(ioStruct.stimIdAll{i}, 'dg') == true, 1, 'last'));
        end
        if sum(contains(ioStruct.stimId, 'rc'))>0
            ioStruct.dataRc{i} = ioStruct.dataAvg{i}(:,find(contains(ioStruct.stimIdAll{i}, 'rc') == true, 1, 'last'));
        end

        ti = [find(contains(ioStruct.stimIdAll{i}, 'b1b2') == true, 1, 'last'), find(contains(ioStruct.stimIdAll{i}, 'b5b6') == true, 1, 'last')];
        ioStruct.dataBlueAvg{i}(:,1) = mean(ioStruct.dataAvg{i}(:,ti), 2, 'omitnan');
        ioStruct.dataBlueStimId{i} = ioStruct.stimIdAll{i}(contains(ioStruct.stimId, 'b'));
        
        if isfield(ioStruct, 'dataAvgDrug') && ~isempty(ioStruct.dataAvgDrug{i})
            ioStruct.dataBlueDrug{i} = ioStruct.dataAvgDrug{i}(:,contains(ioStruct.stimIdAll{i}, 'b'));
            ioStruct.dataBlueAvgDrug{i} = mean(ioStruct.dataAvgDrug{i}(:,ti), 2, 'omitnan');
            if sum(contains(ioStruct.stimId, 'dg'))>0
                ioStruct.dataDgDrug{i} = ioStruct.dataAvgDrug{i}(:,find(contains(ioStruct.stimIdAll{i}, 'dg') == true, 1, 'last'));
            end
            if sum(contains(ioStruct.stimId, 'rc'))>0
                ioStruct.dataRcDrug{i} = ioStruct.dataAvgDrug{i}(:,find(contains(ioStruct.stimIdAll{i}, 'rc') == true, 1, 'last'));
            end
        end

    end


   %ioStruct.dataBlueStimId{i} = ioStruct.stimId(contains(ioStruct.stimId, 'b'));
   ioStruct.mW.blue = [2 4.4 6.6 8.7 10.7 12.7 14.6 16.4 18.3 20 28.3 35.2 40.6 45.1 48.8 52.1 55 57.6 60];



elseif sum(contains(ioStruct.stimId, 'bb')) > 0 
    
    idx = find(contains(ioStruct.kinetics, 'sum') == false); %idx for all kinetics but summation 
    idx2 = find(contains(ioStruct.kinetics, 'sum') == true); %idx for summation

    for j = 1:length(idx) %for all kinetics except summation
        i = idx(j); 
        if isfield(ioStruct, 'dataAvg')
            
            ioStruct.dataBlue{i} = [ioStruct.dataAvg{i}(:, contains(ioStruct.stimIdAll{i}, 'bb')) ioStruct.dataAvg{i}(:, find(string(ioStruct.stimIdAll{i}) == 'br', 1, 'first')) ioStruct.dataAvg{i}(:, find(string(ioStruct.stimIdAll{i}) == 'rb', 1, 'last'))];
            ioStruct.dataRed{i} = [ioStruct.dataAvg{i}(:, contains(ioStruct.stimIdAll{i}, 'rr')) ioStruct.dataAvg{i}(:, find(string(ioStruct.stimIdAll{i}) == 'br', 1, 'last')) ioStruct.dataAvg{i}(:, find(string(ioStruct.stimIdAll{i}) == 'rb', 1, 'first'))];
            ioStruct.dataCombo{i} = ioStruct.dataAvg{i}(:,contains(ioStruct.stimIdAll{i}, 'brbr'));
      
            if sum(contains(ioStruct.stimId, 'dg'))>0
                ioStruct.dataDg{i} = ioStruct.dataAvg{i}(:, contains(ioStruct.stimIdAll{i}, 'dg'));
            end
            
            if sum(contains(ioStruct.stimId, 'rc'))>0
                ioStruct.dataRc{i} = ioStruct.dataAvg{i}(:, contains(ioStruct.stimIdAll{i}, 'rc'));
            end
            
            ti = [find(contains(ioStruct.stimIdAll{i}, 'bb') == true, 1,'first'), find(string(ioStruct.stimIdAll{i}) == 'br', 1, 'first')];
            ti2 = [find(contains(ioStruct.stimIdAll{i}, 'rr') == true, 1,'first'), find(string(ioStruct.stimIdAll{i}) == 'rb', 1, 'first')];
            ioStruct.dataBlueAvg{i} = mean(ioStruct.dataAvg{i}(:,ti),2,'omitnan');
            ioStruct.dataRedAvg{i} = mean(ioStruct.dataAvg{i}(:,ti2),2,'omitnan');

        end
   
        if isfield(ioStruct, 'dataAvgDrug') && ~isempty(ioStruct.dataAvgDrug{i})
            ioStruct.dataBlueDrug{i} = [ioStruct.dataAvgDrug{i}(:, contains(ioStruct.stimIdAll{i}, 'bb')) ioStruct.dataAvgDrug{i}(:, find(string(ioStruct.stimIdAll{i}) == 'br', 1, 'first')) ioStruct.dataAvgDrug{i}(:, find(string(ioStruct.stimIdAll{i}) == 'rb', 1, 'last'))];
            ioStruct.dataRedDrug{i} = [ioStruct.dataAvgDrug{i}(:, contains(ioStruct.stimIdAll{i}, 'rr')) ioStruct.dataAvgDrug{i}(:, find(string(ioStruct.stimIdAll{i}) == 'br', 1, 'last')) ioStruct.dataAvgDrug{i}(:, find(string(ioStruct.stimIdAll{i}) == 'rb', 1, 'first'))];
            ioStruct.dataComboDrug{i} = ioStruct.dataAvgDrug{i}(:,contains(ioStruct.stimIdAll{i}, 'brbr'));
           
            
            ioStruct.dataBlueAvgDrug{i} = mean([ioStruct.dataBlueDrug{i}(:,1) ioStruct.dataBlueDrug{i}(:,3)],2,'omitnan');
            ioStruct.dataRedAvgDrug{i} = mean([ioStruct.dataRedDrug{i}(:,1) ioStruct.dataRedDrug{i}(:,4)],2,'omitnan');
            
            if sum(contains(ioStruct.stimId, 'dg'))>0
                ioStruct.dataDgDrug{i} = ioStruct.dataAvgDrug{i}(:, contains(ioStruct.stimIdAll{i}, 'dg'));
            end
            
            if sum(contains(ioStruct.stimId, 'rc'))>0
                ioStruct.dataRcDrug{i} = ioStruct.dataAvgDrug{i}(:, contains(ioStruct.stimIdAll{i}, 'rc'));
            end
        end
        if isfield(ioStruct, 'dataAvgTtx') && ~isempty(ioStruct.dataAvgTtx{i})
            ioStruct.dataBlueTtx{i} = [ioStruct.dataAvgTtx{i}(:, contains(ioStruct.stimIdAll{i}, 'bb')) ioStruct.dataAvgTtx{i}(:, find(string(ioStruct.stimIdAll{i}) == 'br', 1, 'first')) ioStruct.dataAvgTtx{i}(:, find(string(ioStruct.stimIdAll{i}) == 'rb', 1, 'last'))];
            ioStruct.dataRedTtx{i} = [ioStruct.dataAvgTtx{i}(:, contains(ioStruct.stimIdAll{i}, 'rr')) ioStruct.dataAvgTtx{i}(:, find(string(ioStruct.stimIdAll{i}) == 'br', 1, 'last')) ioStruct.dataAvgTtx{i}(:, find(string(ioStruct.stimIdAll{i}) == 'rb', 1, 'first'))];
            ioStruct.dataComboTtx{i} = ioStruct.dataAvgTtx{i}(:,contains(ioStruct.stimIdAll{i}, 'brbr'));
           
            
            ioStruct.dataBlueAvgTtx{i} = mean([ioStruct.dataBlueTtx{i}(:,1) ioStruct.dataBlueTtx{i}(:,3)],2,'omitnan');
            ioStruct.dataRedAvgTtx{i} = mean([ioStruct.dataRedTtx{i}(:,1) ioStruct.dataRedTtx{i}(:,4)],2,'omitnan');
            
            if sum(contains(ioStruct.stimId, 'dg'))>0
                ioStruct.dataDgTtx{i} = ioStruct.dataAvgTtx{i}(:, contains(ioStruct.stimIdAll{i}, 'dg'));
            end
            
            if sum(contains(ioStruct.stimId, 'rc'))>0
                ioStruct.dataRcTtx{i} = ioStruct.dataAvgTtx{i}(:, contains(ioStruct.stimIdAll{i}, 'rc'));
            end
        end
       
        ioStruct.dataBlueStimId{i} = [ioStruct.stimIdAll{i}(contains(ioStruct.stimIdAll{i}, 'bb')) ioStruct.stimIdAll{i}(find(string(ioStruct.stimIdAll{i}) == 'br', 1, 'first'))  ioStruct.stimIdAll{i}(find(string(ioStruct.stimIdAll{i}) == 'rb', 1, 'last'))];
        ioStruct.dataRedStimId{i} = [ioStruct.stimIdAll{i}(contains(ioStruct.stimIdAll{i}, 'rr')) ioStruct.stimIdAll{i}(find(string(ioStruct.stimIdAll{i}) == 'br', 1, 'last'))  ioStruct.stimIdAll{i}(find(string(ioStruct.stimIdAll{i}) == 'rb', 1, 'first'))];

    end
    for j = 1:length(idx2) %for summation data because you don't want the same averages
        i = idx2(j); 

        ioStruct.dataBlue{i} = [ioStruct.dataAvg{i}(:, contains(ioStruct.stimIdAll{i}, 'bb')) ioStruct.dataAvg{i}(:, find(string(ioStruct.stimIdAll{i}) == 'br', 1, 'first')) ioStruct.dataAvg{i}(:, find(string(ioStruct.stimIdAll{i}) == 'rb', 1, 'last'))];
        ioStruct.dataRed{i} = [ioStruct.dataAvg{i}(:, contains(ioStruct.stimIdAll{i}, 'rr')) ioStruct.dataAvg{i}(:, find(string(ioStruct.stimIdAll{i}) == 'br', 1, 'last')) ioStruct.dataAvg{i}(:, find(string(ioStruct.stimIdAll{i}) == 'rb', 1, 'first'))];
        ioStruct.dataCombo{i} = ioStruct.dataAvg{i}(:,contains(ioStruct.stimIdAll{i}, 'brbr'));

        if isfield(ioStruct, 'dataAvg')
            
            % ioStruct.dataBlue{i} = [ioStruct.dataAvg{i}(:, contains(ioStruct.stimIdAll{i}, 'bb')) ioStruct.dataAvg{i}(:, find(string(ioStruct.stimIdAll{i}) == 'br', 1, 'first')) ioStruct.dataAvg{i}(:, find(string(ioStruct.stimIdAll{i}) == 'rb', 1, 'last'))];
            % ioStruct.dataRed{i} = [ioStruct.dataAvg{i}(:, contains(ioStruct.stimIdAll{i}, 'rr')) ioStruct.dataAvg{i}(:, find(string(ioStruct.stimIdAll{i}) == 'br', 1, 'last')) ioStruct.dataAvg{i}(:, find(string(ioStruct.stimIdAll{i}) == 'rb', 1, 'first'))];
            % ioStruct.dataCombo{i} = ioStruct.dataAvg{i}(:,contains(ioStruct.stimIdAll{i}, 'brbr'));
            
            if sum(contains(ioStruct.stimId, 'dg'))>0
                ioStruct.dataDg{i} = ioStruct.dataAvg{i}(:,find(contains(ioStruct.stimIdAll{i}, 'dg') == true, 1, 'last'));
            end
            if sum(contains(ioStruct.stimId, 'rc'))>0
                ioStruct.dataRc{i} = ioStruct.dataAvg{i}(:,find(contains(ioStruct.stimIdAll{i}, 'rc') == true, 1, 'last'));
            end

            %ti = [find(contains(ioStruct.stimIdAll{i}, 'bb') == true, 1,'last')];
            %ti2 = [find(contains(ioStruct.stimIdAll{i}, 'rr') == true, 1,'last')];
            ioStruct.dataBlueAvg{i} = ioStruct.dataBlue{i}(:,2);
            ioStruct.dataRedAvg{i} = ioStruct.dataRed{i}(:,2);
            %ioStruct.dataBlueAvg{i} = mean(ioStruct.dataAvg{i}(:,ti),2,'omitnan');
            %ioStruct.dataRedAvg{i} = mean([ioStruct.dataAvg{i}(:,1) ioStruct.dataAvg{i}(:,4)],2,'omitnan');

        end
   
        if isfield(ioStruct, 'dataAvgDrug') && ~isempty(ioStruct.dataAvgDrug{i})
            ioStruct.dataBlueDrug{i} = [ioStruct.dataAvgDrug{i}(:, contains(ioStruct.stimIdAll{i}, 'bb')) ioStruct.dataAvgDrug{i}(:, find(string(ioStruct.stimIdAll{i}) == 'br', 1, 'first')) ioStruct.dataAvgDrug{i}(:, find(string(ioStruct.stimIdAll{i}) == 'rb', 1, 'last'))];
            ioStruct.dataRedDrug{i} = [ioStruct.dataAvgDrug{i}(:, contains(ioStruct.stimIdAll{i}, 'rr')) ioStruct.dataAvgDrug{i}(:, find(string(ioStruct.stimIdAll{i}) == 'br', 1, 'last')) ioStruct.dataAvgDrug{i}(:, find(string(ioStruct.stimIdAll{i}) == 'rb', 1, 'first'))];
            ioStruct.dataComboDrug{i} = ioStruct.dataAvgDrug{i}(:,contains(ioStruct.stimIdAll{i}, 'brbr'));
           
            ti = find(contains(ioStruct.stimIdAll{i}, 'bb') == true, 1, 'last');
            ti2 = find(contains(ioStruct.stimIdAll{i}, 'rr') == true, 1, 'last');
            ioStruct.dataBlueAvgDrug{i} = ioStruct.dataAvgDrug{i}(:, ti);
            ioStruct.dataRedAvgDrug{i} = ioStruct.dataAvgDrug{i}(:, ti2);

            if sum(contains(ioStruct.stimId, 'dg'))>0
                ioStruct.dataDgDrug{i} = ioStruct.dataAvgDrug{i}(:,find(contains(ioStruct.stimIdAll{i}, 'dg') == true, 1, 'last'));
            end
            if sum(contains(ioStruct.stimId, 'rc'))>0
                ioStruct.dataRcDrug{i} = ioStruct.dataAvgDrug{i}(:,find(contains(ioStruct.stimIdAll{i}, 'rc') == true, 1, 'last'));
            end
        end
        
        if isfield(ioStruct, 'dataAvgTtx') && ~isempty(ioStruct.dataAvgTtx{i})
            ioStruct.dataBlueTtx{i} = [ioStruct.dataAvgTtx{i}(:, contains(ioStruct.stimIdAll{i}, 'bb')) ioStruct.dataAvgTtx{i}(:, find(string(ioStruct.stimIdAll{i}) == 'br', 1, 'first')) ioStruct.dataAvgTtx{i}(:, find(string(ioStruct.stimIdAll{i}) == 'rb', 1, 'last'))];
            ioStruct.dataRedTtx{i} = [ioStruct.dataAvgTtx{i}(:, contains(ioStruct.stimIdAll{i}, 'rr')) ioStruct.dataAvgTtx{i}(:, find(string(ioStruct.stimIdAll{i}) == 'br', 1, 'last')) ioStruct.dataAvgTtx{i}(:, find(string(ioStruct.stimIdAll{i}) == 'rb', 1, 'first'))];
            ioStruct.dataComboTtx{i} = ioStruct.dataAvgTtx{i}(:,contains(ioStruct.stimIdAll{i}, 'brbr'));
           
            ti = find(contains(ioStruct.stimIdAll{i}, 'bb') == true, 1, 'last');
            ti2 = find(contains(ioStruct.stimIdAll{i}, 'rr') == true, 1, 'last');
            ioStruct.dataBlueAvgTtx{i} = ioStruct.dataAvgTtx{i}(:, ti);
            ioStruct.dataRedAvgTtx{i} = ioStruct.dataAvgTtx{i}(:, ti2);

            if sum(contains(ioStruct.stimId, 'dg'))>0
                ioStruct.dataDgTtx{i} = ioStruct.dataAvgTtx{i}(:,find(contains(ioStruct.stimIdAll{i}, 'dg') == true, 1, 'last'));
            end
            if sum(contains(ioStruct.stimId, 'rc'))>0
                ioStruct.dataRcAvgTtx{i} = ioStruct.dataAvgTtx{i}(:,find(contains(ioStruct.stimIdAll{i}, 'rc') == true, 1, 'last'));
            end
        end       
        ioStruct.dataBlueStimId{i} = [ioStruct.stimIdAll{i}(contains(ioStruct.stimIdAll{i}, 'bb')) ioStruct.stimIdAll{i}(find(string(ioStruct.stimIdAll{i}) == 'br', 1, 'first'))  ioStruct.stimIdAll{i}(find(string(ioStruct.stimIdAll{i}) == 'rb', 1, 'last'))];
        ioStruct.dataRedStimId{i} = [ioStruct.stimIdAll{i}(contains(ioStruct.stimIdAll{i}, 'rr')) ioStruct.stimIdAll{i}(find(string(ioStruct.stimIdAll{i}) == 'br', 1, 'last'))  ioStruct.stimIdAll{i}(find(string(ioStruct.stimIdAll{i}) == 'rb', 1, 'first'))];

    end
    
    ioStruct.mWRed = [1 2.4 3.7 5.1 6.4 7.8 9.1 10.5 11.8 13.1 19.8 26.1 31.5 35.9 39.6 43 46 48.6 51];
    ioStruct.mWBlue = [2 4.4 6.6 8.7 10.7 12.7 14.6 16.4 18.3 20 28.3 35.2 40.6 45.1 48.8 52.1 55 57.6 60];

else 
    fprintf('cannot find proper stimId information')
    return
end

save([path, mouse, id, '_ioStruct'], 'ioStruct', 'mouse', 'id', 'path', 'date');

%% Plotting Histograms of Responses
if his
    if sum(contains(ioStruct.stimId, 'b1')) > 0
        for i = 1:length(files_name)
            figure(i);
            histogram(ioStruct.dataBlue{i}, 'FaceColor', 'b'); hold on;
            histogram(ioStruct.dataDg{i}, 'FaceColor', 'g');
            histogram(ioStruct.dataRc{i}, 'FaceColor', 'y');
            legend('Blue', 'DG', 'RC');
            ylabel('Number of Responses'); xlabel(extractBefore(files_name(i), '2')); title(['Responses to Red/Blue LED Light', mouse, id]);
            h1.Normalization = 'probability'; h2.Normalization = 'probability';h3.Normalization = 'probability';
            h1.BinWidth = 0.3; h2.BinWidth = 0.3; h3.BinWidth = 0.3;
            hold off; 
        end
        
    elseif sum(contains(ioStruct.stimId), 'bb') > 0
        for i = 1:length(files)
%             figure(1);subplot(1,2,i)
%             histogram(ioStruct.dataAll{i}.combo(:,:), 'FaceColor', [.5 0 .5]); hold on;
%             histogram(ioStruct.dataAll{i}.blue(:,:), 'FaceColor', 'b');
%             histogram(ioStruct.dataAll{i}.red(:,:), 'FaceColor', 'r'); 
%             legend('Combined', 'Blue', 'Red');
%             xlabel('Number of Responses'); ylabel('Amplitude (pA)'); title(['Responses to Red/Blue LED Light', mouse, id]);
%             h1.Normalization = 'probability'; h2.Normalization = 'probability';h3.Normalization = 'probability';
%             h1.BinWidth = 0.3; h2.BinWidth = 0.3; h3.BinWidth = 0.3;
        end
    end
end

%% Plot Averages
if plt 
    if sum(contains(ioStruct.stimId, 'b1')) > 0
        for i = 1:length(files_name)
            figure(length(findobj('type', 'figure'))+1)
            plot(ioStruct.dataAvgStim{i}(1:length(ioStruct.dataBlue{i}),1)*10, ioStruct.dataBlueAvg{i}(:,1), 'b-o', 'MarkerSize', 8, 'LineWidth',2); hold on; 
            ti = [find(contains(ioStruct.stimIdAll{i}, 'b1b2') == true, 1, 'first'), find(contains(ioStruct.stimIdAll{i}, 'b3') == true, 1), find(contains(ioStruct.stimIdAll{i}, 'b4') == true, 1), find(contains(ioStruct.stimIdAll{i}, 'b5b6') == true, 1, 'first')];
            plot(ioStruct.dataAvgStim{i}(1:length(ioStruct.dataBlue{i}),1)*10, ioStruct.dataBlue{i}(:, ti), 'b.', 'MarkerSize', 10);
            if isfield(ioStruct, 'dataDg') && ~isempty(ioStruct.dataDg{i})
                plot(ioStruct.dataAvgStim{i}(1:length(ioStruct.dataDg{i}),1)*10, ioStruct.dataDg{i}(:,1), 'g-o', 'MarkerSize', 8, 'LineWidth',2);
            end
            if isfield(ioStruct, 'dataRc') && ~isempty(ioStruct.dataRc{i})
                plot(ioStruct.dataAvgStim{i}(1:length(ioStruct.dataRc{i}),1)*10, ioStruct.dataRc{i}(:,1), 'y-o', 'MarkerSize', 8, 'LineWidth',2);
            end
            
            plot(ioStruct.dataAvgStim{i}(1:length(ioStruct.dataBlue{i}),1)*10, ioStruct.dataBlue{i}(:, ti), 'b.', 'MarkerSize', 10);
            title(string(files_name(i))); legend(string(pathway{string(pathway)=='blue',2})); 
            xlabel('Stimulation Intensity (%)'); ylabel(string(extractBefore(files_name(i), '_i')));
           
        end
    elseif sum(contains(ioStruct.stimId, 'bb')) > 0
        for i = 1:length(files_name) 
            figure(length(findobj('type', 'figure'))+1)
            plot(ioStruct.dataAvgStim{i}(1:length(ioStruct.dataRed{i}),1)*10, ioStruct.dataRedAvg{i}, 'r-o', 'MarkerSize', 8, 'LineWidth', 2); hold on;
            plot(ioStruct.dataAvgStim{i}(1:length(ioStruct.dataRed{i}),1)*10, ioStruct.dataRed{i}(:,[1,4]), 'r.', 'MarkerSize', 10, 'LineWidth', 2); 
            plot(ioStruct.dataAvgStim{i}(1:length(ioStruct.dataBlue{i}),1)*10, ioStruct.dataBlueAvg{i}, 'b-o', 'MarkerSize', 8, 'LineWidth', 2);
            plot(ioStruct.dataAvgStim{i}(1:length(ioStruct.dataBlue{i}),1)*10, ioStruct.dataBlue{i}(:,[1,3]), 'b.', 'MarkerSize', 10, 'LineWidth', 2);
            plot(ioStruct.dataAvgStim{i}(1:length(ioStruct.dataCombo{i}),1)*10, ioStruct.dataCombo{i}(:,1), 'm-o', 'MarkerSize', 8, 'LineWidth', 2);
            if isfield(ioStruct, 'dataDg') && ~isempty(ioStruct.dataDg{i})
                plot(ioStruct.dataAvgStim{i}(1:length(ioStruct.dataDg{i}),1)*10, ioStruct.dataDg{i}(:,1), 'g-o', 'MarkerSize', 8, 'LineWidth', 2);
            end
            if isfield(ioStruct, 'dataRc') && ~isempty(ioStruct.dataRc{i})
                plot(ioStruct.dataAvgStim{i}(1:length(ioStruct.dataRc{i}),1)*10, ioStruct.dataRc{i}(:,1), 'y-o', 'MarkerSize', 8, 'LineWidth', 2);
            end
            xlabel('Stimulation Intensity (%)'); ylabel(string(extractBefore(files_name(i), '_i'))); title(['IO Curve @', ' ', mouse, ' ', id]);
        end
    end
end
%% Save 
save([path, mouse, id, '_ioStruct'], 'ioStruct', 'mouse', 'id', 'path', 'date');
clamp = ioStruct.clamp;

if isfield(ioStruct, 'ei')
    path2 = ['/Volumes/KO_Portable/MEC_LEC_CA3/analysis/',clamp,'/EI/']; cd(path2);
    eiAvg = ioStruct.ei;
    save([date, mouse, id, '_ei'], 'eiAvg');
end 
if isfield(ioStruct, 'pprAvg')
    path3 = ['/Volumes/KO_Portable/MEC_LEC_CA3/analysis/',clamp,'/PPR/']; cd(path3);
    pprAvg = ioStruct.pprAvg;
    save([date, mouse, id, '_ppr'], 'pprAvg');
end
path4 = ['/Volumes/KO_Portable/MEC_LEC_CA3/analysis/',clamp, '/IO/']; cd(path4);
save([date, mouse, id, '_ioStruct'], 'ioStruct', 'mouse', 'id', 'path', 'date');


