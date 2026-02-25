%% Working Draft of Open ABF Files for Ephys

close all; clear all;

mouse = 'KO209';
date = '20240625';
id = 'S3C3';;

%path = ['~/src/basu/data_converted/', date,'/',date,'_',mouse,'_',id,'/']; cd(path);
path = ['/Volumes/KO_Portable/MEC_LEC_CA3/data_converted/', date, '_', mouse,'/',date,'_',mouse,'_',id,'/']; cd(path);
num_files = length(dir([path, '2*.abf']));
files = dir([path, '2*abf']);

%% What would you like to do?
all_trace = 0; %do you want all traces w/ all channels plotted
organize = 0; %do you want to assign labels to your data?

ind_trace = 0; %do you want to look at individual traces 0 = no
file = 0; %which file do you want to look at individual traces

%% Generate Data Structure

if ~exist([path, mouse, id, '_dataStruct.mat'], 'file')
    dataStruct = [];
    for ii = 1:num_files
        file_name = [[path, files(ii).name]];
        [d, si, h] = abf2load(file_name);
        dataStruct{ii}.parameters = h; dataStruct{ii}.name = file_name; dataStruct{ii}.date = date;
        for i = 1:size(d,2)
            c = reshape(d(:,i,:),[size(d,1), size(d,3)]);
            dataStruct{ii}.chan{i} = c;
            if all_trace
                figure(ii); sgtitle([date files(ii).name]);
                subplot(1, size(d,2), i);
                plot(c);title(['channel ' num2str(i)])
            end
        end
        if organize 
            answer = inputdlg('Enter Protocol','Identify Traces', [1 50]);
            dataStruct{ii}.type = answer{1};
        end
        dataStruct{ii}.chanPar{1} = ['Recording']; dataStruct{ii}.chanPar{2} = ['Blue']; dataStruct{ii}.chanPar{3} = ['Red'];
    end
    save([path, mouse, id, '_dataStruct'], 'dataStruct', 'mouse', 'id', 'path', 'date');

else
    fprintf('data structure already exists');
end

%% Look at individual Traces

if ind_trace
    % Does the data structure already exist? 
    if ~exist([path, mouse, id, '_dataStruct.mat'], 'file')
        fprintf('data structure does not exist in current path')
    else
        load([path, mouse, id, '_dataStruct.mat']);
    end
    
    %Generate individual traces for all files
    if file == 0
        num_files = length(dir([path, '*.abf']));

        for j = 1:num_files;
            dataStruct{j}.chan{1};
            for jj = 1:data_struct{j}.parameters.lActualEpisodes;
                figure(j); sgtitle([dataStruct{j}.date 'file ' num2str(j)]);
                subplot(2, round(dataStruct{j}.parameters.lActualEpisodes/2), jj);
                plot(dataStruct{j}.chan{1}(:,jj));title(['Sweep ' num2str(jj)]);
            end
        end

    end
    
    %Generate individual traces for specified file
    if file > 0
       for j = 1:dataStruct{file}.parameters.lActualEpisodes;
       figure(j); sgtitle([dataStruct{file}.date 'file ' num2str(j)]);
       plot(dataStruct{file}.chan{1}(:,j));title(['Sweep ' num2str(j)])
       end
    end
end

%% Trash Code %%
%d = d(~startsWith({d.name}, '2P')); 
%d = d(~ismember({d.name},{'.','..',}));
%num_files = length(dir([['E:\data_converted\', date,'\'] '*.abf']))
%file_num = '2P_2021_12_02_0003'
%file_name = ['E:\data_converted\',date,'\',file_num,'.abf']
%[d, si, h] = abfload(file_name);
%path = ['E:\data_converted\', date,'\']
% d = dir(path);
% d = d(contains({d.name}, mouse))
% dir_num = length(d)

% for j = 1:length(d)
%     d(j).files = dir([path, d(j).name, '\', '*abf'])
% end
% = dir([['E:\data_converted\', date,'\'] '*.abf'])