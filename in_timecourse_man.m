%parent = dir(pwd); parent = parent([parent(:).isdir]);  
clear all; close all;

addpath '/Volumes/KO_Portable/MATLAB/EphysAnalysis_2025' ;
%do you want to make the remake the structure then set runOrg to 1
runOrg = 1;
tcPlot = 1;
loadTc = 1;
indvPlot = 1;
%where are the files you want to put into structure located
path = '/Volumes/KO_Portable/MEC_LEC_CA3/data_converted' ;
%path = 'R:\basulab\basulabspace\KO\MEC_LEC_CA3\data_converted';

%where are the files where the analysis structure is located
path2 = '/Volumes/KO_Portable/MEC_LEC_CA3/analysis/VC/IN' ;
%path2 = 'R:\basulab\basulabspace\KO\MEC_LEC_CA3\analysis\VC\IN';
%start in VC, IN directory and load the relevant structure
% 
% if loadTc 
%     cd(path2)
%     load('vip_tc.mat'); 
%     load('som_tc.mat'); 
%     load('pv_tc.mat'); 
%     load('cck_tc.mat');
% end
% 
% if exist('vip_tc.mat')
%     load('vip_tc.mat')
% else
%     vip_tc.amp = []
%     vip_tc.id = []
%     vip_tc.cno = []
%     vip_tc.time = []
%     vip_tc.ptwy = []
% end
% 
% if exist('cck_tc.mat')
%     load('cck_tc.mat')
% else
%     cck_tc.amp = []
%     cck_tc.id = []
%     cck_tc.cno = []
%     cck_tc.time = []
%     cck_tc.ptwy = []
% end
% 
% if exist('som_tc.mat')
%     load('som_tc.mat')
% else
%     som_tc.amp = []
%     som_tc.id = []
%     som_tc.cno = []
%     som_tc.time = []
%     som_tc.ptwy = []
% end
% 
% if exist('pv_tc.mat')
%     load('pv_tc.mat')
% else
%     pv_tc.amp = []
%     pv_tc.id = []
%     pv_tc.cno = []
%     pv_tc.time = []
%     pv_tc.ptwy = []
% end
%% 


if runOrg 
    cd(path); 
    pv_tc = []; som_tc.amp = []; vip_tc.amp = []; cck_tc. amp = [];
    pv_tc.amp = []; som_tc.amp = []; vip_tc.amp = []; cck_tc.amp = [];
    temp = dir(path); parent = temp(contains({temp(:).name}, {'KO'}));
    for i = 1:length(parent)
        cd(parent(i).name); 
        child = dir(pwd); child = child([child(:).isdir]); 
        child = child(~ismember({child(:).name},{'.','..'}));
%         if isfile('pathway')
%             pathway = readcell('pathway'); 
%         else
%             fprintf('pathway txt file does not exist'); pwd
%             fileID = fopen('pathway', 'w'); ptwy_temp = inputdlg({'blue', 'red'});
%             ptwy_temp = [{['blue_',ptwy_temp{1}]}; {['red_', ptwy_temp{2}]}];
%             fprintf(fileID,  '%s\n', ptwy_temp{:});
%             pathway = readcell('pathway')
%         end
        for ii = 1:length(child) 
            cd(child(ii).name);
            if isfile('pathway.txt')
                pathway = readcell('pathway.txt'); 
            else
                fprintf('pathway txt file does not exist'); pwd
                fileID = fopen('pathway', 'w'); ptwy_temp = inputdlg({'blue', 'red'});
                ptwy_temp = [{['blue_',ptwy_temp{1}]}; {['red_', ptwy_temp{2}]}];
                fprintf(fileID,  '%s\n', ptwy_temp{:});
                pathway = readcell('pathway');
            end
            if sum(cell2mat(pathway(string(pathway)=='in',2)))~= 0 && isfolder('baseline') && cell2mat(pathway(string(pathway)=='include',2)) ~= 0
               %in = readcell('in'); in_type = cell2mat(in(contains(string(in), 'in'),2));
               in_type = cell2mat(pathway(matches(string(pathway), 'in'),2));
               cd baseline; files = dir(['amp_baseline*']);
               %tc = readcell(files.name, "FileType", "spreadsheet"); 
               tc = readcell(files.name,'Delimiter', ',');
               id = string(extractBetween(files.name, 'KO', '.'));
               cno = cell2mat(pathway(matches(string(pathway), 'cno'),2));
               time_vec = -(cno-1):1:(length(tc)-cno-1); 
               % if isfile('pathway')
               %      pathway = readcell('pathway'); 
               %  else
               %      fprintf('pathway txt file does not exist'); pwd
               %      fileID = fopen('pathway', 'w'); ptwy_temp = inputdlg({'blue', 'red'});
               %      ptwy_temp = [{['blue_',ptwy_temp{1}]}; {['red_', ptwy_temp{2}]}];
               %      fprintf(fileID,  '%s\n', ptwy_temp{:});
               %      pathway = readcell('pathway')
               %  end

               if contains(in_type, 'vip') && cno > 0
                    idx = length(vip_tc.amp) + 1;
                    vip_tc.amp{idx} = tc(:, 2:end);
                    vip_tc.in(idx) = string(in_type);
                    vip_tc.id{idx} = id;
                    vip_tc.cno(idx) = cno;
                    vip_tc.time{idx} = time_vec; %vip_tc.time{2,idx} = time_vec/0.25;
                    vip_tc.ptwy{idx,1} = string(pathway{string(pathway)=='blue',2}); vip_tc.ptwy{idx, 2} = string(pathway{string(pathway)=='red',2});
               elseif contains(in_type, 'cck') && cno > 0
                    idx = length(cck_tc.amp) + 1;
                    cck_tc.amp{idx} = tc(:, 2:end);
                    cck_tc.in(idx) = string(in_type);
                    cck_tc.id{idx} = id;
                    cck_tc.cno(idx) = cno;
                    cck_tc.time{idx} = time_vec; %cck_tc.time{2,idx} = time_vec/0.25;
                    cck_tc.ptwy{idx,1} = string(pathway{string(pathway)=='blue',2}); cck_tc.ptwy{idx, 2} = string(pathway{string(pathway)=='red',2});
               elseif contains(in_type, 'pv') && cno > 0
                    idx = length(pv_tc.amp) + 1;
                    pv_tc.amp{idx} = tc(:, 2:end);
                    pv_tc.in(idx) = string(in_type);
                    pv_tc.id{idx} = id;
                    pv_tc.cno(idx) = cno;
                    pv_tc.time{idx} = time_vec; %pv_tc.time{2,idx} = time_vec/0.25;
                    pv_tc.ptwy{idx,1} = string(pathway{string(pathway)=='blue',2}); pv_tc.ptwy{idx, 2} = string(pathway{string(pathway)=='red',2});
               elseif contains(in_type, 'som') && cno > 0
                    idx = length(som_tc.amp) + 1;
                    som_tc.amp{idx} = tc(:, 2:end);
                    som_tc.in(idx) = string(in_type);
                    som_tc.id{idx} = id;
                    som_tc.cno(idx) = cno;
                    som_tc.time{idx} = time_vec; %som_tc.time{2,idx} = time_vec/0.25;
                    som_tc.ptwy{idx,1} = string(pathway{string(pathway)=='blue',2}); som_tc.ptwy{idx, 2} = string(pathway{string(pathway)=='red',2});

               end
               cd ..
            end
            cd .. 
        end
        cd ..
    end
    cd(path2);
    save(['vip_tc'], 'vip_tc'); save(['pv_tc'], 'pv_tc'); save(['som_tc'], 'som_tc'); save(['cck_tc'], 'cck_tc');
end

%% Plotting

if tcPlot
tc_all = who('*_tc');

    for i = 1:length(tc_all)
        temp = [];
        temp = eval(tc_all{i});
        figure(i)
        sgtitle(string(tc_all(i)))
        in_tc{i} = []; in_tc{i}.id = temp.id; in_tc{i}.in = temp.in; in_tc{i}.cno = temp.cno; in_tc{i}.ptwy = temp.ptwy;
        m_temp = NaN(max((cellfun(@(x) length(x), temp.time) - temp.cno)+21), length(temp.time)); m2_temp = m_temp; m3_temp = m_temp;
        for ii = 1:length(temp.amp)
            idx = temp.cno(ii);
            in_tc{i}.time{ii} = temp.time{ii}((idx-20):end);
            in_tc{i}.ampKey{ii} = string(cellfun(@(x) extractAfter(x, 'line_'), temp.amp{ii}(1,:), 'UniformOutput',false));
            in_tc{i}.amp{ii} = temp.amp{ii}((idx-19):end, :);
            in_tc{i}.normVal{ii} = mean(cell2mat(temp.amp{ii}((idx-19):idx, :)));
            in_tc{i}.ampNorm{ii} = cell2mat(in_tc{i}.amp{ii})./in_tc{i}.normVal{ii};
            in_tc{i}.ampBox{ii} = [movmean(in_tc{i}.ampNorm{ii}(1:20,:),4); movmean(in_tc{i}.ampNorm{ii}(21:end,:),4)];
           
            m_temp(1:length(in_tc{i}.ampBox{ii}),ii) = in_tc{i}.ampBox{ii}(:,1);
            m2_temp(1:length(in_tc{i}.ampNorm{ii}),ii) = in_tc{i}.ampNorm{ii}(:,1);
            m3_temp(1:length(in_tc{i}.amp{ii}),ii) = cell2mat(in_tc{i}.amp{ii}(:,1));
            
            subplot(3,1,1); plot(temp.time{ii}, cell2mat(temp.amp{ii}(2:end,1))); hold on;
            title('Raw Amp');
            subplot(3,1,2); plot(in_tc{i}.time{ii}*0.25, in_tc{i}.ampNorm{ii}(:,1)); hold on;
            title('Norm Amp');
            subplot(3,1,3); plot(in_tc{i}.time{ii}*0.25, in_tc{i}.ampBox{ii}(:,1)); hold on;
            title('BoxCar Amp');
            %subplot(1,2,3); 
        end
        
        %% YOU ARE HERE
        in_tc{i}.mAmpBoxMat = m_temp;
        [u, ~,g] = unique(string(in_tc{i}.ptwy(:,1)),'stable');
        in_tc{i}.mAmpPtwy = u;
        in_tc{i}.mAmpBox = splitapply(@(x) mean(m_temp, 2, 'omitnan'), m_temp, g');
        in_tc{i}.stdAmpBox = splitapply(@(x) rowStdOrNaN(x, 2), m_temp, g');
        in_tc{i}.nAmpBox = splitapply(@(x) sum(~isnan(x),2), m_temp, g');
        %in_tc{i}.semAmpBox =  in_tc{i}.stdAmpBox./sqrt(cellfun(@(x) nnz(~isnan(x)), num2cell(m_temp,2)));
        in_tc{i}.semAmpBox = in_tc{i}.stdAmpBox./sqrt(in_tc{i}.nAmpBox);


        in_tc{i}.mAmpNormMat = m2_temp;
        in_tc{i}.mAmpNorm = splitapply(@(x) mean(x, 2, 'omitnan'), m2_temp, g');
        in_tc{i}.stdAmpNorm = splitapply(@(x) rowStdOrNaN(x, 2), m2_temp, g');
        in_tc{i}.nAmpNorm = splitapply(@(x) sum(~isnan(x),2), m2_temp, g');
        in_tc{i}.semAmpNorm = in_tc{i}.stdAmpNorm./sqrt(in_tc{i}.nAmpNorm);
        %in_tc{i}.semAmpNorm = std(m2_temp,0, 2, 'omitnan')./sqrt(cellfun(@(x) nnz(~isnan(x)), num2cell(m2_temp,2)));

        in_tc{i}.AmpIndv = m3_temp;
        in_tc{i}.mAmpIndv = [mean(m3_temp(10:20,:)); mean(m3_temp(55:65,:), 'omitnan')];

        if sum(isnan(in_tc{i}.mAmpIndv),'all') > 0
           [row col] = find(isnan(in_tc{i}.mAmpIndv));
           if row == 2
               tempIdx = sum(~isnan(m3_temp(:,col)));
               in_tc{i}.mAmpIndv(row, col) = mean(m3_temp(tempIdx-10:tempIdx, col));
           end
        end
        %in_tc{i}.mAmpIndv = [splitapply(@(x) mean(x, 'omitnan'), m3_temp(10:20,:), g');splitapply(@(x) mean(x, 2, 'omitnan'), m3_temp(60:70,:), g')]
        in_tc{i}.mAmpAll = splitapply(@(x) mean(x, 2, 'omitnan'), in_tc{i}.mAmpIndv, g');
        in_tc{i}.stdAmpAll = splitapply(@(x) std(x,0,2, 'omitnan'), in_tc{i}.mAmpIndv, g');
        in_tc{i}.nAmpAll = splitapply(@(x) sum(~isnan(x),2), in_tc{i}.mAmpIndv, g');
        in_tc{i}.semAmpAll = in_tc{i}.stdAmpAll ./ sqrt(in_tc{i}.nAmpAll);
        %in_tc{i}.mAmpAll = splitapply(@(x, y) [mean(reshape(x, 1, []),'omitnan');  mean(reshape(y, 1, []),'omitnan')], m3_temp(10:20,:), m3_temp(60:70,:), g');
        %in_tc{i}.stdAmpAll = splitapply(@(x, y) [std(reshape(x, 1, []),'omitnan');  std(reshape(y, 1, []),'omitnan')], m3_temp(10:20,:), m3_temp(60:70,:), g');
        %in_tc{i}.nAmpAll = splitapply(@(x, y) [nnz(~isnan(reshape(x, 1, []))); nnz(~isnan(reshape(y, 1, [])))], m3_temp(10:20,:), m3_temp(60:70,:), g');
        %in_tc{i}.semAmpAll = in_tc{i}.stdAmpAll./sqrt(in_tc{i}.nAmpAll);

        [a, b] = max(cellfun(@(x) length(x), in_tc{i}.time)); in_tc{i}.mTimeAmp = in_tc{i}.time{b};

        figure(i)
        tempIdx = find(u == "mec");
        subplot(3,1,2); errorbar(in_tc{i}.mTimeAmp*0.25, in_tc{i}.mAmpNorm(:,tempIdx), in_tc{i}.semAmpNorm(:,tempIdx), 'k');
        subplot(3,1,3); errorbar(in_tc{i}.mTimeAmp*0.25, in_tc{i}.mAmpBox(:,tempIdx), in_tc{i}.semAmpBox(:,tempIdx), 'k');
    end

end

% for i = 1:length(in_tc)
%     in_tc{i}.mAmpBoxIndv = [mean(in_tc{i}.mAmpBoxMat(10:20,:)); mean(in_tc{i}.mAmpBoxMat(40:50,:))];
%     in_tc{i}.mAmpBoxAll = mean(in_tc{i}.mAmpBoxIndv, 2);
%     in_tc{i}.semAmpBoxAll = std(in_tc{i}.mAmpBoxIndv, 0 , 2)./sqrt(length(in_tc{i}.mAmpBoxIndv));
% end

% for i = 1:length(in_tc)
%     in_tc{i}.mAmpIndv = [mean(in_tc{i}.amp(10:20,:)); mean(in_tc{i}.amp(40:50,:))];
%     in_tc{i}.mAmpAll = mean(in_tc{i}.mAmpIndv, 2);
%     in_tc{i}.semAmpAll = std(in_tc{i}.mAmpIndv, 0 , 2)./sqrt(length(in_tc{i}.mAmpIndv));
% end



%% VR average every four, need to be manually changed if add cells

% in_tc{1}.avgAmpNorm = reshape(mean(reshape(in_tc{1}.mAmpNormMat(1:end-3,:), 4,[]),'omitnan'), [], length(in_tc{1}.id));
% in_tc{2}.avgAmpNorm = reshape(mean(reshape(in_tc{2}.mAmpNormMat(:,:), 4,[]),'omitnan'), [], length(in_tc{2}.id));
% in_tc{3}.avgAmpNorm = reshape(mean(reshape(in_tc{3}.mAmpNormMat(1:end-1,:), 4,[]),'omitnan'), [], length(in_tc{3}.id));
% in_tc{4}.avgAmpNorm = reshape(mean(reshape(in_tc{4}.mAmpNormMat(1:152,:), 4,[]),'omitnan'), [], length(in_tc{4}.id));

for i = 1:length(in_tc)
    [row,col] = size(in_tc{i}.mAmpNormMat);
    %Find closest row number divisible by 4 to reshape
    if abs(row - mod(row,4)) <= abs(row + (4 - mod(row,4)))
        idx = row - mod(row,4);
    else
        idx = row + (4 - mod(row,4));
    end
    for ii = 1:col
        in_tc{i}.avgAmpNorm(:,ii) = mean(reshape(in_tc{i}.mAmpNormMat(1:idx,ii), 4,[]), 'omitnan')' ;
    end
end

for i = 1:length(in_tc)
    [u, ~, g] = unique(string(in_tc{i}.ptwy(:,1)), 'stable');
    in_tc{i}.mAvgAmpNorm = splitapply(@(x) mean(x, 2, 'omitnan'), in_tc{i}.avgAmpNorm, g');
    in_tc{i}.sdAvgAmpNorm = splitapply(@(x) rowStdOrNaN(x, 2), in_tc{i}.avgAmpNorm, g');
    in_tc{i}.nAvgAmpNorm =  splitapply(@(x) sum(~isnan(x),2), in_tc{i}.avgAmpNorm, g');
    in_tc{i}.semAvgAmpNorm = in_tc{i}.sdAvgAmpNorm./sqrt(in_tc{i}.nAvgAmpNorm);
    in_tc{i}.AvgAmpNormPtwy = u;
    in_tc{i}.timeAvgAmpNorm = in_tc{i}.mTimeAmp(1:4:end)*0.25;
end
% for i = 1:length(in_tc)
%     in_tc{i}.mAvgAmpNorm = mean(in_tc{i}.avgAmpNorm,2,'omitnan');
%     in_tc{i}.semAvgAmpNorm = std(in_tc{i}.avgAmpNorm,0, 2, 'omitnan')./sqrt(cellfun(@(x) nnz(~isnan(x)), num2cell(in_tc{i}.avgAmpNorm,2)));
%     in_tc{i}.timeAvgAmpNorm = in_tc{i}.mTimeAmp(1:4:end)*0.25;
% end

save(['in_tc'], 'in_tc');



%% Plot Individual TC 

if indvPlot

    figure
    for i = 1:length(in_tc)
        idx = find(string(in_tc{i}.ptwy(:,1)) == 'mec');
        idx2 = find(in_tc{i}.AvgAmpNormPtwy == 'mec');
        subplot(2,2,i)
        plot(in_tc{i}.mTimeAmp*0.25 , in_tc{i}.mAmpBoxMat(:,idx), '.', 'MarkerSize',10); hold on;
        errorbar(in_tc{i}.mTimeAmp*0.25, in_tc{i}.mAmpBox(:,idx2), in_tc{i}.semAmpBox(:,idx2), 'k.-', 'MarkerSize',30)
        xlim([-5 15]); ylim([0.5 1.5]); title(in_tc{i}.in{1}); sgtitle('Box Car Average');
    end
    figure
    for i = 1:length(in_tc)
        idx = find(string(in_tc{i}.ptwy(:,1)) == 'mec');
        idx2 = find(in_tc{i}.mAmpPtwy == 'mec');
        subplot(2,2,i)
        plot(in_tc{i}.mTimeAmp*0.25 , in_tc{i}.mAmpNormMat(:,idx), '.', 'MarkerSize',10); hold on;
        errorbar(in_tc{i}.mTimeAmp*0.25, in_tc{i}.mAmpNorm(:,idx2), in_tc{i}.semAmpNorm(:,idx2), 'k.-', 'MarkerSize',30)
        xlim([-5 15]); ylim([0.5 1.5]); title(in_tc{i}.in{1}); sgtitle('Normalized Average');
    end
    %Before and After CNO Plots
    figure
    for i = 1:length(in_tc)
        idx = find(string(in_tc{i}.ptwy(:,1)) == 'mec');
        idx2 = find(in_tc{i}.mAmpPtwy == 'mec');
        figure
        plot(1:2, in_tc{i}.mAmpIndv(:,idx)', '.-', 'MarkerSize',10,'Color', [0.7 0.7 0.7],'LineWidth',4); hold on;
        errorbar(1:2, in_tc{i}.mAmpAll(:,idx2), in_tc{i}.semAmpAll(:,idx2), 'r.', 'MarkerSize',50, 'LineWidth', 2)
        xlim([0.5 2.5]); ylim([0 2800]);title(in_tc{i}.in{1}); ylabel('IPSC Amplitude (pA)');
        set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
        xticklabels({'','ACSF','', 'CNO'}); xlabel('IPSC amplitude (pA)');
    end

figure
sgtitle('Averaged across min');
subplot(2,2,1)
idx = find(string(in_tc{1}.ptwy(:,1)) == 'mec');
idx2 = find(in_tc{1}.mAmpPtwy == 'mec');
plot(in_tc{1}.timeAvgAmpNorm(1:end-1), in_tc{1}.avgAmpNorm(:,idx));
hold on; errorbar(in_tc{1}.timeAvgAmpNorm(1:end-1), in_tc{1}.mAvgAmpNorm(:,idx2), in_tc{1}.semAvgAmpNorm(:,idx2), 'k.-', 'MarkerSize', 30);
xlim([-5, 15]); ylim([0.5 1.5]); title(in_tc{1}.in{1}); 


subplot(2,2,2)
idx = find(string(in_tc{2}.ptwy(:,1)) == 'mec');
idx2 = find(in_tc{2}.mAmpPtwy == 'mec');
plot(in_tc{2}.timeAvgAmpNorm(1:end-1), in_tc{2}.avgAmpNorm(:,idx));
hold on; errorbar(in_tc{2}.timeAvgAmpNorm(1:end-1), in_tc{2}.mAvgAmpNorm(:,idx2), in_tc{2}.semAvgAmpNorm(:,idx2), 'k.-', 'MarkerSize', 30);
xlim([-5, 15]); ylim([0.5 1.5]);title(in_tc{2}.in{1});

subplot(2,2,3)
idx = find(string(in_tc{3}.ptwy(:,1)) == 'mec');
idx2 = find(in_tc{3}.mAmpPtwy == 'mec');
plot(in_tc{3}.timeAvgAmpNorm(1:end-1), in_tc{3}.avgAmpNorm(:,idx));
hold on; errorbar(in_tc{3}.timeAvgAmpNorm(1:end-1), in_tc{3}.mAvgAmpNorm(:,idx2), in_tc{3}.semAvgAmpNorm(:,idx2), 'k.-', 'MarkerSize', 30);
xlim([-5, 15]); ylim([0.5 1.5]);title(in_tc{3}.in{1});

subplot(2,2,4)
idx = find(string(in_tc{4}.ptwy(:,1)) == 'mec');
idx2 = find(in_tc{4}.mAmpPtwy == 'mec');
plot(in_tc{4}.timeAvgAmpNorm(1:end-1), in_tc{4}.avgAmpNorm(:,idx));
hold on; errorbar(in_tc{4}.timeAvgAmpNorm(1:end-1), in_tc{4}.mAvgAmpNorm(:,idx2), in_tc{4}.semAvgAmpNorm(:,idx2), 'k.-', 'MarkerSize', 30);
xlim([-5, 15]); ylim([0.5 1.5]);title(in_tc{4}.in{1});

end
%% 

figure
idx = find(string(in_tc{2}.ptwy(:,1)) == 'mec');
idx2 = find(in_tc{2}.mAmpPtwy == 'mec');
yline(100, '--', 'Color', [.7 .7 .7],'LineWidth',2);
hold on; errorbar(in_tc{2}.timeAvgAmpNorm(1:length(in_tc{2}.mAvgAmpNorm)), (100*(in_tc{2}.mAvgAmpNorm(:,idx2))), 100*(in_tc{2}.semAvgAmpNorm(:,idx2)), 'b.', 'MarkerSize', 50,'LineWidth',3);
xlim([-5, 17]); ylim([0 120]); title(in_tc{2}.in{1}); xlabel('time (min)'); ylabel('IPSC amplitude (%)');
set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
%annotation('line', [0.305 .905], [0.86 0.86], 'LineWidth', 5, 'Color',[.7 .7 .7]);
%text(7.8, 115, 'CNO', 'FontSize',14, 'FontWeight', 'bold', 'Color', [0.7 0.7 0.7]);

figure
idx = find(string(in_tc{1}.ptwy(:,1)) == 'mec');
idx2 = find(in_tc{1}.mAmpPtwy == 'mec');
yline(100, '--', 'Color', [.7 .7 .7],'LineWidth',2);
hold on; errorbar(in_tc{1}.timeAvgAmpNorm(1:end-1), (100*(in_tc{1}.mAvgAmpNorm(:,idx2))), 100*in_tc{1}.semAvgAmpNorm(:,idx2), '.', 'MarkerSize', 50, 'Color', "#EDB120",'LineWidth',3);
xlim([-5, 17]); ylim([0 120]);title(in_tc{1}.in{1}); xlabel('time (min)'); ylabel('IPSC amplitude (%)');
set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
%annotation('line', [0.305 .905], [0.86 0.86], 'LineWidth', 5, 'Color',[.7 .7 .7]);
%text(7.8, 115, 'CNO', 'FontSize',14, 'FontWeight', 'bold', 'Color', [0.7 0.7 0.7]);

figure
idx = find(string(in_tc{3}.ptwy(:,1)) == 'mec');
idx2 = find(in_tc{3}.mAmpPtwy == 'mec');
yline(100, '--', 'Color', [.7 .7 .7],'LineWidth',2);
hold on; errorbar(in_tc{3}.timeAvgAmpNorm(1:end-1), (100*(in_tc{3}.mAvgAmpNorm(:,idx2))), 100*in_tc{3}.semAvgAmpNorm(:,idx2), '.', 'MarkerSize', 50, 'Color',"#D95319", 'LineWidth',3);
xlim([-5, 17]); ylim([0 120]);title(in_tc{3}.in{1}); xlabel('time (min)'); ylabel('IPSC amplitude (%)');
set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
%annotation('line', [0.305 .905], [0.86 0.86], 'LineWidth', 5, 'Color',[.7 .7 .7]);
%text(7.8, 115, 'CNO', 'FontSize',14, 'FontWeight', 'bold', 'Color', [0.7 0.7 0.7]);

figure
idx = find(string(in_tc{4}.ptwy(:,1)) == 'mec');
idx2 = find(in_tc{4}.mAmpPtwy == 'mec');
yline(100, '--', 'Color', [.7 .7 .7],'LineWidth',2);
hold on; errorbar(in_tc{4}.timeAvgAmpNorm(1:end-1), (100*(in_tc{4}.mAvgAmpNorm(:,idx2))), 100*in_tc{4}.semAvgAmpNorm(:,idx2), '.', 'MarkerSize', 50, 'Color', "#7E2F8E", 'LineWidth', 3);
xlim([-5, 17]); ylim([0 120]);title(in_tc{4}.in{1}); xlabel('time (min)'); ylabel('IPSC amplitude (%)');
set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
%annotation('line', [0.305 .905], [0.86 0.86], 'LineWidth', 5, 'Color',[.7 .7 .7]);
%text(7.8, 115, 'CNO', 'FontSize',14, 'FontWeight', 'bold', 'Color', [0.7 0.7 0.7]);


%% Plot VR Data for LEC 

load('lec_tc.mat');

figure
plot(1:2, lecCCK.plot,'.-', 'MarkerSize',10,'Color', [0.7 0.7 0.7],'LineWidth',4); hold on;
errorbar(1:2, [lecCCK.mACSF, lecCCK.mCNO], [lecCCK.semACSF, lecCCK.semCNO], '.', 'MarkerSize',50, 'LineWidth', 2, 'Color','#127B2E')
xlim([0.5 2.5]); ylim([0 4500]); ylabel('IPSC Amplitude (pA)'); title('cck');
set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
xticklabels({'','ACSF','', 'CNO'}); xlabel('IPSC amplitude (pA)');

figure
plot(1:2, lecPV.plot,'.-', 'MarkerSize',10,'Color', [0.7 0.7 0.7],'LineWidth',4); hold on;
errorbar(1:2, [lecPV.mACSF, lecPV.mCNO], [lecPV.semACSF, lecPV.semCNO], '.', 'MarkerSize',50, 'LineWidth', 2, 'Color','#127B2E')
xlim([0.5 2.5]); ylim([0 2500]); ylabel('IPSC Amplitude (pA)'); title('pv');
set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
xticklabels({'','ACSF','', 'CNO'}); xlabel('IPSC amplitude (pA)');

figure
plot(1:2, lecSST.plot,'.-', 'MarkerSize',10,'Color', [0.7 0.7 0.7],'LineWidth',4); hold on;
errorbar(1:2, [lecSST.mACSF, lecSST.mCNO], [lecSST.semACSF, lecSST.semCNO], '.', 'MarkerSize',50, 'LineWidth', 2, 'Color','#127B2E')
xlim([0.5 2.5]); ylim([0 3200]); ylabel('IPSC Amplitude (pA)'); title('sst');
set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
xticklabels({'','ACSF','', 'CNO'}); xlabel('IPSC amplitude (pA)');

figure
plot(1:2, lecVIP.plot,'.-', 'MarkerSize',10,'Color', [0.7 0.7 0.7],'LineWidth',4); hold on;
errorbar(1:2, [lecVIP.mACSF, lecVIP.mCNO], [lecVIP.semACSF, lecVIP.semCNO], '.', 'MarkerSize',50, 'LineWidth', 2, 'Color','#127B2E')
xlim([0.5 2.5]); ylim([0 2500]); ylabel('IPSC Amplitude (pA)'); title('vip');
set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);
xticklabels({'','ACSF','', 'CNO'}); xlabel('IPSC amplitude (pA)');


%%Plot TC data for LEC, VR Data
figure
yline(100, '--', 'Color', [.7 .7 .7],'LineWidth',2);
hold on; errorbar([-5:1:15], lecCCK.tc(6:end), lecCCK.semTc(6:end), '.', 'MarkerSize', 50, 'Color','#127B2E', 'LineWidth', 3);
xlim([-5 15]); ylim([0 120]);title('cck'); xlabel('time (min)'); ylabel('IPSC amplitude (%)');
set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);

figure
yline(100, '--', 'Color', [.7 .7 .7],'LineWidth',2);
hold on; errorbar([-5:1:16], lecPV.tc(25:46), lecPV.semTc(25:46), '.', 'MarkerSize', 50, 'Color','#127B2E', 'LineWidth', 3);
xlim([-5 16]); ylim([0 120]);title('pv'); xlabel('time (min)'); ylabel('IPSC amplitude (%)');
set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);


figure
yline(100, '--', 'Color', [.7 .7 .7],'LineWidth',2);
hold on; errorbar([-5:1:21], lecSST.tc(12:38), lecSST.semTc(12:38), '.', 'MarkerSize', 50, 'Color','#127B2E', 'LineWidth', 3);
xlim([-5 16]); ylim([0 120]);title('sst'); xlabel('time (min)'); ylabel('IPSC amplitude (%)');
set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);

figure
yline(100, '--', 'Color', [.7 .7 .7],'LineWidth',2);
hold on; errorbar([-5:1:14], lecVIP.tc(6:end), lecVIP.semTc(6:end), '.', 'MarkerSize', 50, 'Color','#127B2E', 'LineWidth', 3);
xlim([-5 16]); ylim([0 120]);title('vip'); xlabel('time (min)'); ylabel('IPSC amplitude (%)');
set(gca,'FontSize', 18,'FontWeight', 'bold', 'LineWidth',3);

%% Trash Code

%contains(string(tc(1,:)),'_rc')


% 
% %load files if you are loading new animals into the tc structure
% in = readcell('in'); in_type = cell2mat(in(contains(string(in), 'in'),2));
% pathway = readcell('pathway');
% cd baseline; files = dir(['amp_baseline*csv']);
% tc = readcell(files.name);
% 
% %get all amp values for blue light stim, note have to have the tag
% %"baseline_b" if used something else change tag
% idx = length(vip_tc.amp)
% vip_tc.amp(idx) = {tc(2:end, contains(string(tc(1,:)), 'baseline_b'))}
% 
% 
% %store cell ID
% vip_tc.id(idx) = string(extractAfter(files.name, 'KO'));
% 
% %store IN type, this is for sanity check to make sure code is doing what
% %you want it to
% vip_tc.in(idx) = in_type;
% 
% %store sweep with CNO app
% vip_tc.cno(idx) = cell2mat(in(contains(string(in), 'cno'),2));
% 
% time = -(vip_tc.cno(idx)-1):1:(length(vip_tc.amp{idx})-vip_tc.cno(idx))
% 
% 
% 
% 
% 
% abs amplitude
% 
% save(['vip_tc'], 'vip_tc')
