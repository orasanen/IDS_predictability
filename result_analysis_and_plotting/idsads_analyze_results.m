%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2: Analyze results

% repeat the same analysis for MOCM and LSTM models

clear all

predmethod = 'MOCM'; % which method to analyze (MOCM, LSTM or both (average)).

% Compare individual utterances instead of subject means? default = 0
% note: statistical analyses for individual utterances are not properly
% conducted, as there is no correction for repeated measures per subject. 
% use only pooledstats = 0 for proper comparison.
set(0,'DefaultAxesFontSize',16);
set(0,'DefaultAxesFontSize',16);
pooledstats = 0; 
plot_basicf0 = 0; % show also raw F0 tests?

% Define result file here
 
filename = 'results/results_for_the_first_manuscript_submission/results_20-Oct-2017 00:04:48_ManyBabies_usesyllables1_framesize_100.mat'; % syllabic-frame
%filename = 'results/results_for_the_first_manuscript_submission/results_19-Oct-2017 21:07:23_ManyBabies_usesyllables0_framesize_100.mat' % fixed-frame

curdir = fileparts(which('predictability_IDS_ADS_main'));

load(filename);

if(contains(filename,'usesyllables1'))
    usesyllables = 1;
else
    usesyllables = 0;
end
tmp1 = strfind(filename,'framesize');
tmp2 = strfind(filename,'.mat');
framesize = str2num(filename(tmp1+10:tmp2-1))./1000;

% Pick results from the correct model

if(strcmp(predmethod,'LSTM'))
    F0prob = F0prob_LSTM; %Replace results with LSTM?
elseif(strcmp(predmethod,'MOCM'))
    F0prob = F0prob_MOMC; % Replace results with LSTM?
elseif(strcmp(predmethod,'both'))
    F0prob = cell(length(F0prob_LSTM),1);
    for k = 1:length(F0prob)
        F0prob{k} = (F0prob_LSTM{k}+F0prob_MOMC{k})/2;
    end
else
    error('Unknown method for sequence prediction.');
end

%% ------------------------------------------------------------------------
% Part 1: Main effects and storage to .csv files for further analysis
%--------------------------------------------------------------------------

% Actual durations of the utterances
if(usesyllables)
    utterance_length = cellfun(@max,bounds_orig_syllable_t)-cellfun(@min,bounds_orig_syllable_t);
else
    utterance_length = cellfun(@max,bounds_t)-cellfun(@min,bounds_t);
end
    

number_of_frames = cellfun(@length,F0prob);


% Indices to IDS and ADS utterances
ids_i = cellfind(METADATA(:,1),'IDS');
ads_i = cellfind(METADATA(:,1),'ADS');

% Print basic stats
fprintf('Total %d utterances (%d IDS, %d ADS).\n',length(F0prob),length(ids_i),length(ads_i));
fprintf('Total syllables: %d.\n',sum(cellfun(@length,bounds_orig_syllable_t)-1));
fprintf('Mean IDS duration: %0.2f (+- %0.2f) s.\n',mean(utterance_length(ids_i)),std(utterance_length(ids_i)));
fprintf('Mean ADS duration: %0.2f (+- %0.2f) s.\n',mean(utterance_length(ads_i)),std(utterance_length(ads_i)));

% probs
all_SD = cellfun(@nanstd,F0prob);
all_mean = cellfun(@nanmean,F0prob);
all_max = cellfun(@nanmax,F0prob);
all_min = cellfun(@nanmin,F0prob);
all_maxmin = cellfun(@nanmax,F0prob)-cellfun(@nanmin,F0prob);

% raw F0 in Hz

tmp = F0_raw_orig;
for k = 1:length(tmp)
   tmp{k}(tmp{k} == 0) = NaN;
end

all_SD_orig = cellfun(@nanstd,tmp);
all_mean_orig = cellfun(@nanmean,tmp);
all_min_orig = cellfun(@nanmin,tmp);
all_max_orig = cellfun(@nanmax,tmp);
all_maxmin_orig = cellfun(@nanmax,tmp)-cellfun(@nanmin,tmp);

style_labels = zeros(length(utterance_length),1);
style_labels(ads_i) = 1;
style_labels(ids_i) = 2;

subject_labels = str2num(strvcat(METADATA(:,2)));

% Create a matrix of utterance level features
UTT_STATS = [style_labels all_SD all_mean all_max all_min all_maxmin all_SD_orig all_mean_orig all_max_orig all_min_orig all_maxmin_orig log(utterance_length) subject_labels];

% Get subject level average descriptors

uq_subjects = unique(subject_labels);
SUBJ_STATS = zeros(length(uq_subjects)*2,size(UTT_STATS,2));
am_data = zeros(length(uq_subjects).*2,1);
n_utt_per_talker = zeros(length(uq_subjects),1);
n = 1;
for s = 1:length(uq_subjects)
    i1 = find(subject_labels == uq_subjects(s));
    n_utt_per_talker(s) = length(i1);
    i2 = intersect(i1,ads_i);
    i3 = intersect(i1,ids_i);
    
    if(length(i2) > 0 && length(i3) > 0)
       
        SUBJ_STATS(n,:) = mean(UTT_STATS(i2,:),1);
        SUBJ_STATS(n+1,:) = mean(UTT_STATS(i3,:),1);
        %am_data(n) = sum(utterance_length(i2));
        %am_data(n+1) = sum(utterance_length(i3));
        am_data(n) = sum(number_of_frames(i2),1);
        am_data(n+1) = sum(number_of_frames(i3),1);
                
        n = n+2;
    end
end

SUBJ_STATS = SUBJ_STATS(1:n-1,:);
am_data = am_data(1:n-1,:);

fprintf('Average %0.2f (+-%0.2f) utterances pe talker.\n',mean(n_utt_per_talker),std(n_utt_per_talker));

SUBJ_STATS = [SUBJ_STATS am_data]; % append with amount of data per speaker and style

% Apepend utterance-level stats also with the amount of samples behind each
% probability estimate (on average)
am_data_full = zeros(size(UTT_STATS,1),1);
for k = 1:size(UTT_STATS,1)
    i1 = find(SUBJ_STATS(:,1) == UTT_STATS(k,1)); % matching style
    i2 = find(SUBJ_STATS(:,end-1) == UTT_STATS(k,end)); % matching speaker
    i = intersect(i1,i2);    
    am_data_full(k) = SUBJ_STATS(i,end);    
end
UTT_STATS = [UTT_STATS am_data_full];

% Write utterance-level and subject-level stats to .csv files

% utterance stats
headers = {'IDSADS','prob_var','prob_mean','prob_max','prob_min','prob_maxmin','F0_var','F0_mean','F0_max','F0_min','F0_maxmin','log_dur','babyID','amount_of_data'};
commaHeader = [headers;repmat({','},1,numel(headers))]; %insert commaas
commaHeader = commaHeader(:)';
textHeader = cell2mat(commaHeader); %cHeader in text with commas
csv_filename = sprintf('%s/results/results_utterances_%s_%s_usesyllables%d_framesize_%d_%s.csv',curdir,datestr(now),dataset,usesyllables,framesize*1000,predmethod);

fid = fopen(csv_filename,'w'); % write header first
fprintf(fid,'%s\n',textHeader) ;
fclose(fid);
dlmwrite(csv_filename,UTT_STATS,'-append'); % append data matrix

% subject stats
headers = {'IDSADS','prob_var','prob_mean','prob_max','prob_min','prob_maxmin','F0_var','F0_mean','F0_max','F0_min','F0_maxmin','log_dur','babyID','amount_of_data'};
commaHeader = [headers;repmat({','},1,numel(headers))]; %insert commaas
commaHeader = commaHeader(:)';
textHeader = cell2mat(commaHeader); %cHeader in text with commas

csv_filename = sprintf('%s/results/results_subjs_%s_%s_usesyllables%d_framesize_%d_%s.csv',curdir,datestr(now),dataset,usesyllables,framesize*1000,predmethod);

fid = fopen(csv_filename,'w');
fprintf(fid,'%s\n',textHeader);
fclose(fid);
dlmwrite(csv_filename,SUBJ_STATS,'-append');

% Normalize likelihoods by removing the predicted values from the amount of
% matching training data

SPSS_before_norm = UTT_STATS;
D_before_norm = SUBJ_STATS;

for f = 2:6
    amount_of_data = UTT_STATS(:,end);    
    coeffs = polyfit(amount_of_data,UTT_STATS(:,f),1);
    prediction = coeffs(1)*amount_of_data+coeffs(2);
    UTT_STATS(:,f) = UTT_STATS(:,f)-prediction;
    
    amount_of_data = SUBJ_STATS(:,end);    
    coeffs = polyfit(amount_of_data,SUBJ_STATS(:,f),1);
    prediction = coeffs(1)*amount_of_data+coeffs(2);
    SUBJ_STATS(:,f) = SUBJ_STATS(:,f)-prediction;    
end

% Plot results


if(plot_basicf0)
    imaghandle = figure('Position',[200 200 657 558]);clf;
    subplot(3,1,1);
else
    imaghandle = figure('Position',[200 200 657 558*2/3]);clf;
    subplot(2,1,1);
end

descriptor_plot_order = [3 2 4 5 6]; % the order to plot descriptors from result matrices 

if(~pooledstats)
    bardata = [mean(D_before_norm(2:2:end,descriptor_plot_order))' mean(D_before_norm(1:2:end,descriptor_plot_order))'];
    bardevs = [std(D_before_norm(2:2:end,descriptor_plot_order))' std(D_before_norm(1:2:end,descriptor_plot_order))'];
    bardevs = bardevs./sqrt(length(uq_subjects)); % Get standard errors
else
    bardata = [mean(SPSS_before_norm(ids_i,descriptor_plot_order))' mean(SPSS_before_norm(ads_i,descriptor_plot_order))'];
    bardevs = [std(SPSS_before_norm(ids_i,descriptor_plot_order))' std(SPSS_before_norm(ads_i,descriptor_plot_order))'];
    bardevs(:,1) = bardevs(:,1)./sqrt(length(ids_i));
    bardevs(:,2) = bardevs(:,2)./sqrt(length(ads_i));
end
bar(bardata);

drawstds(imaghandle,1:5,bardata(:,1),bardevs(:,1),0.095,2,'red',-0.1375);
drawstds(imaghandle,1:5,bardata(:,2),bardevs(:,2),0.095,2,'red',+0.1375);
grid;
set(gca,'XTickLabel',{'mean','SD','max','min','range'});
title('F0 predictability P(q_s|q_s_-_1, ..., q_s_-_m)');
ylabel('likelihood');

if(plot_basicf0)
    subplot(3,1,2);
else
    subplot(2,1,2);
end

if(~pooledstats)
    
    bardata = [mean(SUBJ_STATS(2:2:end,descriptor_plot_order))' mean(SUBJ_STATS(1:2:end,descriptor_plot_order))'];
    bardevs = [std(SUBJ_STATS(2:2:end,descriptor_plot_order))' std(SUBJ_STATS(1:2:end,descriptor_plot_order))'];
    bardevs = bardevs./sqrt(length(uq_subjects)); % Get standard errors
else
    bardata = [mean(UTT_STATS(ids_i,descriptor_plot_order))' mean(UTT_STATS(ads_i,descriptor_plot_order))'];
    bardevs = [std(UTT_STATS(ids_i,descriptor_plot_order))' std(UTT_STATS(ads_i,descriptor_plot_order))'];
    bardevs(:,1) = bardevs(:,1)./sqrt(length(ids_i));
    bardevs(:,2) = bardevs(:,2)./sqrt(length(ads_i));
end
bar(bardata);

drawstds(imaghandle,1:5,bardata(:,1),bardevs(:,1),0.095,2,'red',-0.1375);
drawstds(imaghandle,1:5,bardata(:,2),bardevs(:,2),0.095,2,'red',+0.1375);
grid;
set(gca,'XTickLabel',{'mean','SD','max','min','range'});
%leg2 = legend({'IDS','ADS'},'Position',[0.8410 0.5808 0.1141 0.0653]);
title('F0 predictability after controlling for training data');
ylabel('normalized likelihood');


if(plot_basicf0)
subplot(3,1,3);


if(~pooledstats)
    bardata = [mean(SUBJ_STATS(2:2:end,descriptor_plot_order+5))' mean(SUBJ_STATS(1:2:end,descriptor_plot_order+5))'];
    bardevs = [std(SUBJ_STATS(2:2:end,descriptor_plot_order+5))' std(SUBJ_STATS(1:2:end,descriptor_plot_order+5))'];
    bardevs = bardevs./sqrt(length(uq_subjects));
else
    bardata = [mean(UTT_STATS(ids_i,descriptor_plot_order+5))' mean(UTT_STATS(ads_i,descriptor_plot_order+5))'];
    bardevs = [std(UTT_STATS(ids_i,descriptor_plot_order+5))' std(UTT_STATS(ads_i,descriptor_plot_order+5))'];
    bardevs(:,1) = bardevs(:,1)./sqrt(length(ids_i));
    bardevs(:,2) = bardevs(:,2)./sqrt(length(ads_i));
end

bar(bardata);

drawstds(imaghandle,1:5,bardata(:,1),bardevs(:,1),0.095,2,'red',-0.1375);
drawstds(imaghandle,1:5,bardata(:,2),bardevs(:,2),0.095,2,'red',+0.1375);
grid;
set(gca,'XTickLabel',{'F0 mean','F0 SD','F0 max','F0 min','F0 max-min'});
title('standard F0 measures');
ylabel('Hz');
end

colormap([0.2 0.7 0.9;0.9 0.9 0.9])

% Swap the order of all features for statistical testing
orderi = [3 2 4 5 6 8 7 9 10 11];

% Do pair-wise ttests

siglevel = 0.05/10; % Bonferroni normalized (not relevant, see below)
x = 1;
h = zeros(10,1);
p = zeros(10,1);
cohens_d = zeros(10,1);
stats = cell(10,1);
tstat = zeros(10,1);
means = zeros(10,2);
devs = zeros(10,2);

for f = orderi
    if(pooledstats)
        %i = find(UTT_STATS(:,1) == 0);
        x_ads = UTT_STATS(ads_i,f);
        i = find(UTT_STATS(:,1) == 1);
        x_ids = UTT_STATS(ids_i,f);
        [h(x),p(x),~,stats{x}] = ttest2(x_ads,x_ids,siglevel);
    else
        x_ads = SUBJ_STATS(1:2:end,f); % ADS values
        x_ids = SUBJ_STATS(2:2:end,f); % IDS values
        [h(x),p(x),~,stats{x}] = ttest(x_ads,x_ids,siglevel);
                
        %pooled_sd = sqrt(0.5.*(var(x_ads)+var(x_ids)));
        %cohens_d(x) = abs(stats{x}.tstat./sqrt(length(uq_subjects)));
        cohens_d(x) = abs(2*stats{x}.tstat)/sqrt(stats{x}.df);
        
        means(x,1) = mean(x_ids);
        means(x,2) = mean(x_ads);
        devs(x,1) = std(x_ids);
        devs(x,2) = std(x_ads);
        tstat(x) = stats{x}.tstat;
        
        
    end
    x = x+1;
end

% Correct significances to Holm-Bonferroni with the given p-values

[siglevels,h] = holmBonferroni(p,0.05);

% Update plots with markers for significant differences

if(plot_basicf0)
    subplot(3,1,2);
else
    subplot(2,1,2);
end
    

for k = 1:5
    if(h(k)) % check if significant after correction
        H = sigstar([k-0.125,k+0.125],0.03);
        
        v = get(H(1));
        text(k,mean(v.YData)+0.023,sprintf('t = %0.2f',stats{k}.tstat),'HorizontalAlignment','center','FontSize',13);
    end
end
ylim([-0.1 0.13])
if(plot_basicf0)
subplot(3,1,3);
for k = 1:5
    if(h(k+5))
        H = sigstar([k-0.125,k+0.125],0.03);
        
        v = get(H(1));
        text(k,mean(v.YData)+45,sprintf('t = %0.2f',stats{k+5}.tstat),'HorizontalAlignment','center','FontSize',13);
    end
end
end

if(plot_basicf0)
    subplot(3,1,1);
    leg1 = legend({'IDS','ADS'},'Position',[0.8363 0.8733 0.1141 0.0653]);
    subplot(3,1,2);
    leg2 = legend({'IDS','ADS'},'Position',[0.8410 0.5808 0.1141 0.0653]);
    subplot(3,1,3);
    leg3 = legend({'IDS','ADS'},'Position',[0.8427 0.2891 0.1141 0.0653]);
else
    subplot(2,1,1);
    leg1 = legend({'IDS','ADS'},'Position',[0.8363 0.8733 0.1141 0.0653]);
    subplot(2,1,2);
    leg2 = legend({'IDS','ADS'},'Position',[0.8410 0.4 0.1141 0.0653]);   
end
% End of plotting

%% ------------------------------------------------------------------------
% Part 2: Compare sentence styles
%--------------------------------------------------------------------------

nl_i = cellfind(METADATA(:,4),'N/A','exact');
fam_i = cellfind(METADATA(:,4),'familiar','exact');
unfam_i = cellfind(METADATA(:,4),'unfamiliar','exact');

[anno,annofilenames] = parseManyBabiesAnnotation();

descriptors = {'prob_mean','prob_var','F0_mean','F0_var'};

subject_column = cellfind(headers,'babyID');

familiar_measures = zeros(length(uq_subjects),length(descriptors));
unfamiliar_measures = zeros(length(uq_subjects),length(descriptors));
nolabel_measures = zeros(length(uq_subjects),length(descriptors));

h = zeros(length(descriptors),1);
p = zeros(length(descriptors),1);
cohens_d = zeros(length(descriptors),1);
dd_means = zeros(length(descriptors),2);
dd_devs = zeros(length(descriptors),2);
stat = cell(length(descriptors),1);



for dd = 1:length(descriptors)
   f = cellfind(headers,descriptors{dd},'exact');
      
   for s = 1:length(uq_subjects)       
       tmp1 = find(SPSS_before_norm(:,subject_column) == uq_subjects(s));
       tmp2 = find(SPSS_before_norm(:,1) == 2); % IDS only
       tmp3 = intersect(tmp1,tmp2); % correct speaker and style
       
       utts_familiar = intersect(tmp3,fam_i);
       utts_unfamiliar = intersect(tmp3,unfam_i);
       utts_nolabel = intersect(tmp3,nl_i);
       
       familiar_measures(s,dd) = mean(SPSS_before_norm(utts_familiar,f));
       unfamiliar_measures(s,dd) = mean(SPSS_before_norm(utts_unfamiliar,f));       
       nolabel_measures(s,dd) = mean(SPSS_before_norm(utts_nolabel,f));       
   end        
   
   % Do stats   
   [h(dd),p(dd),~,stat{dd}] = ttest(familiar_measures(:,dd),unfamiliar_measures(:,dd),0.05);   
   cohens_d(dd) = abs(2*stat{dd}.tstat)/sqrt(stat{dd}.df);   
   dd_means(dd,1) = mean(familiar_measures(:,dd));
   dd_means(dd,2) = mean(unfamiliar_measures(:,dd));   
   dd_devs(dd,1) = std(familiar_measures(:,dd));
   dd_devs(dd,2) = std(unfamiliar_measures(:,dd));   
end

% Normalize for multiple comparisons
[siglevels,h] = holmBonferroni(p,0.05);

tmp = find(h); % anything significant?
if(~isempty(tmp))
    fprintf('Significant differences between familiar and unfamiliar labeling sentences:\n');
    for j = 1:length(tmp)
        fprintf('*%s (p = %0.4f, t = %0.2f, cohen''s d = %0.2f). ',descriptors{tmp(j)},p(tmp(j)),stat{tmp(j)}.tstat,cohens_d(tmp(j)));        
        fprintf('Familiar: %0.2f (+- %0.2f), unfamiliar: %0.2f (+- %0.2f).\n.',dd_means(tmp(j),1),dd_devs(tmp(j),1),dd_means(tmp(j),2),dd_devs(tmp(j),2));
    end
else
    fprintf('Utterance category differences are not significant.\n');
end

%% ------------------------------------------------------------------------
% Part 2.5: Compare sentence styles only for keywords
%--------------------------------------------------------------------------

nl_i = cellfind(METADATA(:,4),'N/A','exact');
fam_i = cellfind(METADATA(:,4),'familiar','exact');
unfam_i = cellfind(METADATA(:,4),'unfamiliar','exact');

[anno,annofilenames] = parseManyBabiesAnnotation();

%descriptors = {'prob_mean','prob_var','F0_mean','F0_var'};
descriptors = {'mean','std','max'};

subject_column = cellfind(headers,'babyID');

familiar_measures = zeros(length(uq_subjects),length(descriptors));
unfamiliar_measures = zeros(length(uq_subjects),length(descriptors));
familiar_measures_raw = zeros(length(uq_subjects),length(descriptors));
unfamiliar_measures_raw = zeros(length(uq_subjects),length(descriptors));
nolabel_measures = zeros(length(uq_subjects),length(descriptors));

h = zeros(length(descriptors),1);
p = zeros(length(descriptors),1);
cohens_d = zeros(length(descriptors),1);
dd_means = zeros(length(descriptors),2);
dd_devs = zeros(length(descriptors),2);
stat = cell(length(descriptors),1);

for dd = 1:length(descriptors)
    f = cellfind(headers,descriptors{dd},'exact');
    
    for s = 1:length(uq_subjects)
        tmp1 = find(SPSS_before_norm(:,subject_column) == uq_subjects(s));
        tmp2 = find(SPSS_before_norm(:,1) == 2); % IDS only
        tmp3 = intersect(tmp1,tmp2); % correct speaker and style
        
        utts_familiar = intersect(tmp3,fam_i);
        utts_unfamiliar = intersect(tmp3,unfam_i);
        utts_nolabel = intersect(tmp3,nl_i);
        
        cc = 1;
        for j = 1:length(utts_familiar)
            [~,dah,~] = fileparts(filenames{utts_familiar(j)});
            tmp = cellfind(annofilenames,dah);
            
            o1 = anno{tmp}.kw_onset_time(1);
            o2 = anno{tmp}.kw_offset_time(1);
            bb = bounds_t{utts_familiar(j)};
            
            first_frame = max(1,find(bb > o1,1)-1);
            if(isempty(first_frame))
                first_frame = 1;
            end
            last_frame = min(find(bb > o2,1),length(F0prob{utts_familiar(j)}));
            if(isempty(last_frame))
                last_frame = length(bb)-1;
            end
            
            y = F0prob{utts_familiar(j)}(first_frame:last_frame);
            familiar_measures(s,dd) = familiar_measures(s,dd)+eval(sprintf('%s(y)',descriptors{dd}))./length(utts_familiar);
            y = F0_raw_orig{utts_familiar(j)}(round(100*o1:100*o2));
            % Listen to make sure all is alright
            % [x,fs] = audioread(filenames{utts_familiar(j)});
            %soundsc(x(round(o1*fs:o2*fs)),fs);
            %pause;
            y(y == 0) = [];
            if(~isempty(y))
                familiar_measures_raw(s,dd) = familiar_measures_raw(s,dd)+eval(sprintf('%s(y)',descriptors{dd}));
                cc = cc+1;
            end
            
        end
        familiar_measures_raw(s,dd) = familiar_measures_raw(s,dd)./cc;
        
        
        cc = 1;
        for j = 1:length(utts_unfamiliar)
            [~,dah,~] = fileparts(filenames{utts_unfamiliar(j)});
            tmp = cellfind(annofilenames,dah);
            
            o1 = anno{tmp}.kw_onset_time(1);
            o2 = anno{tmp}.kw_offset_time(1);
            bb = bounds_t{utts_unfamiliar(j)};
            
            first_frame = max(1,find(bb > o1,1)-1);
            if(isempty(first_frame))
                first_frame = 1;
            end
            last_frame = min(find(bb > o2,1),length(F0prob{utts_unfamiliar(j)}));
            if(isempty(last_frame))
                last_frame = length(bb)-1;
            end
            
            y = F0prob{utts_unfamiliar(j)}(first_frame:last_frame);
            unfamiliar_measures(s,dd) = unfamiliar_measures(s,dd)+eval(sprintf('%s(y)',descriptors{dd}))./length(utts_unfamiliar);
            y = F0_raw_orig{utts_unfamiliar(j)}(round(100*o1:min(100*o2,length(F0_raw_orig{utts_unfamiliar(j)}))));
            
            % Make sure everything is right
            
            %[x,fs] = audioread(filenames{utts_unfamiliar(j)});
            %soundsc(x(round(o1*fs:o2*fs)),fs);
            %pause;
            %figure(6);plot(y);pause;
            y(y == 0) = [];
            if(~isempty(y))
                unfamiliar_measures_raw(s,dd) = unfamiliar_measures_raw(s,dd)+eval(sprintf('%s(y)',descriptors{dd}));
                cc = cc+1;
            end
        end
        unfamiliar_measures_raw(s,dd) = unfamiliar_measures_raw(s,dd)./cc;        
    end
end

h = [];
p = [];
testname = {};
for k = 1:length(descriptors)
    [h(k),p(k)] = ttest(familiar_measures(:,k),unfamiliar_measures(:,k));    
    testname{k} = sprintf('likelihood %s',descriptors{k});
end
for k = 1:length(descriptors)
    [h(k+length(descriptors)),p(k+length(descriptors))] = ttest(familiar_measures_raw(:,k),unfamiliar_measures_raw(:,k));    
    testname{k+length(descriptors)} = sprintf('F0 %s',descriptors{k});
end



% Normalize for multiple comparisons
[siglevels,h] = holmBonferroni(p,0.05);

tmp = find(h); % anything significant?
if(~isempty(tmp))
    fprintf('Significant differences between familiar and unfamiliar labeling sentences in IDS:\n');
else
    fprintf('No significant differences between familiar and unfamiliar labeling sentences in IDS.\n')
end
   


%% ------------------------------------------------------------------------
% Part 3: Compare against style typicality ratings (from https://osf.io/re95x/)
%--------------------------------------------------------------------------

load saves/ADIDratings_MT.mat

% Get ratings and pair them with current data
rating = ones(N,1).*NaN;
for k = 1:length(file)
    ff = file{k};
    i = cellfind(filenames,ff);
    rating(i) = idsness(k);
end

a = find(~isnan(rating));
a = intersect(union(ids_i,ads_i),a);

p = [];
fprintf('Correlations w.r.t. human ratings on IDS-likeness:\n');
for j = 2:12
    [rho,p(j-1)] = corr(SPSS_before_norm(a,j),rating(a),'type','Spearman');
    fprintf('%s: r = %0.3f (p = %0.3f).\n',headers{j},rho,p(j-1));
end

[siglevel,h] = holmBonferroni(p,0.05);

DAT = [SPSS_before_norm(a,:) rating(a)];

headers = {'IDSADS','prob_var','prob_mean','prob_max','prob_min','prob_maxmin','F0_var','F0_mean','F0_max','F0_min','F0_maxmin','log_dur','babyID','amount_of_data','human_rating'};
commaHeader = [headers;repmat({','},1,numel(headers))]; %insert commaas
commaHeader = commaHeader(:)';
textHeader = cell2mat(commaHeader); %cHeader in text with commas

csv_filename = sprintf('%s/results/results_ratings_%s_%s_usesyllables%d_framesize_%d_%s.csv',curdir,datestr(now),dataset,usesyllables,framesize*1000,predmethod);

fid = fopen(csv_filename,'w');
fprintf(fid,'%s\n',textHeader);
fclose(fid);
dlmwrite(csv_filename,DAT,'-append');

