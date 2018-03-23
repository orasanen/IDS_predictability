%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 2: Analyze results (this is only for main effects in case of varying
% distributions of IDS in the training data). Use idsads_analyze_results.m 
% for the basic analysis with fixed IDS/ADS proportion.

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
plot_basicf0 = 0;

% Define result file
filename = '/results/varying_IDS_proportion_results_for_the_revised_submission/results_20-Mar-2018 08:57:41_ManyBabies_usesyllables1_framesize_100_proportions.mat';

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

PROPDAT = cell(size(F0prob_MOMC,2),1);

for propiter = 1:size(F0prob_MOMC,2)
    
    
    F0prob = cell(size(F0prob_MOMC,1),1);
    if(strcmp(predmethod,'LSTM'))
        for k = 1:length(F0prob)
        F0prob{k} = F0prob_LSTM{k,propiter}; %Replace results with LSTM?
        end
    elseif(strcmp(predmethod,'MOCM'))
        for k = 1:length(F0prob)
            F0prob{k} = F0prob_MOMC{k,propiter}; % Replace results with LSTM?
        end
    elseif(strcmp(predmethod,'both'))
        F0prob = cell(length(F0prob_LSTM),1);
        for k = 1:length(F0prob)
            F0prob{k} = (F0prob_LSTM{k,propiter}+F0prob_MOMC{k,propiter})/2;
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
            am_data(n) = sum(number_of_frames(i2),1);
            am_data(n+1) = sum(number_of_frames(i3),1);
            
            n = n+2;
        end
    end
    
    SUBJ_STATS = SUBJ_STATS(1:n-1,:);
    am_data = am_data(1:n-1,:);
        
    PROPDAT{propiter} = SUBJ_STATS;
    
end

id_proportions = 0:0.25:1;
Nprops = length(id_proportions);

propmeans = zeros(Nprops,2);
propmean_devs = zeros(Nprops,2);

propdev = zeros(Nprops,2);
propdev_devs = zeros(Nprops,2);

propmax = zeros(Nprops,2);
propmax_devs = zeros(Nprops,2);

for prop = 1:Nprops   
    a = find(PROPDAT{prop}(:,1) == 1);
    propmeans(prop,1) = mean(PROPDAT{prop}(a,3));    
    propmean_devs(prop,1) = std(PROPDAT{prop}(a,3))./sqrt(length(a));  
    
    propdev(prop,1) = mean(PROPDAT{prop}(a,2));
    propdev_devs(prop,1) = std(PROPDAT{prop}(a,2))./sqrt(length(a));  
    
    propmax(prop,1) = mean(PROPDAT{prop}(a,4));
    propmax_devs(prop,1) = std(PROPDAT{prop}(a,4))./sqrt(length(a));  
    
    b = find(PROPDAT{prop}(:,1) == 2);
    propmeans(prop,2) = mean(PROPDAT{prop}(b,3));        
    propmean_devs(prop,2) = std(PROPDAT{prop}(b,3))./sqrt(length(b));    
    
    propdev(prop,2) = mean(PROPDAT{prop}(b,2));        
    propdev_devs(prop,2) = std(PROPDAT{prop}(b,2))./sqrt(length(b));    
    
    propmax(prop,2) = mean(PROPDAT{prop}(b,4));
    propmax_devs(prop,2) = std(PROPDAT{prop}(b,4))./sqrt(length(b));  
        
end

h = figure(1);clf;
subplot(1,3,1);
bar(fliplr(propmeans));
grid;

set(gca,'XTickLabel',id_proportions*100)
xlabel('proportion of IDS in training data (%)');
ylabel('mean predictability');

drawstds(h,(1:Nprops)-0.15,propmeans(:,2),propmean_devs(:,2),0.1,2,'red');
drawstds(h,(1:Nprops)+0.15,propmeans(:,1),propmean_devs(:,1),0.1,2,'red');

legend({'IDS','ADS'},'Location','NorthEast')

subplot(1,3,2);
bar(fliplr(propdev));
grid;

set(gca,'XTickLabel',id_proportions*100)
xlabel('proportion of IDS in training data (%)');
ylabel('SD of predictability');

drawstds(h,(1:Nprops)-0.15,propdev(:,2),propdev_devs(:,2),0.1,2,'red');
drawstds(h,(1:Nprops)+0.15,propdev(:,1),propdev_devs(:,1),0.1,2,'red');
legend({'IDS','ADS'},'Location','NorthEast')
subplot(1,3,3);
bar(fliplr(propmax));
grid;

set(gca,'XTickLabel',id_proportions*100)
xlabel('proportion of IDS in training data (%)');
ylabel('max predictability');

drawstds(h,(1:Nprops)-0.15,propmax(:,2),propmax_devs(:,2),0.1,2,'red');
drawstds(h,(1:Nprops)+0.15,propmax(:,1),propmax_devs(:,1),0.1,2,'red');
colormap([0.2 0.7 0.9;0.9 0.9 0.9])
legend({'IDS','ADS'},'Location','NorthEast')



