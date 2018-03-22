% Do temporal analyses reported in the first version of the paper

filename = 'results/results_for_the_first_manuscript_submission/results_20-Oct-2017 00:04:48_ManyBabies_usesyllables1_framesize_100.mat'; % syllabic-frame

load(filename);

predmethod = 'MOCM';

if(strcmp(predmethod,'LSTM'))
    F0prob = F0prob_LSTM;
elseif(strcmp(predmethod,'MOCM'))
    F0prob = F0prob_MOMC;
elseif(strcmp(predmethod,'both'))
    F0prob = cell(length(F0prob_LSTM),1);
    for k = 1:length(F0prob)
        F0prob{k} = (F0prob_LSTM{k}+F0prob_MOMC{k})/2;
    end
else
    error('Unknown method for sequence prediction.');
end


%% -----------------------------------------------------------------------
% Part 1, compare long and short IDS and ADS utterances as a function of
% the relative temporal position in the utterances.
%-------------------------------------------------------------------------

for short = 0:1
    
    % number of syllables per utterance
    nsyls = cellfun(@length,bounds_orig_syllable_t)-1;
    
    % pointers to IDS and ADS data
    ids_i = cellfind(METADATA(:,1),'IDS');
    ads_i = cellfind(METADATA(:,1),'ADS');
        
    subject_labels = str2num(strvcat(METADATA(:,2)));
    
    % Get subject means for all descriptors
    uq_subjects = unique(subject_labels);
    
    % Define length categories
    if(short)
        minsyls = 0;
        maxsyls = 11;
    else
        minsyls = 12;
        maxsyls = Inf;
    end
    
    ids_i_selected = intersect(intersect(ids_i,find(nsyls >= minsyls)),find(nsyls <= maxsyls));
    ads_i_selected = intersect(intersect(ads_i,find(nsyls >= minsyls)),find(nsyls <= maxsyls));
            
    mm_ads = mean(nsyls(ads_i_selected));
    mm_ids = mean(nsyls(ids_i_selected));
        
    % Get average F0 and likelihood trajectories as a function of position
    % in the utterances
    
    fprintf('%d IDS, %d ADS that satisfy length constraints.\n',length(ids_i_selected),length(ads_i_selected));
    
    % how many relative positions to analyze
    if(short)
        target_length = 10;
    else
        target_length = 10;
    end
        
    cumprob_ids = zeros(length(uq_subjects),target_length);
    cumprob_ads = zeros(length(uq_subjects),target_length);
    cumf0_ids = zeros(length(uq_subjects),target_length);
    cumf0_ads = zeros(length(uq_subjects),target_length);
        
    % IDS
    for k = 1:length(uq_subjects)
        i = find(subject_labels == uq_subjects(k)); % k:th subject
        i2 = intersect(i,ids_i_selected); % IDS utterances that match length and subject
                 
        for j = 1:length(i2) % go through all matching utterances
            
            tmp = resample_by_interpolation(F0prob{i2(j)},target_length); % resample to fixed length
                                                
            % get corresponding F0 contours (from segment boundary information)
            bb = round(bounds_t{i2(j)}.*100);            
            tmp2 = F0_raw{i2(j)}(max(1,bb(1)):min(bb(end),length(F0_raw{i2(j)})));            
            tmp2 = resample_by_interpolation(tmp2,target_length);
            
            cumprob_ids(k,:) = cumprob_ids(k,:)+(tmp./length(i2));
            cumf0_ids(k,:) = cumf0_ids(k,:)+(tmp2./length(i2));
        end
    end
       
    % ADS (same as above for ADS utterances)
    for k = 1:length(uq_subjects)
        i = find(subject_labels == uq_subjects(k));
        i2 = intersect(i,ads_i_selected);
        for j = 1:length(i2)
            tmp = resample_by_interpolation(F0prob{i2(j)},target_length);
            
            bb = round(bounds_t{i2(j)}.*100);
            tmp2 = F0_raw{i2(j)}(max(1,bb(1)):min(bb(end),length(F0_raw{i2(j)})));            
            tmp2 = resample_by_interpolation(tmp2,target_length);
            
            cumprob_ads(k,:) = cumprob_ads(k,:)+(tmp./length(i2));
            cumf0_ads(k,:) = cumf0_ads(k,:)+(tmp2./length(i2));
        end
    end
    
    % Plot results
    
    h = figure;clf;
    subplot(2,1,1);
    plot(mean(cumprob_ids),'LineWidth',2,'Color','blue');
    hold on;
    plot(mean(cumprob_ads),'LineWidth',2,'Color','red');
    
    devi_ids = std(cumprob_ids)./sqrt(11);
    devi_ads = std(cumprob_ads)./sqrt(11);
    
    means_ids = mean(cumprob_ids);
    means_ads = mean(cumprob_ads);
    
    drawstds(h,1:target_length,mean(cumprob_ids),devi_ids,0.25,2,'blue');
    drawstds(h,1:target_length,mean(cumprob_ads),devi_ads,0.25,2,'red');
    grid;
    
    xlabel('position in the utterance');
    ylabel('F0 likelihood');
    ticknames = {};
    ticknames{1} = 'onset';
    ticknames{target_length} = 'offset';
    for k = 2:1:target_length-1
        ticknames{k} = sprintf('%d/%d',k,target_length);
    end
    set(gca,'XTick',1:1:target_length);
    set(gca,'XTickLabel',ticknames);
    legend({'IDS','ADS'},'Location','NorthEast');
    
    xlim([0.5 target_length+0.5])
        
    subplot(2,1,2)
        
    error_ids = std((cumf0_ids))./sqrt(11);
    info_ids = mean((cumf0_ids));
    
    error_ads = std((cumf0_ads))./sqrt(11);
    info_ads = mean((cumf0_ads));
    
    hold on;
    plot(info_ids,'LineWidth',2,'Color','blue');
    plot(info_ads,'LineWidth',2,'Color','red');
    grid;
    xlabel('position in the utterance');
    ylabel('z-score normalized log(F0)');
    ticknames = {};
    ticknames{1} = 'onset';
    ticknames{target_length} = 'offset';
    for k = 2:1:target_length-1
        ticknames{k} = sprintf('%d/%d',k,target_length);
    end
    set(gca,'XTick',1:1:target_length);
    set(gca,'XTickLabel',ticknames);
    
    xlim([0.5 target_length+0.5])
    
    errpatch_ads = shadedErrorBar([],info_ads,error_ads);
    
    errpatch_ads.patch.FaceColor = [1 0 0];
    errpatch_ads.patch.FaceAlpha = 0.3;
    errpatch_ads.mainLine.Color = [1 0 0];
    
    errpatch_ids = shadedErrorBar([],info_ids,error_ids);
    
    errpatch_ids.patch.FaceColor = [0 0 1];
    errpatch_ids.patch.FaceAlpha = 0.3;
    errpatch_ids.mainLine.Color = [0 0 1];
    
    plot(info_ids,'LineWidth',2,'Color','blue');
    plot(info_ads,'LineWidth',2,'Color','red');
    
    legend({'IDS','ADS'},'Location','SouthWest');
    ylim([-0.5 0.5])
  
    
    subplot(2,1,1);
    [h,p,ci,stat] = ttest(cumprob_ids,cumprob_ads,0.05/target_length);
    [siglevels,h] = holmBonferroni(p,0.05);
    mm = mean([mean(cumprob_ids);mean(cumprob_ads)]);
    for k = 1:length(h)        
        tmp = ylim;        
        %if(h(k)) % dont plot significances, just t-stats        
            %text(k,means_ads(k)+devi_ads(k),'*','HorizontalAlignment','center','FontSize',24);
        %end
        text(k,means_ads(k)+devi_ads(k)+0.017,sprintf('t = %0.2f\n',abs(stat.tstat(k))),'HorizontalAlignment','center','FontSize',14);
        
    end
    ylim([tmp(1) max(means_ads)+max(devi_ads)+0.03]);
    
    halfway = target_length/2;
        
    if(short)
        lengthname = 'short';
    else
        lengthname = 'long';
    end
    
    [r,p] = corr((1:length(means_ids))',means_ids','type','Pearson');
    fprintf('Correlation during full %s IDS: r = %0.4f (p = %0.4f).',lengthname,r,p);
    if(p < 0.05) fprintf(' (*)\n'); else fprintf('\n'); end
    
    [r,p] = corr((1:length(means_ids(1:halfway)))',means_ids(1:halfway)','type','Pearson');
    fprintf('Correlation during first 1/2 of %s IDS: r = %0.4f (p = %0.4f).',lengthname,r,p);
    if(p < 0.05) fprintf(' (*)\n'); else fprintf('\n'); end
    
    [r,p] = corr((1:length(means_ids(halfway+1:end)))',means_ids(halfway+1:end)','type','Pearson');
    fprintf('Correlation during last 1/2 of %s IDS: r = %0.4f (p = %0.4f).',lengthname,r,p);
    if(p < 0.05) fprintf(' (*)\n'); else fprintf('\n'); end
    
    [r,p] = corr((1:length(means_ads(1:end)))',means_ads(1:end)','type','Pearson');
    fprintf('Correlation during full %s ADS: r = %0.4f (p = %0.4f).',lengthname,r,p);
    if(p < 0.05) fprintf(' (*)\n'); else fprintf('\n'); end
    
    [r,p] = corr((1:length(means_ads(1:halfway)))',means_ads(1:halfway)','type','Pearson');
    fprintf('Correlation during first 1/2 of %s ADS: r = %0.4f (p = %0.4f).',lengthname,r,p);
    if(p < 0.05) fprintf(' (*)\n'); else fprintf('\n'); end
    
    [r,p] = corr((1:length(means_ads(halfway+1:end)))',means_ads(halfway+1:end)','type','Pearson');
    fprintf('Correlation during last 1/2 of %s ADS: r = %0.4f (p = %0.4f).',lengthname,r,p);
    if(p < 0.05) fprintf(' (*)\n'); else fprintf('\n'); end
    
    
    if(short)
        title('short')
    else
        title('long')
    end
    
end


%% -----------------------------------------------------------------------
% Part 2, measure correlation between utterance length and different
% descriptors
%-------------------------------------------------------------------------

% likelihoods

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

npos = target_length;
all_mean_short = zeros(length(F0prob),1);
all_mean_normal = zeros(length(F0prob),1);
all_mean_end = zeros(length(F0prob),1);
relpos = zeros(length(F0prob),npos);

nsyls = cellfun(@length,bounds_orig_syllable_t)-1;

for k = 1:length(F0prob)   
    relpos(k,:) = resample_by_interpolation(F0prob{k},npos);
end

figure(123);clf;
for k = 1:4
    
    if(k == 1)
        dat = all_mean(ids_i);
        dat_raw = all_mean_orig(ids_i);
        header = 'mean';
    elseif(k == 2)
        dat = all_SD(ids_i);
        dat_raw = all_SD_orig(ids_i);        
        header = 'SD';
    elseif(k == 3)
        dat = all_max(ids_i);
        dat_raw = all_max_orig(ids_i);        
        header = 'max';
    elseif(k == 4)
        dat = all_maxmin(ids_i);
        dat_raw = all_maxmin_orig(ids_i);               
        header = 'range';
    end    
    
    subplot(4,2,1+(k-1)*2);
    
    scatter(nsyls(ids_i),dat)
    
    [r,p] = corr(nsyls(ids_i),dat);
    a = polyfit(nsyls(ids_i),dat,1);
    
    x = min(nsyls(ids_i)):max(nsyls(ids_i));
    y = a(1).*x+a(2);
    hold on;
    plot(x,y,'Color','red','LineWidth',2)
    title(sprintf('F0 likelihood %s',header));
    ylabel(header);
    
    text(20,0,sprintf('r = %0.3f (p = %0.3f)',r,p));
    
    subplot(4,2,2+(k-1)*2);
            
    scatter(nsyls(ids_i),dat_raw)
    
    [r,p] = corr(nsyls(ids_i),dat_raw);
    a = polyfit(nsyls(ids_i),dat_raw,1);
    
    x = min(nsyls(ids_i)):max(nsyls(ids_i));
    y = a(1).*x+a(2);
    hold on;
    plot(x,y,'Color','red','LineWidth',2)
    title(sprintf('raw F0 %s',header));
    %title(sprintf('r = %0.3f (p = %0.3f)',r,p));
    text(20,0,sprintf('r = %0.3f (p = %0.3f)',r,p));
    
end

subplot(4,2,7);
xlabel('number of syllables');
subplot(4,2,8);
xlabel('number of syllables');


% Same thing for ads

valid_inds = intersect(ads_i,find(nsyls <= 100));

figure(124);clf;
for k = 1:4
    
    if(k == 1)
        dat = all_mean(valid_inds);
        dat_raw = all_mean_orig(valid_inds);
        header = 'mean';
    elseif(k == 2)
        dat = all_SD(valid_inds);
        dat_raw = all_SD_orig(valid_inds);        
        header = 'SD';
    elseif(k == 3)
        dat = all_max(valid_inds);
        dat_raw = all_max_orig(valid_inds);        
        header = 'max';
    elseif(k == 4)
        dat = all_maxmin(valid_inds);
        dat_raw = all_maxmin_orig(valid_inds);               
        header = 'range';
    end
    
    
    subplot(4,2,1+(k-1)*2);
    
    scatter(nsyls(valid_inds),dat)
    
    [r,p] = corr(nsyls(valid_inds),dat);
    a = polyfit(nsyls(valid_inds),dat,1);
    
    x = min(nsyls(valid_inds)):max(nsyls(valid_inds));
    y = a(1).*x+a(2);
    hold on;
    plot(x,y,'Color','red','LineWidth',2)
    title(sprintf('F0 likelihood %s',header));
    ylabel(header);
    
    text(20,0,sprintf('r = %0.3f (p = %0.3f)',r,p));
    
    subplot(4,2,2+(k-1)*2);
            
    scatter(nsyls(valid_inds),dat_raw)
    
    [r,p] = corr(nsyls(valid_inds),dat_raw);
    a = polyfit(nsyls(valid_inds),dat_raw,1);
    
    x = min(nsyls(valid_inds)):max(nsyls(valid_inds));
    y = a(1).*x+a(2);
    hold on;
    plot(x,y,'Color','red','LineWidth',2)
    title(sprintf('raw F0 %s',header));
    %title(sprintf('r = %0.3f (p = %0.3f)',r,p));
    text(20,0,sprintf('r = %0.3f (p = %0.3f)',r,p));
    
end

subplot(4,2,7);
xlabel('number of syllables');
subplot(4,2,8);
xlabel('number of syllables');


% Get linear correlation between predictabiliy at each position (IDS or
% ADS) and the number of syllables in the corresponding utterances.
for j = 1:npos
    r(j) = corr(nsyls(ids_i),relpos(ids_i,j));
    r2(j) = corr(nsyls(ads_i),relpos(ads_i,j));
end
figure;plot(r,'LineWidth',2)
hold on;
plot(r2,'LineWidth',2,'LineStyle','--');
set(gca,'XTickLabel',{'1/5th','2/5th','3/5th','4/5th','5/5th'});
set(gca,'XTick',1:5)
grid;
ylabel('correlation r');
xlabel('relative position in utterance');
legend({'IDS','ADS'},'Location','NorthWest');

set(gca,'XTick',1:1:npos);
set(gca,'XTickLabel',ticknames);







