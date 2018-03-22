
clear all

filename1 = 'results/results_for_the_first_manuscript_submission/results_20-Oct-2017 00:04:48_ManyBabies_usesyllables1_framesize_100.mat'; % syllabic-frame
filename2 = 'results/results_for_the_first_manuscript_submission/results_19-Oct-2017 21:07:23_ManyBabies_usesyllables0_framesize_100.mat'; % fixed-frame

LL_per_frame = zeros(4,1);
LL_per_timeunit = zeros(4,1);
LL = zeros(4,1);

all_MOCM_syl = [];
all_LSTM_syl = [];
all_MOCM_uni = [];
all_LSTM_uni = [];

names = {'MOCM syllable','LSTM syllable','MOCM uniform','LSTM uniform'};

load(filename1)
durs = cellfun(@max,bounds_t)-cellfun(@min,bounds_t);
for k = 1:length(F0prob_MOMC)
    
    all_MOCM_syl = [all_MOCM_syl;F0prob_MOMC{k}];
    all_LSTM_syl = [all_LSTM_syl;F0prob_LSTM{k}];
    
    LL(1) = LL(1)+sum(log(F0prob_MOMC{k}));
    LL(2) = LL(2)+sum(log(F0prob_LSTM{k}));
    
   LL_per_frame(1) = LL_per_frame(1)+sum(F0prob_MOMC{k})./length(F0prob_MOMC{k})/length(F0prob_MOMC);
   LL_per_frame(2) = LL_per_frame(2)+sum(F0prob_LSTM{k})./length(F0prob_LSTM{k})/length(F0prob_LSTM);
      
   LL_per_timeunit(1) = LL_per_timeunit(1)+sum(F0prob_MOMC{k})./durs(k)/length(F0prob_MOMC);
   LL_per_timeunit(2) = LL_per_timeunit(2)+sum(F0prob_LSTM{k})./durs(k)/length(F0prob_LSTM);
   
end

load(filename2)
durs = cellfun(@max,bounds_t)-cellfun(@min,bounds_t);
for k = 1:length(F0prob_MOMC)
    
        
    all_MOCM_uni = [all_MOCM_uni;F0prob_MOMC{k}];
    all_LSTM_uni = [all_LSTM_uni;F0prob_LSTM{k}];
    
    LL(3) = LL(3)+sum(log(F0prob_MOMC{k}));
    LL(4) = LL(4)+sum(log(F0prob_LSTM{k}));
   LL_per_frame(3) = LL_per_frame(3)+sum(F0prob_MOMC{k})./length(F0prob_MOMC{k})/length(F0prob_MOMC);
   LL_per_frame(4) = LL_per_frame(4)+sum(F0prob_LSTM{k})./length(F0prob_LSTM{k})/length(F0prob_LSTM);
   
   LL_per_timeunit(3) = LL_per_timeunit(3)+sum(F0prob_MOMC{k})./durs(k)/length(F0prob_MOMC);
   LL_per_timeunit(4) = LL_per_timeunit(4)+sum(F0prob_LSTM{k})./durs(k)/length(F0prob_LSTM);
end

for k = 1:length(LL_per_frame)
   fprintf('%s LL per frame: %0.3f\n',names{k},LL_per_frame(k)); 
end


% Look at the correlations





  