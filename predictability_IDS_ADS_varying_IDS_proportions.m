% Scripts used to test intonation (F0) predictability on ManyBabies data 
% (The ManyBabies Consortium, 2017; https://osf.io/re95x/)
% as reported in :
%
%   Rasanen, O. Kakouros, S. & Soderstrom, M.: 
%   "Is infant-directed speech interesting because it is surprising? – 
%   Linking properties of IDS to statistical learning and attention at the 
%   prosodic level". Submitted for publication.
%
% Any use or derivation of this code should cite the above manuscript. 
%
% Note: requires YAAPT (Zahorian & Hu, 2008) for F0 estimation. 
% Version 4.0 was used for the results reported in the manuscript.
% Download from: http://ws2.binghamton.edu/zahorian/yaapt.htm
%
% Note 2: requires ManyBabies IDS/ADS stimuli
% Download the IDS and ADS samples from https://osf.io/xbv95/files/
% --> OSF Storage --> "IDS Clips Soderstrom" and "ADS Clips Soderstrom",
% and add the corersponding paths to the filedir1 and filedir2 variables below. 
%
% Note 3: the package includes oscillator-based syllabification algorithm
% published in Rasanen, Doyle & Frank (2018): Pre-linguistic
% segmentation of speech into syllable-like units. Cognition, 171, 130–150. 
% Direct download is also available at 
% https://github.com/orasanen/thetaOscillator
% 
% Note 4: the package includes MATLAB implementation of a mixed-order 
% Markov chain model (MOCM) (Saul & Pereira, 1997) by O. Rasanen. The
% scripts can be found in the subfolder /MixedOrderMarkov
% 
% Code by Okko Rasanen (okko.rasanen@aalto.fi). Last update: 22.3.2018
%-----
% THE MAIN DIFFERENCE TO predictability_IDS_ADS_main.m: 
% this script runs the main experiment five times by varying the
% proportion of IDS speech in the training data from 0 to 100% with 25%
% steps. This experiment was added to a revised version of the paper and is
% reported in section 5.2 of the manuscript.
% 
% To analyze the results, use idsads_analyze_distribution_proportions.m
% -----

clear all

% Parameters (defaults used in the manuscript main results):

polyorder = 2; % order of the polynomial fit to F0 tracks (1 or 2; default = 2)
Q_values = [6,12,24]; % codebook sizes used for F0 track quantization (default = [6, 12, 24])
total_folds = 10; % how many folds to use in training/testing (default = 10) 
n_repeats = 5; % how many times the experiment is repeated and averaged (default = 5)
mocm_order = 5; % order of the mixed-order Markov chain model (default = 5)
em_iters = 30; % how many EM iters to train mixed-order Markov chain models (default = 30) 
usesyllables = 1; % use syllables or fixed-frame windowing (default = 1)
framesize = 0.1; % segment length un fixed-frame windowing (seconds; default = 0.1) 
doLSTM = 0; % Do LSTM Analysis in addition to MOCM? (needs Python + Keras with Theano/Tensorflow Backend)  

pythondir = '/Users/orasanen/anaconda/bin/python'; % location of python executable (for LSTM)
filedir1 = '/Users/orasanen/speechdb/ManyBabies/ADS Clips Soderstrom/'; % location of ManyBabies ADS data
filedir2 = '/Users/orasanen/speechdb/ManyBabies/IDS Clips Soderstrom/'; % location of ManyBabies IDS data

dataset = 'ManyBabies'; 

% Get location of the current scripts and add paths

tmp = which('predictability_IDS_ADS_main.m');
[curdir,b,c] = fileparts(tmp);
addpath([curdir '/thetaOscillator']);
addpath([curdir '/thetaOscillator/gammatone']);
addpath([curdir '/MixedOrderMarkov']);
addpath([curdir '/aux']);
addpath([curdir '/result_analysis_and_plotting']);
 
% Get metadata and filenames for the ManyBabies recordings
% (add your own paths to the arguments)

[METADATA,filenames,METADATA_HEADERS] = loadManyBabiesData(filedir1,filedir2);

% Extract utterance F0 trajectories with YAAPT (if not already done before)
%
% Note: You need to download YAAPT from 
% http://ws2.binghamton.edu/zahorian/yaapt.htm or use some alternative
% pitch tracker. Original experiments conducted with YAAPT 4.0.0.

F0_filename = [curdir sprintf('/saves/F0_trajectories_and_metadata_%s.mat',dataset)];
if(~exist(F0_filename,'file'))
    F0_raw = cell(length(filenames),1);
    
    for k = 1:length(filenames)
        [x,fs] = audioread(filenames{k});
        
        x = resample(x,16000,fs);
        fs = 16000;
        
        ExtrPrm.f0_min = 120;   % Limit F0 min to 120 Hz
        ExtrPrm.f0_max = 600;   % Limit F0 max to 600 Hz
        [F0_raw{k}, numfrms, frmrate] = yaapt(x, fs,1,ExtrPrm,0,1);
        
        procbar(k,length(filenames));
    end
    
    save(F0_filename,'F0_raw','METADATA','filenames','METADATA_HEADERS');
else
    load(F0_filename,'F0_raw','METADATA','filenames','METADATA_HEADERS');
end

N = length(F0_raw);

F0_raw_orig = F0_raw;

% Do sonority-based syllabification with algorithm described in Rasanen,
% Doyle & Frank (in press). Included in the package
    
syllable_filename = [curdir sprintf('/saves/syllables_%s.mat',dataset)];
if(~exist(syllable_filename,'file'))
    [bounds,bounds_t] = thetaseg(filenames);
    save(syllable_filename,'bounds','bounds_t');
else
    load(syllable_filename,'bounds','bounds_t');
end

bounds_orig_syllable = bounds;
bounds_orig_syllable_t = bounds_t;

% Run YAAPT pitch fix and unvoiced interpolation using the following
% settings (edit to ptch_fix or change the command file)
%
% pitch_half = 2;
% pitch_half_sens = 3;
% pitch_double = 2;
% pitch_double_sens = 3.5;
% interp = 1;
% smooth_factor = 5;
% smooth = 0;
% extrap = 1;
% ptch_typ = 100;

for k = 1:N    
    a = find(F0_raw_orig{k} > 0);
    F0_raw{k} = ptch_fix(F0_raw{k});
    F0_raw{k} = medfilt1(F0_raw{k},3); % small medfilt to remove outliers
    F0_raw_orig{k}(a) = F0_raw{k}(a); % update to non-interpolated tracks
end

% Take logarithm and z-score normalize all utterance-level F0 contours

for k = 1:N
    F0_raw{k} = log(F0_raw{k});
    F0_raw{k} = F0_raw{k}-mean(F0_raw{k});
    F0_raw{k} = F0_raw{k}./std(F0_raw{k});         
end
  
% Divide syllables into N sub-segments for modeling (default N = 2)

slices_per_syllable = 2;

for k = 1:length(bounds)
    bb = bounds{k};
    
    all_nb = [];
    for jj = 1:length(bb)-1 % go through all bounds
        
        b1 = bb(jj);
        b2 = bb(jj+1);
        
        nb = [];
        nb(1) = b1;
        y = b1:b2;
        
        len = floor(length(y)/(slices_per_syllable));
        
        for r = 1:slices_per_syllable-1
            nb(r+1) = nb(r)+len;
        end
        nb(end+1) = b2;
        all_nb = [all_nb;nb'];
    end
    all_nb(slices_per_syllable+1:slices_per_syllable+1:end-slices_per_syllable+1) = [];
    
    bounds{k} = all_nb;
    bounds_t{k} = bounds{k}./1000;
end

% Use uniformly spaced segments in addition to the syllables

bounds_t_uni = cell(size(bounds_t));
bounds_uni = cell(size(bounds_t));

for k = 1:length(bounds_t)
    
    beg = bounds_t{k}(1);
    ending = bounds_t{k}(end);
    
    bb = (beg:framesize:ending)';
    
    bounds_t_uni{k} = bb;
    bounds_uni{k} = round(bb*1000);
end

% Use segmentation method determined in hyperparams
if(~usesyllables)
   bounds_t = bounds_t_uni;
   bounds = bounds_uni;
end
    

% Remove segments that precede or follow the utterance (as defined by
% voicing).
for k = 1:N
  
    utt_onset = find(F0_raw{k} ~= F0_raw{k}(1),1);
    utt_offset = length(F0_raw{k})-find(fliplr(F0_raw{k}) ~= F0_raw{k}(1),1)-1;
    
    tmp1 = find(bounds_t{k} < utt_onset/100);
    if(isempty(tmp1))
        tmp1 = 1;
    end
        
    tmp2 = find(bounds_t{k} > utt_offset/100);
    if(isempty(tmp2))
        tmp2 = length(bounds_t{k});
    end
    
    bounds_t{k} = bounds_t{k}(tmp1(end):tmp2(1));
    bounds{k} = bounds{k}(tmp1(end):tmp2(1));
        
    bo = bounds_orig_syllable_t{k};
    bo(bo < bounds_t{k}(1)) = [];
    bo(bo > bounds_t{k}(end)) = [];
    
    % Fix original syllable boundaries as well
    
    bounds_orig_syllable{k} = round(bo*100);
    bounds_orig_syllable_t{k} = bo;            
end
 

% Parametrize F0 during each segment using polyfit
% Also store parametrized & reconstructed and dc-corrected F0 trajectories 
% for later visualization.

F0PARAMS = cell(N,1); % polyfit parameters for each utterance
CONCAT = cell(N,1); % reconstructed F0 contours

for k = 1:N
    F0PARAMS{k} = zeros(length(bounds{k})-1,polyorder);
    CONCAT{k} = zeros(size(F0_raw{k}));
        
    counter = 1;
    for j = 2:length(bounds{k})
        b1 = max(1,floor(bounds_t{k}(j-1)*100));
        b2 = ceil(bounds_t{k}(j)*100);
        
        fy = F0_raw{k}(b1:min(b2,length(F0_raw{k})));
        
        %if(sum(fy == mode(F0_raw{k})) == length(fy))
        %    disp('something went wrong');
           % figure(6);plot(F0_raw{k});
           % pause;
        %end
        
        if(length(fy) <= polyorder)
            error('segment is shorter or equal to order of polyfit');
        end
        
        a = polyfit((1:length(fy))',fy',polyorder);
        
        % Pitch changes in this segment
        
        F0PARAMS{k}(counter,:) = a(1:polyorder);
        counter = counter+1;
                
        % Create also concatenation of the syllable-wise F0 fits for
        % visualization purposes.
                        
        x = (1:(b2-b1+1))';
        if(polyorder == 2)
            y = a(1).*x.^2+a(2).*x+a(3);
        else
            y = a(1).*x+a(2);
        end
        CONCAT{k}(b1:b2) = y;
    end
    F0PARAMS{k} = F0PARAMS{k}(1:counter-1,:);
end

% Remove utterances with less than 3 valid segments

segments_per_utterance = cellfun(@length,F0PARAMS);
a = find(segments_per_utterance < 3);

F0_raw(a) = [];
F0_raw_orig(a) = [];
F0PARAMS(a) = [];
bounds(a) = [];
bounds_orig_syllable(a) = [];
bounds_orig_syllable_t(a) = [];
bounds_t(a) = [];
METADATA(a,:) = [];
filenames(a) = [];
CONCAT(a) = [];

% Get duration of each utterance (from first and last syllable boundaries)
durs = cellfun(@max,bounds_t)-cellfun(@min,bounds_t);

% Pointers to IDS and ADS utterances
ids_i = cellfind(METADATA(:,1),'IDS');
ads_i = cellfind(METADATA(:,1),'ADS');


N = length(F0_raw); % Number of utterances
 
padlength = max(10,mocm_order);  % constant zero padding at utterance onsets

id_proportions = [0:0.25:1];

F0prob_MOMC = cell(N,length(id_proportions));
F0prob_LSTM = cell(N,length(id_proportions));

% Pre-allocate memory for LSTM likelihood contours since they have to be
% reconstructed frame-by-frame.
for k = 1:length(F0prob_LSTM)
   F0prob_LSTM{k} = zeros(length(F0PARAMS{k})+padlength,1); 
end

totrain = floor((1-1/total_folds)*N); % How many training signals in each fold
counts = zeros(N,length(id_proportions)); 

% Run the experiment n_repeats times

for iter = 1:n_repeats
    
    % Define utterances that are allowed for training (all IDS + enough ADS
    % to match IDS in length).
    ads_order = ads_i(randperm(length(ads_i)));
    ads_cumdur = cumsum(durs(ads_order));
    ids_order = ids_i(randperm(length(ids_i)));
    ids_cumdur = cumsum(durs(ids_order));
            
    % Get performance for different proportions of IDS in the training data
    
    total_dur_ids = sum(durs(ids_i));
    total_dur_ads = sum(durs(ads_i));
        
    for propiter = 1:length(id_proportions)
        
        idprop = id_proportions(propiter);        
        maxdur = sum(durs(ids_i)); % maximum duration of data
        
        ads_to_include = [];
        ids_to_include = [];
        
        
        % Get the desired amount of IDS and ADS speech (in seconds).
        
        a = find(cumsum(durs(ids_order)) >= maxdur.*idprop,1);
        if(isempty(a))
            a = length(ids_order);
        elseif(idprop == 0)
            a = 0;
        end
        
        ids_to_include = ids_order(1:a);
        
        b = find(cumsum(durs(ads_order)) >= maxdur.*(1-idprop),1);
        if(idprop == 1)
            b = 0;
        end
        
        ads_to_include = ads_order(1:b);
        
        valid_traininds = sort([ids_to_include;ads_to_include],'ascend');
        
        fprintf('Total dur IDS: %0.2f min. Total durs ADS: %0.2f min. Total dur: %0.2f min.\n',sum(durs(ids_to_include))/60,sum(durs(ads_to_include))/60,sum(durs(valid_traininds))/60);
                        
        % Randomize the order of the utterances for each repeat of the
        % experiment
        rng(123+iter,'twister');
        
        allinds = (randperm(N))';
        
        % Do N-fold evaluation
        
        for fold = 1:total_folds
            fprintf('Experiment %d out of %d. Fold %d out of %d.\n',iter,n_repeats,fold,total_folds);
            
            % Select training and testing samples
            % (rotation of allinds-vector is at the end of the n-fold loop)
            traininds = allinds(1:totrain);
            testinds = allinds(totrain+1:end);
                                    
            % Combine all training F0 into one long signal
            % (result does not depend on the order due to the padding).
            allf0_train = zeros(sum(cellfun(@length,F0PARAMS(traininds))),size(F0PARAMS{1},2));
            
            wloc = 1;
            for k = traininds'
                len = length(F0PARAMS{k});
                % Check that utterance is valid (IDS or within the valid ADS
                % pool matched in duration)
                if(length(intersect(k,valid_traininds)))
                    allf0_train(wloc:wloc+len-1,:) = F0PARAMS{k};
                    
                    allf0_train(wloc+len:wloc+len+padlength-1) = 0;
                    wloc = wloc+len+padlength;
                end
            end
            allf0_train = allf0_train(1:wloc-padlength-1,:);
            
            % Quantize F0 parameters
            
            qiter = 1;
            for Q = Q_values
                fprintf('Codebook size %d.\n',Q);
                % Get VQ codebook for training data
                [data_q,centroids] = kmeans(allf0_train,Q);
                
                % Get VQ sequences for all utterances
                F0_Q = cell(size(F0PARAMS));
                for k = 1:length(F0PARAMS)
                    if(~isempty(F0PARAMS{k}))
                        if(size(F0PARAMS{k},2) > 2)
                            F0PARAMS{k} = F0PARAMS{k}';
                        end
                        % add padding
                        tmp = [zeros(padlength,size(F0PARAMS{k},2));F0PARAMS{k}];
                        % get nearest centroid with euclidean distance
                        d = pdist2(tmp,centroids);
                        [~,F0_Q{k}] = min(d,[],2);
                    end
                end
                
                % Train mixed-order Markov chain model on the training data
                [M,lambda] = getMixedOrderModel(data_q,mocm_order,em_iters);
                
                % Get F0 likelihoods on the test utterances using MOMC
                for k = testinds'
                    if(~isempty(F0prob_MOMC{k,propiter}))
                        F0prob_MOMC{k,propiter} = F0prob_MOMC{k,propiter}+getSeqLikelihood(F0_Q{k},M,lambda);
                        counts(k,propiter) = counts(k,propiter)+1;
                    else
                        F0prob_MOMC{k,propiter} = getSeqLikelihood(F0_Q{k},M,lambda);
                        counts(k,propiter) = counts(k,propiter)+1;
                    end
                end
                
                % Gather data into fixed-length sequences for LSTM training.
                
                if(doLSTM)
                    
                    seqlen = padlength;
                    
                    train_in = zeros(999999,seqlen); % discrete sequence
                    train_out = zeros(999999,Q); % distributed softmax target
                    y = 1;
                    for k = 1:length(traininds)
                        seq = F0_Q{traininds(k)};
                        for j = 1:length(seq)-padlength
                            train_in(y,:) = seq(j:j+padlength-1);
                            train_out(y,:) = zeros(1,Q);
                            train_out(y,seq(j+padlength)) = 1;
                            y = y+1;
                        end
                    end
                    
                    train_in = train_in(1:y-1,:); % discrete sequence
                    train_out = train_out(1:y-1,:); % distributed softmax target
                    
                    test_in = zeros(999999,seqlen);
                    test_out = zeros(999999,Q);
                    test_utt_ind = zeros(999999,2);
                    y = 1;
                    for k = 1:length(testinds)
                        seq = F0_Q{testinds(k)};
                        for j = 1:length(seq)-padlength
                            test_in(y,:) = seq(j:j+padlength-1);
                            test_out(y,:) = zeros(1,Q);
                            test_out(y,seq(j+padlength)) = 1;
                            test_utt_ind(y,1) = testinds(k);
                            test_utt_ind(y,2) = j;
                            y = y+1;
                        end
                    end
                    
                    test_in = test_in(1:y-1,:);
                    test_out = test_out(1:y-1,:);
                    test_utt_ind = test_utt_ind(1:y-1,:);
                    
                    % Run LSTM in Python with Keras DNN library (tested with Theano
                    % backend).
                    
                    save saves/lstm_traindata.mat train_in train_out test_in test_out
                    
                    % Format the python command  to run with system().
                    % (replace  first part with your own path to python or add to system path)
                    ss = [pythondir ' ' curdir '/idsads_LSTM.py' ' ' curdir '/saves/'];
                    
                    [status,result] = system(ss);
                    
                    % load LSTM likelihood outputs
                    load saves/preds_out preds
                    
                    % Put LSTM results back together into original segmental sequences
                    for k = 1:size(preds,1)
                        signal = test_utt_ind(k,1);
                        position = test_utt_ind(k,2)+padlength;
                        trueval = find(test_out(k,:) == 1);
                        prob = preds(k,trueval);
                        F0prob_LSTM{signal,propiter}(position) = F0prob_LSTM{signal,propiter}(position)+prob;
                        
                    end
                end
                qiter = qiter+1;
            end
            % Rotate training/testing sets between the folds
            allinds = circshift(allinds,length(testinds));
        end
    end
end

% Remove zero paddings and get average likelihoods by dividing with the 
% number of times the signal was measured in the test set

for propiter = 1:length(id_proportions)
    for k = 1:N
        if(~isempty(F0prob_MOMC{k,propiter}))
            F0prob_MOMC{k,propiter} = F0prob_MOMC{k,propiter}(padlength+1:end);
            F0prob_MOMC{k,propiter} = F0prob_MOMC{k,propiter}./counts(k,propiter);
            F0prob_LSTM{k,propiter} = F0prob_LSTM{k,propiter}(padlength+1:end);
            F0prob_LSTM{k,propiter} = F0prob_LSTM{k,propiter}./counts(k,propiter);
        end
    end
end

% Save results in a file

result_filename = sprintf('%s/results/results_%s_%s_usesyllables%d_framesize_%d_proportions.mat',curdir,datestr(now),dataset,usesyllables,framesize*1000);
save(result_filename,'F0prob_MOMC','F0prob_LSTM','F0_raw_orig','METADATA','F0_raw','F0PARAMS','bounds','bounds_t','bounds_orig_syllable','bounds_orig_syllable_t','filenames','N')

