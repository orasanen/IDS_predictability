function anno_converted = convertTextGridAnno(anno)
% function anno_converted = convertTextGridAnno(anno)
%
% Get word info from TextGrid tiers into more usable format 

anno_converted = cell(length(anno),1);

for signal = 1:length(anno)
    
    
    dat = anno{signal};
    a = fieldnames(dat);
    
    % Get structs within struct to extract TextGrid tiers from annotation
    tmp2 = zeros(length(a),1);
    for k = 1:length(a)
        tmp = getfield(dat,a{k});
        tmp2(k) = isstruct(tmp);
    end
    
    % Get utterances and timestamps for each tier
    
    counter_word = 1;
    counter_utterance = 1;
    talker_id = cell(0,0);
    t_onset = zeros(0,0);
    t_offset = zeros(0,0);
    word = cell(0,0);
    
    t_onset_utterance = zeros(0,0);
    t_offset_utterance = zeros(0,0);
    utterance = cell(0,0);
    talker_id_utterance = cell(0,0);
    
    tiernames = a(tmp2 == 1);
    for k = 1:length(tiernames)
        curtier = dat.(tiernames{k});
        
        utterances = curtier.val;
        onsets = curtier.start_at;
        offsets = curtier.end_at;
        
        for j = 1:length(utterances)
            if(~isempty(utterances{j}))
                s = utterances{j};
                s_split = strsplit(s);
                n_words = length(s_split);
                utt_dur = offsets(j)-onsets(j);
                utterance{counter_utterance} = s;
                t_onset_utterance(counter_utterance) = onsets(j);
                t_offset_utterance(counter_utterance) = offsets(j);
                talker_id_utterance{counter_utterance} = tiernames{k};
                counter_utterance = counter_utterance+1;
                
                word_dur = utt_dur/n_words;
                for n = 1:n_words
                    t_onset(counter_word) = onsets(j)+(n-1)*word_dur;
                    t_offset(counter_word) = onsets(j)+(n)*word_dur;
                    word{counter_word} = s_split{n};
                    talker_id{counter_word} = tiernames{k};
                    counter_word = counter_word+1;
                end
            end
        end
    end
    
    % Sort into temporal order
    [t_onset,i] = sort(t_onset,'ascend');
    t_offset = t_offset(i);
    word = word(i);
    talker_id = talker_id(i);
    
    [t_onset_utterance,i] = sort(t_onset_utterance,'ascend');
    t_offset_utterance = t_offset_utterance(i);
    utterance = utterance(i);
    talker_id_utterance = talker_id_utterance(i);
    
    anno_converted{signal}.word = word;
    anno_converted{signal}.t_onset_word = t_onset;
    anno_converted{signal}.t_offset_word = t_offset;
    anno_converted{signal}.talker_id = talker_id;
    
    anno_converted{signal}.utterance = utterance;
    anno_converted{signal}.t_onset_utterance = t_onset_utterance;
    anno_converted{signal}.t_offset_utterance = t_offset_utterance;
    anno_converted{signal}.talker_id_utterance = talker_id_utterance;
    
    anno_converted{signal}.filename = dat.filename;
    
end