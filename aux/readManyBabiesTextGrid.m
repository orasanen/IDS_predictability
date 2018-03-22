function anno = readManyBabiesTextGrid(filename)

fid = fopen(filename,'r');
data = textscan(fid,'%s','Delimiter','\n');
fclose(fid);

data = data{1};

k = 1;
while(~strcmp(data{k},'intervals [2]:'))
    k = k+1;
end
tmp = data{k+1};
utt_onset_time = str2num(tmp(8:end));
tmp = data{k+2};
utt_offset_time = str2num(tmp(8:end));
tmp = data{k+3};
utt_ortho = tmp(9:end-2);

i = find(strcmp(data,'text = "kw" '));

kw_onset_time = zeros(length(i),1);
kw_offset_time = zeros(length(i),1);
for j = 1:length(i)
    tmp = data{i(j)-2};
    kw_onset_time(j) = str2num(tmp(8:end));
    tmp = data{i(j)-1};
    kw_offset_time(j) = str2num(tmp(8:end));
end

anno = struct;
anno.utt_onset_time = utt_onset_time;
anno.utt_offset_time = utt_offset_time;
anno.kw_onset_time = kw_onset_time;
anno.kw_offset_time = kw_offset_time;