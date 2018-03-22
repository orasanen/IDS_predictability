function [METADATA,filenames,METADATA_HEADERS] = loadManyBabiesData(filedir_ADS,filedir_IDS)
%function [METADATA,filenames,METADATA_HEADERS] = loadManyBabiesData(filedir_ADS,filedir_IDS)
%
% Get signal filenames and related metadata for ManyBabies corpus
%
% Uses IDS Clips Soderstrom and ADS Clips Soderstrom from https://osf.io/re95x/

% Add correct path here
if nargin <1
    filedir_ADS = '/Users/orasanen/speechdb/ManyBabies/ADS Clips Soderstrom/';
end
if nargin <2
    filedir_IDS = '/Users/orasanen/speechdb/ManyBabies/IDS Clips Soderstrom/';
end

if(~exist(filedir_ADS,'dir'))
    error('Directory: %s does not exist. Must specify location of ManyBabies / ADS Clips Soderstrom',filedir_ADS);
end

if(~exist(filedir_IDS,'dir'))
    error('Directory: %s does not exist. Must specify location of ManyBabies / IDS Clips Soderstrom',filedir_IDS);
end

a = dir(filedir_ADS);
a = removenondir(a); % Remove non-directory files (OSX indexing etc.)

b = dir(filedir_IDS);
b = removenondir(b); % Remove non-directory files (OSX indexing etc.)

totfiles = 0;
% Go through each subdirectory to get number of files
for rootdir = 1:length(a)
    subdir = [filedir_ADS a(rootdir).name '/'];        
    c = dir([subdir '*.wav']);    
    totfiles = totfiles+length(c);        
end

for rootdir = 1:length(b)
    subdir = [filedir_IDS b(rootdir).name '/'];        
    c = dir([subdir '*.wav']);    
    totfiles = totfiles+length(c);        
end


METADATA = cell(totfiles,5);
filenames = cell(totfiles,1);
%% Get actual filenames and stuff

signal = 1;
for ss = 1:2
    for rootdir = 1:length(a)
        
        if(ss == 1)
            subdir = [filedir_ADS a(rootdir).name '/'];
        else
            subdir = [filedir_IDS b(rootdir).name '/'];
        end
        c = dir([subdir '*.wav']);
        
        for fileind = 1:length(c) % go through wavs
            filename = c(fileind).name;
            
            tmp = strfind(filename,',');
            if(tmp) % labeling event?
                label = filename(1:tmp-1);
            else
                label = 'N/A';
            end
            
            tmp2 = strfind(filename,'_');
            if(isempty(tmp2))
               tmp2 = strfind(filename,' at'); 
               atcase = 1;
            else
                atcase = 0;
            end
            BabyID = filename(tmp2-5:tmp2-1);
            
            
            tmp3 = strfind(filename,'.wav');
            if(~atcase)
                onset_time = filename(tmp2+1:tmp3-1);
            else
                onset_time = filename(tmp2+4:tmp3-1);
            end
            
            if(ss == 1)
                style = 'ADS';
            else
                style = 'IDS';
            end
                        
            if(~isempty(strfind(lower(subdir),'unfamiliar')))
                item_type = 'unfamiliar';
            elseif(~isempty(strfind(lower(subdir),'familiar')))
                item_type = 'familiar';
            else
                item_type = 'N/A';
            end
                
            METADATA{signal,1} = style;
            METADATA{signal,2} = BabyID;
            METADATA{signal,3} = label;
            METADATA{signal,4} = item_type;
            METADATA{signal,5} = onset_time;
            
            filenames{signal} = [subdir filename];
            signal = signal+1;
        end
    end
end

METADATA_HEADERS = {'style','BabyID','object label','familiar?','onset_time'};

METADATA = strtrim(METADATA); % remove potential leading and trailing whitespaces