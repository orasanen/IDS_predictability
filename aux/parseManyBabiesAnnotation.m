function [anno,annofilenames] = parseManyBabiesAnnotation(annodir)

if nargin <1
    annodir = '/Users/orasanen/speechdb/ManyBabies/manual_annotations/';
end



% Add correct path here

filedir_ADS = [annodir '/ADS Clips Soderstrom/'];

filedir_IDS = [annodir '/IDS Clips Soderstrom/'];


a = dir(filedir_ADS);
a = removenondir(a); % Remove non-directory files (OSX indexing etc.)

b = dir(filedir_IDS);
b = removenondir(b); % Remove non-directory files (OSX indexing etc.)

totfiles = 0;
% Go through each subdirectory to get number of files
for rootdir = 1:length(a)
    subdir = [filedir_ADS a(rootdir).name '/'];        
    c = dir([subdir '*.TextGrid']);    
    totfiles = totfiles+length(c);        
end

for rootdir = 1:length(b)
    subdir = [filedir_IDS b(rootdir).name '/'];        
    c = dir([subdir '*.TextGrid']);    
    totfiles = totfiles+length(c);        
end


%METADATA = cell(totfiles,5);
%filenames = cell(totfiles,1);
%% Get actual filenames and stuff
anno = {};
annofilenames = {};

signal = 1;
for ss = 2:2
    for rootdir = 1:length(b)
        
        if(ss == 1)
            subdir = [filedir_ADS a(rootdir).name '/'];
        else
            subdir = [filedir_IDS b(rootdir).name '/'];
        end
        c = dir([subdir '*.TextGrid']);
        
        for fileind = 1:length(c) % go through TextGrids
            filename = c(fileind).name;
            
            anno{signal} = readManyBabiesTextGrid([subdir filename]);
            anno{signal}.filename = filename;
            annofilenames{signal} = filename;
                        
            signal = signal+1;
        end
    end
end

anno = anno';
annofilenames = annofilenames';

for k = 1:length(annofilenames)
    s = annofilenames{k};
    tmp = strfind(s,'__');
    s(tmp:tmp+1) = ', ';
    tmp2 = strfind(s,'_');
    s(tmp2(1)) = ' ';
    s(tmp2(end)) = '.';
    if(length(tmp2) == 4)
        s(tmp2(2)) = ' ';
        s(tmp2(3)) = ' ';
    end
    tmp3 = strfind(s,'.TextGrid');
    s(tmp3:end) = [];
    annofilenames{k} = s;    
end
