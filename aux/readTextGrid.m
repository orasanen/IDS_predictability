function [tga] = readTextGrid(filename)
%READTEXTGRID Reads TextGrid files based on two primary formats:
% [A]: Standard TextGrid tier structure (LONGFT)
% [B]: Table TextGrid structure delimited by commas (SHORTFT) - Note: a
% standard TextGrid file can be converted to Table format using Praat
%
% Input: filename of TextGrid file
% 
% Version 1.0
% Created on 01.04.2016
% SK
%

%
% Can also use a script for TextGrid conversions between different formats:
% (not available in this version)
% [status,result] = system('/Applications/Praat.app/Contents/MacOS/Praat --open "convert_to_table.praat"')
%

% Define file type
LONGFT = true;
SHORTFT = false;


% [A] Read long-text file type of TextGrid
if LONGFT
    
    % Read file into array
    fid = fopen(filename,'r');
    data = textscan(fid,'%s','Delimiter','\n');
    nd = length(data{1});
    
    % Preprocess data (if using double byte unicode character set)
    if round(length(data{1}{2})/2) == sum(double(data{1}{2})==0)
        for i = 1:nd
            data{1}{i} = data{1}{i}(2:2:end);
        end
    end
    
    % Process header information
    tga.start_t = str2double(regexp(data{1}{4},'\d+\.?\d*','match'));
    tga.end_t = str2double(regexp(data{1}{5},'\d+\.?\d*','match'));
    tga.num_tiers = str2double(regexp(data{1}{7},'\d+\.?\d*','match'));
    
    % Process file line by line
    i = 9;
    while i <= nd
       
        % Start new tier
        if strfind(data{1}{i},'name =')
            tier_name = cell2mat(regexp(data{1}{i}, '(?<=")[^"]+(?=")', 'match'));
            tga.(tier_name).start_tier_at = str2double(regexp(data{1}{i+1},'\d+\.?\d*','match'));
            tga.(tier_name).end_tier_at = str2double(regexp(data{1}{i+2},'\d+\.?\d*','match'));
            tga.(tier_name).tier_numel = str2double(regexp(data{1}{i+3},'\d+\.?\d*','match'));
            
            % Initialize vars
            tga.(tier_name).start_at = zeros(1,tga.(tier_name).tier_numel);
            tga.(tier_name).end_at = zeros(1,tga.(tier_name).tier_numel);
            tga.(tier_name).val = cell(1,tga.(tier_name).tier_numel);
            
            % Tier increment
            ti = 1;
            
            % Skip
            i = i+4;
        end
        
        % Get interval increment
        if strfind(data{1}{i},'intervals [')
            ti = str2double(regexp(data{1}{i},'\d+\.?\d*','match'));
        end
        
        % Add interval information
        if strfind(data{1}{i},'xmin =')
            tga.(tier_name).start_at(ti) = str2double(regexp(data{1}{i},'\d+\.?\d*','match'));
        elseif strfind(data{1}{i},'xmax =')
            tga.(tier_name).end_at(ti) = str2double(regexp(data{1}{i},'\d+\.?\d*','match'));
        elseif strfind(data{1}{i},'text =')
            tga.(tier_name).val{ti} = cell2mat(regexp(data{1}{i}, '(?<=")[^"]+(?=")', 'match'));
        end
        
        % Add code here
        
        % Increment
        i = i + 1;
        
    end
    
    fclose(fid);
    
end %[/A] 



% [B] Read short-text (Table format) file type of TextGrid
if SHORTFT
    
    % Read file into array
    fid = fopen(filename,'r');
    data = textscan(fid,'%s','Delimiter','\n');
    nd = length(data{1});
    
    % Preprocess data (if using double byte unicode character set)
    if round(length(data{1}{2})/2) == sum(double(data{1}{2})==0)
        for i = 1:nd
            data{1}{i} = data{1}{i}(2:2:end);
        end
    end
    
    % Iterate through data
    for i = 1:nd
        %procbar(i,nd);
        C(i,:) = strsplit(data{1}{i},',');
    end
    
    % Find unique annotation classes
    tier_names = unique(C(:,2));
    n_tiers = length(tier_names);
    
    % Iterate through tiers
    for i = 1:n_tiers
        
        inds = find(~cellfun(@(x)(isempty(x)),strfind(C(:,2),tier_names{i})));
        tga.(tier_names{i}).start_t = str2double(C(inds,1));
        tga.(tier_names{i}).end_t = str2double(C(inds,4));
        tga.(tier_names{i}).val = C(inds,3);
        
    end
    
    %tga.C = C;
    fclose(fid);
    
end %[/B] 


end % end of function