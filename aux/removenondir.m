function a = removenondir(a)

% Remove non-directory files (OSX indexing etc.)
toremove = [];
for k = 1:length(a)
    if(~a(k).isdir || strcmp(a(k).name,'.') || strcmp(a(k).name,'..'))
        toremove = [toremove k];    
    end
end
a(toremove) = [];
