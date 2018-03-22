function anno_out = mergeAnno(anno1,anno2)

assert(length(anno1) == length(anno2),'Both annotations must contain the same number of signals')




anno_out = anno1; % Use anno1 as reference and move fields from anno2 to anno1

for k = 1:length(anno1)
    a1 = fieldnames(anno1{k});
    a2 = fieldnames(anno2{k});
   
    len = min(length(anno1{k}.filename),length(anno2{k}.filename));
    assert(strcmp(anno1{k}.filename(1:len),anno2{k}.filename(1:len)),'signal filenames are not matching')
    
    for j = 1:length(a2)
       if(~isfield(anno1{k},a2{j}))
          anno_out{k}.(a2{j}) = anno2{k}.(a2{j}); 
       end
    end     
end








