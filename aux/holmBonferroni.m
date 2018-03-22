function [siglevels,h,adjusted_p] = holmBonferroni(p,alpha)
% function [siglevels,h] = holmBonferroni(p,alpha)

if nargin <2
    alpha = 0.05;
end

if(size(p,1) > 1 && size(p,2) > 1)
    error('p must be a vector');
end
if(size(p,2) > 1)
    p = p';
end

[psort,i] = sort(p,'ascend'); 

m = length(p);

siglev = zeros(size(p));
for k = 1:length(i)    
   siglev(k) =  alpha/(m+1-k);
end

adjusted_p = zeros(size(p));

for k = 1:length(i)
    j = (1:k)';
    adjusted_p(k) = min(1,max((m-j+1).*psort(j)));
end


h = zeros(size(p));
a = find(psort > siglev,1);
if(~isempty(a))
    h(i(1:a-1)) = 1;
else
    h(:) = 1;
end

siglevels = zeros(length(siglev),1);
adjusted_ps = zeros(length(siglev),1);
for k = 1:length(siglev)
    siglevels(i(k)) = siglev(k);
    adjusted_ps(i(k)) = adjusted_p(k);    
end

adjusted_p = adjusted_ps;




