function [M,lambda] = formatEM(seq,lags,maxval)
% function [M,lambda] = formatEM(seq,lags,maxval)
% Initializes a mixed-order Markov model from lagged bi-gram transition 
% frequencies for input "seq" and lags 1:m. 
%
% Optionally, size of the discrete alphabet (maxval) can be provided 
% if it differs from the maximum value of "seq".

if nargin <3
    maxval = max(seq);
end


m = length(lags);
alphsize = maxval;

lambda = ones(m,alphsize);
M = zeros(alphsize,alphsize,m);

% Estimate original values of M from frequencies

for t = m+1:length(seq)
    for k = 1:m
        M(seq(t-k),seq(t),k) = M(seq(t-k),seq(t),k)+1; 
    end
end

for k = 1:m
    for x = 1:alphsize
        M(x,:,k) = M(x,:,k)./sum(M(x,:,k));
    end
end

% Initialize lambda as 1/m

for k = 1:m-1
    lambda(k,:) = lambda(k,:)./sum(lambda(k,:));
end

