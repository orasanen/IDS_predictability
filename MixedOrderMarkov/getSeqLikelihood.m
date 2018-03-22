function LL = getSeqLikelihood(seq,M,lambda)
% function LL = getSeqLikelihood(seq,M,lambda)
%
% Calculates the element-by-element likelihoods of input sequence "seq" 
% given the existing mixed-order Markov model specified by "M" and "lambda". 

m = size(lambda,1);
LL = zeros(length(seq),1);

for t = m+1:length(seq)

    LL(t) = getP(seq(1:t),M,lambda);
 
end