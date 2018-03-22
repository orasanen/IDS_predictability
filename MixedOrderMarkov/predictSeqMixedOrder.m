function [prob,prediction,P] = predictSeqMixedOrder(seq,M,lambda)
% [prob,prediction] = predictSeqMixedOrder(seq,M,lambda)
% 
% Predicts the most likely element following seq, given the existing
% mixed-order Markov model specified in M and lambda. Also outputs 
% the predictive distribution P.

P = zeros(size(M,1),1);

for alph = 1:size(M,1)

    seqn = [seq;alph];    
    [P(alph),deno] = getP(seqn,M,lambda);
end

[prob,prediction] = max(P);