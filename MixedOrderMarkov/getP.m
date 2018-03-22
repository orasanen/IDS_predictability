function [P,deno] = getP(seq,M,lambda)
% function [P,deno] = getP(seq,M,lambda)
%
% Computes probability P(w[t]|w[t-1],w[t-2],...,w[t-m]) using an existing
% mixed-order Markov model specified by M and lambda.

m = size(lambda,1);

P = 0;

for k = 1:m
    if(k > 1)
        deno = 1-lambda(1,seq(end-1));
        for j = 2:k-1
            deno = deno*(1-lambda(j,seq(end-j)));
        end
    else
        deno = 1;
    end
    P = P+lambda(k,seq(end-k))*M(seq(end-k),seq(end),k)*deno;
end



