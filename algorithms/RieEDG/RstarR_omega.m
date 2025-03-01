function RstarR_X = RstarR_omega(X,samples)
n=size(X,1);
D_in = ones(n,1)*diag(X)'+diag(X)*ones(1,n)-2*X;
D_hat = P_omega(D_in,samples);  
JXJ = -1/2*(D_hat - 1/n*ones(n,1)*sum(D_hat,2)' - 1/n*sum(D_hat,1)'*ones(1,n) + 1/n^2*sum(sum(D_hat))*ones(n));
P_JXJ = P_omega(JXJ,samples);
RstarR_X = P_JXJ - diag(sum(P_JXJ,1));

return