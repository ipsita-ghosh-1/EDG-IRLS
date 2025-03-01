function Gram = kappa_EDG(Dist)
%kappa_EDG Given a matrix of pairwise squared distances 
% 'Dist', computes the associated gram matrix 'Gram'.
n = size(Dist,1);
Gram = -0.5.*(eye(n)-ones(n,n)./n)...
    *Dist*(eye(n)-ones(n,n)./n);
end