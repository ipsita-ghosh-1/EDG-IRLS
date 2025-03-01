function Dist = kappainv_EDG(Gram)
%kappainv_EDG Summary of this function goes here
%   Detailed explanation goes here
n = size(Gram,1);
Dist = 1.*(ones(n,1)*diag(Gram')'+ ...
    diag(Gram)*ones(1,n)-Gram-Gram);
end