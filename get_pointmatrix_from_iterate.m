function P_c = get_pointmatrix_from_iterate(X_c,r)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[UG1,Lambda] = eig(X_c.Gam1(1:r,1:r));
P_c = X_c.U(:,1:r)*UG1*diag(sqrt(max(diag(Lambda),0)));
end